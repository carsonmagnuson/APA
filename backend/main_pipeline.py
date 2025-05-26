import logging
import os
from typing import Dict, Any, Optional

# Assuming these custom modules are in the same directory or PYTHONPATH
from data_fetching import pull_model_organism_orthologs, pull_OG_fasta
from sequence_processing import (
    fasta_to_seqrecord,
    get_model_organism_genes,
    select_genes_from_seqrecord,
    select_big_taxa_seq,
    select_biggest_k_seq,
    pad_sequence_lengths,
    remove_stop_codons_in_multiple,
    convert_organism_id_to_names, # Added this step based on original app.py
    format_ids_and_create_alignment
)
from phylogenetic_analysis import construct_phylo_tree, run_codeml_positive_selection
from file_utils import write_phylip_manual, write_newick_tree_with_header, CACHE_DIR

# --- Configure Logging ---
# This basic configuration is for when the script is run directly.
# If this module is imported, the calling application should configure logging.
if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(name)s - %(module)s.%(funcName)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

logger = logging.getLogger(__name__)

# --- Constants for default file names and directories ---
DEFAULT_OUTPUT_PHYLIP_FILE = 'output.phylip'
DEFAULT_OUTPUT_NEWICK_FILE = 'tree.nwk'
DEFAULT_PAML_OUTPUT_DIR = "paml_output" # Relative to where PAML runs

def run_phylogenetic_pipeline(
    gene_query: str,
    target_gene_orthodb_id: str,
    ortholog_group_id: str,
    k_sequences: int = 10,
    output_phylip_file: str = DEFAULT_OUTPUT_PHYLIP_FILE,
    output_newick_file: str = DEFAULT_OUTPUT_NEWICK_FILE,
    paml_output_dir: str = DEFAULT_PAML_OUTPUT_DIR,
    cleanup_paml_dir: bool = True,
    rooted_tree: bool = False
) -> Optional[Dict[str, Any]]:
    """
    Executes the full phylogenetic analysis pipeline.

    This pipeline fetches ortholog data, processes sequences, constructs a
    phylogenetic tree, and runs CODEML for positive selection analysis.

    Args:
        gene_query: The gene identifier or name to search for (e.g., "WBGene00004963").
        target_gene_orthodb_id: The OrthoDB ID of the primary gene of interest
                                (e.g., "6239_0:000672"). This gene's sequence
                                will be included in the analysis.
        ortholog_group_id: The OrthoDB ortholog group ID to fetch sequences from
                           (e.g., "430340at2759").
        k_sequences: The total number of sequences desired for the analysis,
                     including the target gene. Defaults to 10.
        output_phylip_file: Path to save the generated PHYLIP alignment file.
        output_newick_file: Path to save the generated Newick tree file.
        paml_output_dir: Directory for CODEML output files.
        cleanup_paml_dir: Whether to remove the PAML output directory after analysis.
        rooted_tree: Whether to construct a rooted (UPGMA) or unrooted (NJ) tree.
                     Defaults to False (unrooted NJ tree).

    Returns:
        A dictionary containing the parsed CODEML results if successful,
        otherwise None.
    """
    logger.info(f"Starting phylogenetic analysis pipeline for gene query: {gene_query}, target ID: {target_gene_orthodb_id}, OG: {ortholog_group_id}")

    # --- 1. Data Fetching ---
    logger.info("Step 1: Fetching data...")
    model_organism_data = pull_model_organism_orthologs(gene_query)
    if not model_organism_data:
        logger.error("Failed to pull model organism orthologs. Aborting pipeline.")
        return None
    
    # Extract the list of orthologs in model organisms
    orthologs_list = model_organism_data.get("orthologs_in_model_organisms")
    if not orthologs_list:
        logger.error("'orthologs_in_model_organisms' key not found in fetched data. Aborting pipeline.")
        return None

    og_fasta_data = pull_OG_fasta(ortholog_group_id)
    if not og_fasta_data:
        logger.error(f"Failed to pull FASTA data for ortholog group {ortholog_group_id}. Aborting pipeline.")
        return None

    all_ortholog_records = fasta_to_seqrecord(og_fasta_data)
    if not all_ortholog_records:
        logger.error("Failed to convert FASTA data to SeqRecord objects. Aborting pipeline.")
        return None

    # --- 2. Sequence Selection and Processing ---
    logger.info("Step 2: Processing sequences...")
    model_organism_gene_ids = get_model_organism_genes(orthologs_list) # Pass the list
    if not model_organism_gene_ids:
        logger.warning("No model organism gene IDs extracted. Proceeding with available records.")
    
    model_organism_sequences = select_genes_from_seqrecord(model_organism_gene_ids, all_ortholog_records)
    if not model_organism_sequences:
        logger.warning("No sequences selected for model organisms. This might affect taxa selection.")

    # Filter for one sequence per major taxon by selecting the longest
    representative_taxa_ids = select_big_taxa_seq(model_organism_sequences)
    representative_taxa_sequences = select_genes_from_seqrecord(representative_taxa_ids, all_ortholog_records)
    if not representative_taxa_sequences:
        logger.warning("No representative taxa sequences selected. This might affect k-selection.")
        # Fallback to using model_organism_sequences if representative_taxa_sequences is empty
        # and k_sequences > 1 (we still need the target sequence)
        if k_sequences > 1:
            logger.info("Falling back to use all model organism sequences for k-selection due to empty representative set.")
            representative_taxa_sequences = model_organism_sequences


    # Select top K sequences including the target gene
    if target_gene_orthodb_id not in all_ortholog_records:
        logger.error(f"Target gene ID {target_gene_orthodb_id} not found in the pulled ortholog group records. Aborting pipeline.")
        return None

    # If representative_taxa_sequences is empty but k > 1, we might still want to proceed with just the target
    # or a very small set. select_biggest_k_seq should handle this.
    # If representative_taxa_sequences is empty and k=1, it should just return the target.
    if not representative_taxa_sequences and k_sequences > 1:
         logger.warning("Representative taxa sequences are empty. Selection for k_sequences might be limited.")
         # Provide an empty dict for records_to_filter if it's truly empty to avoid errors in select_biggest_k_seq
         # if it expects a dict. The function should handle this gracefully.

    selected_k_sequences = select_biggest_k_seq(
        k_sequences,
        representative_taxa_sequences, # This might be empty, select_biggest_k_seq should handle it
        target_gene_orthodb_id,
        all_ortholog_records
    )

    if not selected_k_sequences or len(selected_k_sequences) < 2 : # PAML needs at least 2 sequences
        logger.error(f"Not enough sequences selected ({len(selected_k_sequences)}) for analysis (need at least 2). Aborting pipeline.")
        return None

    padded_sequences = pad_sequence_lengths(selected_k_sequences)
    processed_sequences = remove_stop_codons_in_multiple(padded_sequences)
    
    # Convert IDs to organism names (optional step, can be commented out if not desired)
    # This step was in the original app.py's main logic
    named_sequences = convert_organism_id_to_names(processed_sequences)

    # Format IDs for CODEML and create alignment object
    alignment = format_ids_and_create_alignment(named_sequences) # Use named_sequences here
    
    if not alignment or len(alignment) < 2:
        logger.error("Alignment creation failed or resulted in an alignment with less than 2 sequences. Aborting pipeline.")
        return None

    # --- 3. File Output for PAML ---
    logger.info("Step 3: Writing files for PAML...")
    # Ensure output directories exist if not using default CACHE_DIR for these
    output_dir_phylip = os.path.dirname(output_phylip_file)
    if output_dir_phylip and not os.path.exists(output_dir_phylip):
        os.makedirs(output_dir_phylip, exist_ok=True)
    
    output_dir_newick = os.path.dirname(output_newick_file)
    if output_dir_newick and not os.path.exists(output_dir_newick):
        os.makedirs(output_dir_newick, exist_ok=True)

    write_phylip_manual(alignment, output_phylip_file)
    
    # --- 4. Phylogenetic Tree Construction ---
    logger.info("Step 4: Constructing phylogenetic tree...")
    tree = construct_phylo_tree(alignment, rooted=rooted_tree)
    if tree:
        write_newick_tree_with_header(tree, output_newick_file)
    else:
        logger.error("Phylogenetic tree construction failed. Aborting CODEML analysis.")
        return None

    # --- 5. Run CODEML Analysis ---
    logger.info("Step 5: Running CODEML Analysis...")
    codeml_results = run_codeml_positive_selection(
        tree_filepath=output_newick_file,
        seqfile_filepath=output_phylip_file,
        output_dir=paml_output_dir,
        cleanup_working_dir=cleanup_paml_dir
    )

    if codeml_results:
        logger.info("CODEML analysis completed successfully.")
        # Log the main results (e.g., lnL values)
        for model, data in codeml_results.items():
            logger.info(f"  Model {model}: lnL = {data.get('lnL')}")
            # You might want to log more details or specific parameters here
    else:
        logger.error("CODEML analysis failed or produced no results.")
    
    logger.info("Phylogenetic analysis pipeline finished.")
    return codeml_results

if __name__ == "__main__":
    # Example parameters based on your test.py and provided files
    # From WBGene00004963_search.json, for C. elegans (6239_0):
    # gene_id is "CELE_C17D12.6", param is "6239_0:000672"
    # One of the ortholog groups at Eukaryota level is "430340at2759"
    
    # Ensure the cache directory exists for fetching to work smoothly if not already created
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)
        logger.info(f"Created cache directory for main execution: {CACHE_DIR}")

    gene_to_analyze = "WBGene00004963"
    target_gene_id_orthodb = "6239_0:000672" # This is the OrthoDB ID for the C. elegans gene
    ortholog_group = "430340at2759" # An example ortholog group ID
    num_sequences_for_analysis = 10

    # Define output file paths (can be customized)
    phylip_file = os.path.join(DEFAULT_PAML_OUTPUT_DIR, "alignment.phylip")
    newick_file = os.path.join(DEFAULT_PAML_OUTPUT_DIR, "tree.nwk")
    
    # Create PAML output directory if it doesn't exist, as PAML will write files there
    if not os.path.exists(DEFAULT_PAML_OUTPUT_DIR):
        os.makedirs(DEFAULT_PAML_OUTPUT_DIR, exist_ok=True)

    logger.info(f"--- Running Full Phylogenetic Pipeline for {gene_to_analyze} ---")
    results = run_phylogenetic_pipeline(
        gene_query=gene_to_analyze,
        target_gene_orthodb_id=target_gene_id_orthodb,
        ortholog_group_id=ortholog_group,
        k_sequences=num_sequences_for_analysis,
        output_phylip_file=phylip_file,
        output_newick_file=newick_file,
        paml_output_dir=DEFAULT_PAML_OUTPUT_DIR, # PAML output will go here
        cleanup_paml_dir=False, # Set to False for debugging to inspect PAML files
        rooted_tree=False # Use NJ tree
    )

    if results:
        logger.info("Pipeline completed. CODEML Results:")
        for model, data in results.items():
            logger.info(f"  Model: {model}")
            logger.info(f"    lnL: {data.get('lnL')}")
            parameters = data.get('parameters', {})
            if parameters:
                logger.info("    Parameters:")
                for param_name, param_value in parameters.items():
                    if isinstance(param_value, dict):
                        logger.info(f"      {param_name}:")
                        for sub_key, sub_value in param_value.items():
                            logger.info(f"        {sub_key}: {sub_value}")
                    else:
                        logger.info(f"      {param_name}: {param_value}")
            # Example: Accessing specific parameters for M2a if interested in positive selection
            if model == "M2" and parameters: # Note: Biopython might label M2a as 'M2'
                if 'site classes' in parameters:
                    for sc in parameters['site classes']:
                        if sc.get('omega', 0) > 1:
                            logger.info(f"      Potential positive selection: w = {sc['omega']:.3f} with proportion p = {sc['proportion']:.3f}")
    else:
        logger.error("Pipeline did not complete successfully or CODEML analysis failed.")

    logger.info("--- Main Pipeline Execution Finished ---")

