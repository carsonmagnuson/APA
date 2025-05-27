import logging
import os
from typing import Dict, Any, Optional, List
import json # Added import for json

# Assuming these custom modules are in the same directory or PYTHONPATH
from data_fetching import (
    pull_model_organism_orthologs,
    pull_OG_fasta,
    pull_ortholog_group_ids_by_gene_and_level # New import
)
from sequence_processing import (
    fasta_to_seqrecord,
    get_model_organism_genes,
    select_genes_from_seqrecord,
    select_big_taxa_seq,
    select_biggest_k_seq,
    pad_sequence_lengths,
    remove_stop_codons_in_multiple,
    convert_organism_id_to_names,
    format_ids_and_create_alignment
)
from phylogenetic_analysis import construct_phylo_tree, run_codeml_positive_selection
from file_utils import write_phylip_manual, write_newick_tree_with_header, CACHE_DIR

# --- Configure Logging ---
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
DEFAULT_PAML_OUTPUT_DIR = "paml_output"

def run_phylogenetic_pipeline(
    ncbi_id: str,
    k_sequences: int = 10,
    ortho_level_tax_id: str = "2759", # Default to Eukaryota (NCBI TaxID for Eukaryota)
    output_phylip_file: str = DEFAULT_OUTPUT_PHYLIP_FILE,
    output_newick_file: str = DEFAULT_OUTPUT_NEWICK_FILE,
    paml_output_dir: str = DEFAULT_PAML_OUTPUT_DIR,
    cleanup_paml_dir: bool = True,
    rooted_tree: bool = False
) -> Optional[Dict[str, Any]]:
    """
    Executes the full phylogenetic analysis pipeline using an NCBI ID.

    Args:
        ncbi_id: The NCBI gene identifier (e.g., "173042" for a C. elegans gene).
        k_sequences: The number of sequences to include in the analysis (including the target).
        ortho_level_tax_id: The NCBI taxonomy ID for the level at which to find ortholog groups.
                            Defaults to "2759" (Eukaryota).
        output_phylip_file: Path for the output PHYLIP alignment file.
        output_newick_file: Path for the output Newick tree file.
        paml_output_dir: Directory for PAML output files.
        cleanup_paml_dir: Whether to remove the PAML output directory after completion.
        rooted_tree: Whether to construct a rooted (UPGMA) or unrooted (NJ) tree.

    Returns:
        A dictionary containing CODEML results or error information.
    """
    logger.info(f"Starting phylogenetic analysis pipeline for NCBI ID: {ncbi_id}, Orthology Level TaxID: {ortho_level_tax_id}")

    # --- 0a. Fetch initial gene information using NCBI ID to get OrthoDB gene ID ---
    logger.info(f"Step 0a: Fetching gene details for NCBI ID: {ncbi_id} using OrthoDB /genesearch...")
    initial_gene_data = pull_model_organism_orthologs(gene_query=ncbi_id)

    if not initial_gene_data or "gene" not in initial_gene_data:
        logger.error(f"Failed to fetch or parse initial gene data for NCBI ID '{ncbi_id}'. Aborting pipeline.")
        return {"status": "PIPELINE_ERROR", "message": f"Could not retrieve initial gene data for NCBI ID {ncbi_id} from OrthoDB /genesearch."}

    try:
        target_gene_details = initial_gene_data["gene"]
        # This is the OrthoDB internal ID for the specific gene, e.g., "6239_0:000672"
        target_gene_orthodb_id = target_gene_details.get("gene_id", {}).get("param")
        
        if not target_gene_orthodb_id:
            logger.error(f"Failed to extract OrthoDB gene ID (param) for NCBI ID '{ncbi_id}'.")
            logger.debug(f"Gene details received: {target_gene_details}")
            return {"status": "PIPELINE_ERROR", "message": "Could not derive OrthoDB gene ID (param) from /genesearch result."}
        logger.info(f"Derived OrthoDB gene ID for NCBI ID '{ncbi_id}': '{target_gene_orthodb_id}'")

    except (KeyError, TypeError) as e:
        logger.error(f"Error parsing initial gene data (for OrthoDB gene ID) for NCBI ID '{ncbi_id}': {e}", exc_info=True)
        logger.debug(f"Initial gene data received: {initial_gene_data}")
        return {"status": "PIPELINE_ERROR", "message": f"Error parsing /genesearch data for NCBI ID {ncbi_id}."}

    # --- 0b. Fetch Ortholog Group ID using the /search endpoint ---
    logger.info(f"Step 0b: Fetching Ortholog Group ID for NCBI ID '{ncbi_id}' at level '{ortho_level_tax_id}' using OrthoDB /search...")
    # We use the original ncbi_id for the search, as it's the primary input identifier.
    og_ids_list = pull_ortholog_group_ids_by_gene_and_level(
        gene_identifier=ncbi_id, 
        level_tax_id=ortho_level_tax_id,
        identifier_type="ncbi" # Explicitly state we are using an NCBI ID for the search
    )

    if og_ids_list is None: # Indicates an API error or malformed response
        logger.error(f"Failed to retrieve ortholog group IDs for NCBI ID '{ncbi_id}' at level '{ortho_level_tax_id}'. API call might have failed.")
        return {"status": "PIPELINE_ERROR", "message": f"Could not fetch OG IDs for {ncbi_id} at level {ortho_level_tax_id}."}
    
    if not og_ids_list: # Empty list means no OGs found matching the criteria
        logger.error(f"No ortholog group IDs found for NCBI ID '{ncbi_id}' at level '{ortho_level_tax_id}'.")
        return {"status": "PIPELINE_ERROR", "message": f"No OrthoDB Ortholog Groups found for NCBI ID {ncbi_id} at taxonomic level {ortho_level_tax_id}."}

    # Strategy: Use the first OG ID returned. 
    # OrthoDB's /search can return multiple OGs if the gene is part of several at that level.
    # For now, we assume the first one is the most relevant or a good starting point.
    # More sophisticated selection might be needed if multiple OGs are common.
    ortholog_group_id = og_ids_list[0]
    logger.info(f"Selected Ortholog Group ID: {ortholog_group_id} for NCBI ID '{ncbi_id}' at level '{ortho_level_tax_id}'. (Found {len(og_ids_list)} OGs, using the first one).")
    
    # --- 1. Data Fetching (FASTA sequences for the selected OG) ---
    logger.info(f"Step 1: Fetching FASTA sequences for Ortholog Group ID: {ortholog_group_id}...")
    
    # The `initial_gene_data` (from /genesearch) is still useful for `orthologs_in_model_organisms` if needed for other purposes,
    # but for the main alignment, we use the sequences from the `ortholog_group_id` obtained via /search.
    model_organism_data = initial_gene_data 
    orthologs_list_in_models = model_organism_data.get("orthologs_in_model_organisms") # This is from /genesearch
    if not orthologs_list_in_models:
        logger.warning("'orthologs_in_model_organisms' key not found in /genesearch data. This might affect some downstream logic if it relies on this specific list.")
        # Depending on the logic, this might not be critical if the main analysis uses the broader OG.

    og_fasta_data = pull_OG_fasta(ortholog_group_id) # Now uses the correctly fetched OG ID
    if not og_fasta_data:
        logger.error(f"Failed to pull FASTA data for selected ortholog group {ortholog_group_id}. Aborting pipeline.")
        return {"status": "PIPELINE_ERROR", "message": f"Could not fetch FASTA for OG {ortholog_group_id}."}

    all_ortholog_records = fasta_to_seqrecord(og_fasta_data)
    if not all_ortholog_records:
        logger.error(f"Failed to convert FASTA data for OG {ortholog_group_id} to SeqRecord objects. Aborting pipeline.")
        return {"status": "PIPELINE_ERROR", "message": f"FASTA parsing error for OG {ortholog_group_id}."}

    # Critical Check: Ensure the target_gene_orthodb_id (e.g., "6239_0:000672") is in the pulled FASTA
    if target_gene_orthodb_id not in all_ortholog_records:
        logger.error(f"Target OrthoDB gene ID '{target_gene_orthodb_id}' (derived from NCBI ID {ncbi_id}) "
                     f"not found in the FASTA records for the selected Ortholog Group '{ortholog_group_id}'. "
                     "This may indicate the selected OG or taxonomic level doesn't include the specific input gene, "
                     "or there's an issue with ID mapping. Try a broader taxonomic level for 'ortho_level_tax_id'.")
        if all_ortholog_records:
            logger.debug(f"Sample record IDs from OG '{ortholog_group_id}': {list(all_ortholog_records.keys())[:5]}")
        else:
            logger.debug(f"No records found in all_ortholog_records for OG '{ortholog_group_id}'.")
        return {"status": "PIPELINE_ERROR", 
                "message": f"Target OrthoDB gene ID {target_gene_orthodb_id} not found in selected OG {ortholog_group_id} FASTA. Consider a different 'ortho_level_tax_id'."}

    # --- 2. Sequence Selection and Processing ---
    logger.info("Step 2: Processing sequences...")
    
    # Get gene IDs from the model organisms list (obtained from /genesearch)
    # These are OrthoDB gene IDs (e.g., "species_ncbi_id_0:gene_locus_tag")
    model_organism_gene_ids_from_search = get_model_organism_genes(orthologs_list_in_models or []) # Pass empty list if None
    
    # Select these model organism sequences from the broader OG FASTA data
    # This ensures we are working with sequences that are part of the chosen OG
    model_organism_sequences_from_og = select_genes_from_seqrecord(model_organism_gene_ids_from_search, all_ortholog_records)
    
    if not model_organism_sequences_from_og:
        logger.warning("No sequences for model organisms (from /genesearch result) were found within the selected Ortholog Group. "
                       "The k-selection will proceed using all sequences from the OG.")
        # If no model organism sequences are found in the OG,
        # representative_taxa_sequences will be based on all_ortholog_records.
        # This might happen if the chosen OG level is very specific and doesn't include those model organisms,
        # or if the gene IDs don't match.
    
    # Select one representative sequence per taxon from the model organism sequences found in the OG
    # If model_organism_sequences_from_og is empty, this will also be empty.
    representative_taxa_ids = select_big_taxa_seq(model_organism_sequences_from_og)
    representative_taxa_sequences = select_genes_from_seqrecord(representative_taxa_ids, all_ortholog_records)

    # If representative_taxa_sequences is empty (e.g., no model organisms found in OG or no diversity),
    # or if we need more sequences for k-selection, fall back to using all sequences from the OG.
    # The `select_biggest_k_seq` function will handle picking from this larger pool.
    source_for_k_selection = representative_taxa_sequences if representative_taxa_sequences else all_ortholog_records
    if not representative_taxa_sequences:
        logger.info("No representative taxa sequences derived from model organisms. Using all sequences from the OG for k-selection.")


    selected_k_sequences = select_biggest_k_seq(
        k_sequences,
        source_for_k_selection, # Use representative_taxa_sequences or all_ortholog_records
        target_gene_orthodb_id, # This is the specific OrthoDB gene ID of our input gene
        all_ortholog_records    # This ensures the target gene is always sourced from the full OG dataset
    )

    if not selected_k_sequences or len(selected_k_sequences) < 2:
        logger.error(f"Not enough sequences selected ({len(selected_k_sequences) if selected_k_sequences else 0}) for analysis (need at least 2). Target gene: '{target_gene_orthodb_id}'. Aborting pipeline.")
        return {"status": "PIPELINE_ERROR", "message": "Not enough sequences for analysis after selection."}

    padded_sequences = pad_sequence_lengths(selected_k_sequences)
    processed_sequences = remove_stop_codons_in_multiple(padded_sequences)
    
    # Convert OrthoDB IDs (like "taxonID_instance:geneID") to organism names for tree display
    # This step should happen *before* format_ids_and_create_alignment if you want PAML IDs to be based on organism names.
    # However, format_ids_and_create_alignment will truncate and ensure alphanumeric, so the exact output might vary.
    # The current sequence_processing.convert_organism_id_to_names tries to parse from description.
    # Let's ensure the IDs used for PAML are derived consistently.
    # For now, we'll keep the flow, but be mindful that format_ids_and_create_alignment will create the final PAML IDs.
    named_sequences_for_display = convert_organism_id_to_names(processed_sequences.copy()) # Use a copy if original IDs are needed later

    # Format IDs specifically for PAML (alphanumeric, 10 chars) and create alignment object
    alignment = format_ids_and_create_alignment(processed_sequences) # Pass the records with OrthoDB IDs for PAML ID generation

    if not alignment or len(alignment) < 2:
        logger.error("Alignment creation failed or resulted in an alignment with less than 2 sequences. Aborting pipeline.")
        return {"status": "PIPELINE_ERROR", "message": "Alignment creation failed."}

    # --- 3. File Output for PAML ---
    logger.info("Step 3: Writing files for PAML...")
    # Ensure the output directory for PAML files exists
    if not os.path.exists(paml_output_dir):
        os.makedirs(paml_output_dir, exist_ok=True)
        logger.info(f"Created PAML output directory: {paml_output_dir}")

    # Construct full paths for output files within the paml_output_dir
    # This makes it easier if paml_output_dir is specified as something other than the default.
    final_phylip_path = os.path.join(paml_output_dir, os.path.basename(output_phylip_file))
    final_newick_path = os.path.join(paml_output_dir, os.path.basename(output_newick_file))

    write_phylip_manual(alignment, final_phylip_path)

    # --- 4. Phylogenetic Tree Construction ---
    logger.info("Step 4: Constructing phylogenetic tree...")
    tree = construct_phylo_tree(alignment, rooted=rooted_tree) # Alignment has PAML-formatted IDs
    if tree:
        # The tree written here will have PAML-formatted IDs.
        write_newick_tree_with_header(tree, final_newick_path)
    else:
        logger.error("Phylogenetic tree construction failed. Aborting CODEML analysis.")
        return {"status": "PIPELINE_ERROR", "message": "Phylogenetic tree construction failed."}

    # --- 5. Run CODEML Analysis ---
    logger.info(f"Step 5: Running CODEML Analysis for NCBI ID {ncbi_id} (OG: {ortholog_group_id})...")
    codeml_results = run_codeml_positive_selection(
        tree_filepath=final_newick_path,
        seqfile_filepath=final_phylip_path,
        output_dir=paml_output_dir, # PAML will use this as its working directory
        cleanup_working_dir=cleanup_paml_dir
    )

    # Add some context to the results
    final_pipeline_result = {
        "input_ncbi_id": ncbi_id,
        "target_gene_orthodb_id": target_gene_orthodb_id,
        "ortholog_group_id_used": ortholog_group_id,
        "ortho_level_tax_id_used": ortho_level_tax_id,
        "k_sequences_requested": k_sequences,
        "sequences_in_alignment": len(alignment) if alignment else 0,
        "codeml_analysis": codeml_results
    }


    if codeml_results:
        if "status" in codeml_results and ("FAILED" in codeml_results.get("status", "").upper() or "ERROR" in codeml_results.get("status", "").upper() or "INCOMPLETE" in codeml_results.get("status", "").upper()):
            logger.warning(f"CODEML analysis for {ncbi_id} completed with issues: {codeml_results.get('message')}")
        else:
            logger.info(f"CODEML analysis for {ncbi_id} completed (or Biopython parsing was successful).")
    else:
        logger.error(f"CODEML analysis for {ncbi_id} failed to return results. Check logs for details.")
        # Update the status in final_pipeline_result if codeml_results is None or indicates failure
        if not final_pipeline_result.get("codeml_analysis"):
             final_pipeline_result["codeml_analysis"] = {"status": "CODEML_ERROR", "message": "CODEML analysis did not return results."}


    logger.info(f"Phylogenetic analysis pipeline finished for NCBI ID {ncbi_id}.")
    return final_pipeline_result

if __name__ == "__main__":
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)
        logger.info(f"Created cache directory for main execution: {CACHE_DIR}")

    # Test case: C. elegans gene spe-39, NCBI Gene ID: 173042
    # OrthoDB ID for this gene: 6239_0:000672
    # A known Ortholog Group for this gene at Eukaryota level (2759) could be, for example, '430340at2759'
    # Let's test with a specific NCBI gene ID and a taxonomic level.
    test_ncbi_id = "173042"  # NCBI Gene ID for C. elegans spe-39
    # test_ortho_level = "2759" # Eukaryota
    test_ortho_level = "6231" # Nematoda - a more specific level to test if the target gene is found
    # test_ortho_level = "33208" # Metazoa - another level to test

    num_sequences_for_analysis = 10 # Requesting 10 sequences

    # Define output file paths within the PAML output directory
    paml_run_dir = os.path.join(DEFAULT_PAML_OUTPUT_DIR, f"run_{test_ncbi_id}_lvl{test_ortho_level}")
    
    # Clean up previous run directory if it exists to ensure a fresh test
    # if os.path.exists(paml_run_dir):
    #     import shutil
    #     logger.info(f"Removing previous test directory: {paml_run_dir}")
    #     shutil.rmtree(paml_run_dir)

    if not os.path.exists(paml_run_dir):
        os.makedirs(paml_run_dir, exist_ok=True)
        logger.info(f"Created PAML output directory for this run: {paml_run_dir}")

    phylip_file = os.path.join(paml_run_dir, DEFAULT_OUTPUT_PHYLIP_FILE)
    newick_file = os.path.join(paml_run_dir, DEFAULT_OUTPUT_NEWICK_FILE)

    logger.info(f"--- Running Full Phylogenetic Pipeline for NCBI ID: {test_ncbi_id} at Level: {test_ortho_level} ---")
    results = run_phylogenetic_pipeline(
        ncbi_id=test_ncbi_id,
        k_sequences=num_sequences_for_analysis,
        ortho_level_tax_id=test_ortho_level,
        output_phylip_file=phylip_file,
        output_newick_file=newick_file,
        paml_output_dir=paml_run_dir, # PAML runs within this specific subdirectory
        cleanup_paml_dir=False, # Set to False for debugging, True for production
        rooted_tree=False
    )

    if results:
        logger.info(f"Pipeline completed for NCBI ID {test_ncbi_id}. Overall Status/Results:")
        logger.info(json.dumps(results, indent=2))

        codeml_analysis_results = results.get("codeml_analysis", {})
        if codeml_analysis_results and not ("status" in codeml_analysis_results and "ERROR" in codeml_analysis_results.get("status", "")):
            logger.info("CODEML Analysis Summary:")
            for model, data in codeml_analysis_results.items():
                if model.startswith("M") and isinstance(data, dict): # Check if it's M0, M1a, M2a etc.
                    logger.info(f"  Model: {model}")
                    logger.info(f"    lnL: {data.get('lnL')}")
                    parameters = data.get('parameters', {})
                    if parameters:
                        logger.info("    Parameters:")
                        # Only print key parameters for brevity, or all if needed
                        logger.info(f"      kappa: {parameters.get('kappa')}")
                        logger.info(f"      omega: {parameters.get('omega')}") # For M0
                        if 'site classes' in parameters: # For M1a, M2a
                            logger.info("      Site Classes:")
                            for sc in parameters['site classes']:
                                logger.info(f"        p: {sc.get('proportion')}, w: {sc.get('omega')}")
                                if isinstance(sc.get('omega'), (int, float)) and sc.get('omega', 0) > 1:
                                    logger.info(f"        Potential positive selection: w = {sc.get('omega', 0):.3f} with proportion p = {sc.get('proportion', 0):.3f}")
                        # For M2a BEB results
                        if 'BEB' in parameters and 'positive sites' in parameters['BEB']:
                            logger.info("    BEB Positively Selected Sites:")
                            for site_info in parameters['BEB']['positive sites']:
                                logger.info(f"      Site: {site_info.get('site')} ({site_info.get('aa')}), Pr(w>1): {site_info.get('probability'):.3f}, mean w: {site_info.get('mean w'):.3f} +- {site_info.get('std error w'):.3f}")
        elif codeml_analysis_results:
            logger.warning(f"CODEML analysis had issues: {codeml_analysis_results.get('status')} - {codeml_analysis_results.get('message')}")
            if "results_file_path" in codeml_analysis_results:
                logger.info(f"  CODEML output file: {codeml_analysis_results.get('results_file_path')}")

    else:
        logger.error(f"Pipeline for NCBI ID {test_ncbi_id} did not complete successfully or returned None.")

    logger.info("--- Main Pipeline Execution Finished ---")

