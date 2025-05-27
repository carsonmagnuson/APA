import logging
import os
import json
from typing import Dict, Any, Optional, List

# Assuming these custom modules are in the same directory or PYTHONPATH
from data_fetching import (
    pull_model_organism_orthologs,
    pull_OG_fasta,
    pull_ortholog_group_ids_by_gene_and_level
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
DEFAULT_OUTPUT_PHYLIP_FILENAME = 'output.phylip'
DEFAULT_OUTPUT_NEWICK_FILENAME = 'tree.nwk'
DEFAULT_PAML_BASE_OUTPUT_DIR = "paml_output" # Base directory for all PAML runs
DEFAULT_PAML_RESULTS_FILENAME = "results.out" # Default name for CODEML output

def sanitize_for_path(value: str) -> str:
    """Sanitizes a string to be safely used as a file or directory name component."""
    return "".join(c if c.isalnum() else "_" for c in str(value))

def run_phylogenetic_pipeline(
    ncbi_id: str,
    k_sequences: int = 10,
    ortho_level_tax_id: str = "2759",
    paml_base_output_dir: str = DEFAULT_PAML_BASE_OUTPUT_DIR,
    rooted_tree: bool = False
) -> Optional[Dict[str, Any]]:
    """
    Executes the full phylogenetic analysis pipeline using an NCBI ID.
    Files for PAML will be stored in a run-specific subdirectory.
    """
    
    safe_ncbi_id = sanitize_for_path(ncbi_id)
    safe_ortho_level_tax_id = sanitize_for_path(ortho_level_tax_id)
    
    run_identifier = f"run_{safe_ncbi_id}_lvl{safe_ortho_level_tax_id}"
    run_specific_paml_dir = os.path.join(paml_base_output_dir, run_identifier)
    
    logger.info(f"Starting phylogenetic analysis pipeline for NCBI ID: {ncbi_id}, Orthology Level TaxID: {ortho_level_tax_id}")
    logger.info(f"PAML outputs for this run will be in: {run_specific_paml_dir}")

    # --- 0a. Fetch initial gene information (for target_gene_orthodb_id and model organism list) ---
    logger.info(f"Step 0a: Fetching gene details for NCBI ID: {ncbi_id} using OrthoDB /genesearch...")
    initial_gene_data = pull_model_organism_orthologs(gene_query=ncbi_id)

    if not initial_gene_data or "gene" not in initial_gene_data: # Check if 'gene' key exists
        logger.error(f"Failed to fetch or parse essential gene data for NCBI ID '{ncbi_id}'. Aborting pipeline.")
        return {"status": "PIPELINE_ERROR", "message": f"Could not retrieve essential gene data for NCBI ID {ncbi_id} from OrthoDB /genesearch."}

    try:
        target_gene_details = initial_gene_data["gene"]
        target_gene_orthodb_id = target_gene_details.get("gene_id", {}).get("param")
        if not target_gene_orthodb_id:
            logger.error(f"Failed to extract OrthoDB gene ID (param) for NCBI ID '{ncbi_id}'.")
            return {"status": "PIPELINE_ERROR", "message": "Could not derive OrthoDB gene ID (param) from /genesearch result."}
        logger.info(f"Derived OrthoDB gene ID: '{target_gene_orthodb_id}' for NCBI ID '{ncbi_id}'.")
    except (KeyError, TypeError) as e:
        logger.error(f"Error parsing initial gene data for NCBI ID '{ncbi_id}': {e}", exc_info=True)
        return {"status": "PIPELINE_ERROR", "message": f"Error parsing data from OrthoDB /genesearch for NCBI ID {ncbi_id}."}

    # --- 0b. Fetch Ortholog Group ID using the /search endpoint ---
    logger.info(f"Step 0b: Fetching Ortholog Group ID for NCBI ID '{ncbi_id}' at level '{ortho_level_tax_id}' using OrthoDB /search...")
    og_ids_list = pull_ortholog_group_ids_by_gene_and_level(
        gene_identifier=ncbi_id, level_tax_id=ortho_level_tax_id, identifier_type="ncbi"
    )

    if og_ids_list is None:
        logger.error(f"Failed to retrieve ortholog group IDs for NCBI ID '{ncbi_id}' at level '{ortho_level_tax_id}'. API call might have failed.")
        return {"status": "PIPELINE_ERROR", "message": f"Could not fetch OG IDs for {ncbi_id}."}
    if not og_ids_list:
        logger.error(f"No ortholog group IDs found for NCBI ID '{ncbi_id}' at level '{ortho_level_tax_id}'.")
        return {"status": "PIPELINE_ERROR", "message": f"No OrthoDB Ortholog Groups found for NCBI ID {ncbi_id} at taxonomic level {ortho_level_tax_id}."}
    
    ortholog_group_id = og_ids_list[0]
    logger.info(f"Selected Ortholog Group ID: {ortholog_group_id} for NCBI ID '{ncbi_id}' at level '{ortho_level_tax_id}'. (Found {len(og_ids_list)} OGs, using the first one).")
    
    # --- 1. Data Fetching (FASTA sequences for the selected OG) ---
    logger.info(f"Step 1: Fetching FASTA sequences for Ortholog Group ID: {ortholog_group_id}...")
    og_fasta_data = pull_OG_fasta(ortholog_group_id)
    if not og_fasta_data:
        logger.error(f"Failed to pull FASTA data for selected ortholog group {ortholog_group_id}. Aborting pipeline.")
        return {"status": "PIPELINE_ERROR", "message": f"Could not fetch FASTA for OG {ortholog_group_id}."}

    all_ortholog_records = fasta_to_seqrecord(og_fasta_data)
    if not all_ortholog_records:
        logger.error(f"Failed to convert FASTA data for OG {ortholog_group_id} to SeqRecord objects. Aborting pipeline.")
        return {"status": "PIPELINE_ERROR", "message": f"FASTA parsing error for OG {ortholog_group_id}."}

    if target_gene_orthodb_id not in all_ortholog_records:
        logger.error(f"Target OrthoDB gene ID '{target_gene_orthodb_id}' (derived from NCBI ID {ncbi_id}) "
                     f"not found in the FASTA records for the selected Ortholog Group '{ortholog_group_id}'. "
                     "This might indicate the selected OG or taxonomic level doesn't include the specific input gene, "
                     "or there's an issue with ID mapping. Try a broader taxonomic level for 'ortho_level_tax_id'.")
        return {"status": "PIPELINE_ERROR", 
                "message": f"Target OrthoDB gene ID {target_gene_orthodb_id} not found in selected OG {ortholog_group_id} FASTA."}

    # --- 2. Sequence Selection and Processing ---
    logger.info("Step 2: Processing sequences...")
    
    orthologs_list_in_models = initial_gene_data.get("orthologs_in_model_organisms") 

    if orthologs_list_in_models:
        logger.info("Found 'orthologs_in_model_organisms'. Proceeding with model organism-based selection strategy.")
        model_organism_gene_ids_from_search = get_model_organism_genes(orthologs_list_in_models)
        model_organism_sequences_from_og = select_genes_from_seqrecord(model_organism_gene_ids_from_search, all_ortholog_records)
        
        representative_taxa_ids = select_big_taxa_seq(model_organism_sequences_from_og)
        representative_taxa_sequences = select_genes_from_seqrecord(representative_taxa_ids, all_ortholog_records)
        
        source_for_k_selection = representative_taxa_sequences if representative_taxa_sequences else all_ortholog_records
        if not representative_taxa_sequences:
            logger.info("No representative taxa sequences derived from model organisms. Using all sequences from the OG for k-selection.")
    else:
        logger.warning(f"'orthologs_in_model_organisms' not found for NCBI ID {ncbi_id}. "
                       f"Falling back to selecting from all sequences in OG: {ortholog_group_id}.")
        representative_taxa_ids = select_big_taxa_seq(all_ortholog_records) # Select representatives from all OG sequences
        representative_taxa_sequences = select_genes_from_seqrecord(representative_taxa_ids, all_ortholog_records)
        source_for_k_selection = representative_taxa_sequences if representative_taxa_sequences else all_ortholog_records

    selected_k_sequences = select_biggest_k_seq(
        k_sequences,
        source_for_k_selection, 
        target_gene_orthodb_id,
        all_ortholog_records    
    )

    if not selected_k_sequences or len(selected_k_sequences) < 2:
        logger.error(f"Not enough sequences selected ({len(selected_k_sequences) if selected_k_sequences else 0}) for analysis (need at least 2). Target gene: '{target_gene_orthodb_id}'. Aborting pipeline.")
        return {"status": "PIPELINE_ERROR", "message": "Not enough sequences for analysis after selection."}

    padded_sequences = pad_sequence_lengths(selected_k_sequences)
    processed_sequences = remove_stop_codons_in_multiple(padded_sequences)
    alignment = format_ids_and_create_alignment(processed_sequences)

    if not alignment or len(alignment) < 2:
        logger.error("Alignment creation failed or resulted in an alignment with less than 2 sequences. Aborting pipeline.")
        return {"status": "PIPELINE_ERROR", "message": "Alignment creation failed."}

    # --- 3. File Output for PAML ---
    logger.info(f"Step 3: Writing PAML input files to {run_specific_paml_dir}...")
    if not os.path.exists(run_specific_paml_dir):
        try:
            os.makedirs(run_specific_paml_dir)
            logger.info(f"Created run-specific PAML directory: {run_specific_paml_dir}")
        except OSError as e:
            logger.error(f"Error creating PAML run directory {run_specific_paml_dir}: {e}", exc_info=True)
            return {"status": "PIPELINE_ERROR", "message": f"Could not create PAML run directory {run_specific_paml_dir}."}

    final_phylip_path = os.path.join(run_specific_paml_dir, DEFAULT_OUTPUT_PHYLIP_FILENAME)
    final_newick_path = os.path.join(run_specific_paml_dir, DEFAULT_OUTPUT_NEWICK_FILENAME)

    write_phylip_manual(alignment, final_phylip_path)

    # --- 4. Phylogenetic Tree Construction ---
    logger.info("Step 4: Constructing phylogenetic tree...")
    tree = construct_phylo_tree(alignment, rooted=rooted_tree)
    if tree:
        write_newick_tree_with_header(tree, final_newick_path)
    else:
        logger.error("Phylogenetic tree construction failed. Aborting CODEML analysis.")
        return {"status": "PIPELINE_ERROR", "message": "Phylogenetic tree construction failed."}

    # --- 5. Run CODEML Analysis ---
    logger.info(f"Step 5: Running CODEML for NCBI ID {ncbi_id} (OG: {ortholog_group_id})...")
    codeml_results_dict = run_codeml_positive_selection(
        tree_filepath=final_newick_path,
        seqfile_filepath=final_phylip_path,
        output_dir=run_specific_paml_dir,
        results_filename=DEFAULT_PAML_RESULTS_FILENAME, # Pass the constant here
        cleanup_working_dir=False 
    )

    final_pipeline_result = {
        "input_ncbi_id": ncbi_id,
        "target_gene_orthodb_id": target_gene_orthodb_id,
        "ortholog_group_id_used": ortholog_group_id,
        "ortho_level_tax_id_used": ortho_level_tax_id,
        "k_sequences_requested": k_sequences,
        "sequences_in_alignment": len(alignment) if alignment else 0,
        "paml_run_identifier": run_identifier, 
        "paml_results_filename": codeml_results_dict.get("results_filename", DEFAULT_PAML_RESULTS_FILENAME),
        "codeml_analysis_summary": codeml_results_dict
    }
    
    if codeml_results_dict.get("status", "").startswith("CODEML_SUCCESS"):
        logger.info(f"CODEML analysis for {ncbi_id} completed with status: {codeml_results_dict.get('status')}.")
    else:
        logger.warning(f"CODEML analysis for {ncbi_id} completed with issues: {codeml_results_dict.get('status')} - {codeml_results_dict.get('message')}")

    logger.info(f"Phylogenetic analysis pipeline finished for NCBI ID {ncbi_id}.")
    return final_pipeline_result

if __name__ == "__main__":
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)
        logger.info(f"Created cache directory for main execution: {CACHE_DIR}")
    
    base_paml_dir = DEFAULT_PAML_BASE_OUTPUT_DIR
    if not os.path.exists(base_paml_dir):
        os.makedirs(base_paml_dir, exist_ok=True)
        logger.info(f"Ensured base PAML output directory exists: {base_paml_dir}")

    test_ncbi_ids = ["173042", "173402", "3039"] # C. elegans spe-39, C. elegans lin-39, Human HBA1
    test_levels = ["6231", "6231", "40674"]    # Nematoda, Nematoda, Mammalia
    
    num_sequences_for_analysis = 10

    for test_id, test_level in zip(test_ncbi_ids, test_levels):
        logger.info(f"\n--- Running Full Phylogenetic Pipeline for NCBI ID: {test_id} at Level: {test_level} ---")
        
        results = run_phylogenetic_pipeline(
            ncbi_id=test_id,
            k_sequences=num_sequences_for_analysis,
            ortho_level_tax_id=test_level,
            paml_base_output_dir=base_paml_dir,
            rooted_tree=False
        )

        if results:
            logger.info(f"Pipeline completed for NCBI ID {test_id}. Overall Status/Results:")
            try:
                logger.info(json.dumps(results, indent=2, ensure_ascii=False))
            except Exception as e:
                logger.error(f"Error serializing results to JSON: {e}")
                logger.info(str(results)) 

            codeml_summary = results.get("codeml_analysis_summary", {})
            if codeml_summary:
                logger.info(f"  CODEML Status: {codeml_summary.get('status')}")
                logger.info(f"  CODEML Message: {codeml_summary.get('message')}")
                if "results_file_path_on_server" in codeml_summary:
                    logger.info(f"  CODEML Results File (server path): {codeml_summary.get('results_file_path_on_server')}")
        else:
            logger.error(f"Pipeline for NCBI ID {test_id} did not complete successfully or returned None.")
        logger.info(f"--- Finished test for NCBI ID: {test_id} ---")

    logger.info("--- Main Pipeline Execution Finished ---")

