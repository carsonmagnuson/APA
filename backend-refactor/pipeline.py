from typing import Dict, Any, Optional

def run_pipeline(
    ncbi_gene_id: str,
    k_sequences: int = 10,
    taxanomic_level_id: str = "2759",
    ) -> Optional[Dict[str, Any]]:
    """
    Runs the main pipeline.

    Args:
        ncbi_gene_id: What is the id of the ID of the gene you want to run analysis on?
        k_sequences: How many orthologs you want to compare in your analysis? Defaults to 10.
        taxanomic_level_id: At what taxanomic level is your analysis to occur? Defaults to NCBI taxon ID 2759 of Eukaryota.

    Returns:
        Dictionary of results parsed from results.out, format tbd.

    """

    ortholog_compilation = compile_orthologs(ncbi_gene_id, taxanomic_level_id) # This returns a seqRecord list
    selection_path = select_orthologs(ortholog_compilation, k_sequences) # This generates a fasta file, returns path to file
    phylip_path = run_muscle(selection_path) # This generates a phylip file, returns path to file
    results_path = run_codeml(phylip_path) # This generates a codeml results file, returns path to file
    results = parse_results(results_path) # This returns a dictionary of relevant results from the results file

    return results




    

    
