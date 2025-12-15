from typing import Dict, Any, Optional
from pipe_fun import (
    select_ortholog_group,
    compile_orthologs,
    select_orthologs,
    convert_to_proteins,
    run_muscle,
    run_pal2nal,
    convert_colon2dash,
    run_iqtree,
    run_codeml
)

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

    ortholog_group_id = select_ortholog_group(ncbi_gene_id, taxanomic_level_id) # This returns an ortholog group ID.
    ortholog_compilation = compile_orthologs(ortholog_group_id) # This returns a seqRecord list of orthologs.
    selection_path = select_orthologs(ortholog_compilation, k_sequences, ortholog_group_id) # This generates ortholog CDS fasta, returns path to file.
    proteins_path = convert_to_proteins(selection_path) # This converts to protein fasta, returns path to file.
    aligned_proteins_path = run_muscle(proteins_path) # This runs MUSCLE, returns protein MSA filepath.
    paml_path = convert_colon2dash(run_pal2nal(aligned_proteins_path, selection_path)) # This converts and returns a CDS MSA paml format file.
    fasta_path = convert_colon2dash(run_pal2nal(aligned_proteins_path, selection_path, "fasta")) # Same as above, except fasta format file, for IQTREE3
    treefile_path = run_iqtree(fasta_path) # This runs IQTREE3 on the new fasta, returns a path to the tree file.
    results_path = run_codeml(paml_path, treefile_path) # This runs CODEML on the tree file and paml file, returns a results file path.
    # results = parse_results(results_path) # This returns a dictionary of relevant results from the results file.

    return results_path




if __name__ == '__main__':

    test_ncbi_ids = ["173042", "173402", "3039"] # C. elegans spe-39, C. elegans lin-39, Human HBA1
    test_levels = ["6231", "6231", "40674"]    # Nematoda, Nematoda, Mammalia

    print(run_pipeline(test_ncbi_ids[0]))
 

    
