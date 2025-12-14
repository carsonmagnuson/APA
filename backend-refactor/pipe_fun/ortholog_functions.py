from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from typing import List
from external_tools import run_muscle, run_pal2nal
from conversion_functions import convert_to_proteins

import requests, io, os


def select_ortholog_group(
    ncbi_gene_id: str,
    taxanomic_level_id: str,
    ) -> str:
    """
    Determines an appropriate ortholog group given a specified gene and taxonomic level.

    Args:
        ncbi_gene_id: What is the id of the ID of the gene you want to find orthologs for?
        taxanomic_level_id: At what taxanomic level do you want to focus?

    Returns:
        An OrthoDB ortholog group ID

    """

    orthoDB_API = "https://data.orthodb.org/v12/" 

    # STEP 1: Query OrthoDB for an ortholog group ID and take the first available one.
    ortholog_group_id = requests.get(f"{orthoDB_API}search", params={
        'ncbi': ncbi_gene_id,
        'level': taxanomic_level_id,
        }).json()['data'][0]

    #STEP 2: Return the ID
    return ortholog_group_id

def compile_orthologs(
        ortholog_group_id: str
        ) -> List[SeqRecord]:

    """
    Compiles a list of orthologs from a specified ortholog group

    Args:
        ncbi_gene_id: What is the ID of the ortholog group you want to compile a list from?

    Returns:
        A list of seqRecord objects, with seq, id, name and description -> id is orthoDB_id:orthoDB_OG_id I think, name is the actual name

    """

    orthoDB_API = "https://data.orthodb.org/v12/" 


    # STEP 1: Query OrthoDB for ortholog group fasta data.
    ortholog_fasta_compilation = requests.get(f"{orthoDB_API}fasta", params={
        'id': ortholog_group_id,
        'seqtype': 'cds',
        }).text
    
    # STEP 2: convert fasta text to list of seqRecords and return it (setting the name to the actual name instead of the id again with some magic).
    ortholog_compilation = [SeqRecord(record.seq, id=record.id, name=record.description.split('organism_name":"')[1].split('"')[0], description=record.description) for record in SeqIO.parse(io.StringIO(ortholog_fasta_compilation), "fasta")]
    return ortholog_compilation


def select_orthologs(
        ortholog_compilation: List[SeqRecord],
        k_sequences: int,
        ortholog_group_id: str,
        preselected_sequences: list = [],
        ) -> str:
    """
    Selects orthologs from a compiled seqRecord list for positive selection analysis.

    Args:
        ortholog_compilation: What list of orthologs do you want to look through?
        k_sequences: How many orthologs do you want, total, for your analysis (preselected sequences are included)?
        ortholog_group_id: What is the ID of the ortholog group you're selecting from?
        preselected_sequences: Are there any orthologs (or at least the original gene to which all the others are orthologous) you want to preselect?

    Returns:
        A filepath to an outputted fasta file.

    """

    # STEP 1: Remove the preselected sequences from the ortholog compilation.
    for record in preselected_sequences:
        ortholog_compilation.remove(record)

    # STEP 2: Place ortholog compilation in order of most suitable to least suitable.
    # ortholog_compilation.sort(key = lambda x: -len(x.seq))

    # STEP 3: Extend preselected list with most suitable sequence candidates
    remaining_k = k_sequences - len(preselected_sequences)
    selected_orthologs = preselected_sequences + ortholog_compilation[0:remaining_k]

    # STEP 4: Check for/create a directory for output using the ortholog group id
    output_directory = ortholog_group_id + "_run"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # STEP 5: Write selected list to a fasta file and return the path. Using length of the selected list in the name.
    output_filepath = f"{output_directory}/{len(selected_orthologs)}_selected_orthologs_CDS.fasta"
    SeqIO.write(selected_orthologs, output_filepath, "fasta")
    return output_filepath
 

if __name__ == '__main__':

    test_ncbi_ids = ["173042", "173402", "3039"] # C. elegans spe-39, C. elegans lin-39, Human HBA1
    test_levels = ["6231", "6231", "40674"]    # Nematoda, Nematoda, Mammalia
    
    ortholog_group_id = select_ortholog_group(test_ncbi_ids[0], test_levels[0])
    print(ortholog_group_id)
    ortholog_compilation = compile_orthologs(ortholog_group_id)
    print(ortholog_compilation)
    selection_path = select_orthologs(ortholog_compilation, 10, ortholog_group_id)
    print(selection_path)
    proteins_path = convert_to_proteins(selection_path)
    print(proteins_path)
    aligned_proteins_path = run_muscle(proteins_path)
    print(aligned_proteins_path)
    phylip_path = run_pal2nal(aligned_proteins_path, selection_path)
    print(phylip_path)


    
    

