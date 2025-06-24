from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from typing import List

import requests, io

def compile_orthologs(
    ncbi_gene_id: str,
    taxanomic_level_id: str,
    ) -> List[SeqRecord]:

    """
    Compiles a list of orthologs to the specified gene at the specified taxanomic level.

    Args:
        ncbi_gene_id: What is the id of the ID of the gene you want to compile orthologs for?
        taxanomic_level_id: At what taxanomic level do you want to look?

    Returns:
        A list of seqRecord objects, with seq, id, name and description -> id is orthoDB_id:orthoDB_OG_id I think, name is the actual name

    """

    orthoDB_API = "https://data.orthodb.org/v12/" 

    # STEP 1: Query OrthoDB for an ortholog group ID and take the first available one
    ortholog_group_id = requests.get(f"{orthoDB_API}search", params={
        'ncbi': ncbi_gene_id,
        'level': taxanomic_level_id,
        }).json()['data'][0]

    # STEP 2: Query OrthoDB for ortholog group fasta data
    ortholog_fasta_compilation = requests.get(f"{orthoDB_API}fasta", params={
        'id': ortholog_group_id,
        'seqtype': 'cds',
        }).text
    
    # STEP 3: convert fasta text to list of seqRecords and return it (setting the name to the actual name instead of the id again with some magic)
    ortholog_compilation = [SeqRecord(record.seq, id=record.id, name=record.description.split('organism_name":"')[1].split('"')[0], description=record.description) for record in SeqIO.parse(io.StringIO(ortholog_fasta_compilation), "fasta")]
    return ortholog_compilation
    
if __name__ == '__main__':

    test_ncbi_ids = ["173042", "173402", "3039"] # C. elegans spe-39, C. elegans lin-39, Human HBA1
    test_levels = ["6231", "6231", "40674"]    # Nematoda, Nematoda, Mammalia

    print(compile_orthologs(test_ncbi_ids[0], test_levels[0])[0])

