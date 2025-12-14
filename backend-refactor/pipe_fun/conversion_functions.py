from Bio import SeqIO
from Bio.Seq import Seq
import os
import json

def convert_to_proteins(fasta_filepath: str) -> str:
    """
    Converts CDS fasta to protein sequences fasta. Note: this assumes a +1 reading frame.

    Args:
        fasta_filepath: What is the name of the file you're pulling the CDS from?

    Returns:
        A filepath to an outputted fasta file.

    """

    # STEP 1: Read in the fasta file.
    records = list(SeqIO.parse(fasta_filepath, "fasta"))

    # STEP 2: Convert to protein sequences.
    for record in records:
        record.seq = record.seq.translate(to_stop=True)

    # STEP 3: Write out to a fasta file and return the filepath
    output_filepath = f"{fasta_filepath.split('.')[0]}_coverted2prot.fasta"
    SeqIO.write(records, output_filepath, "fasta")
    return output_filepath
