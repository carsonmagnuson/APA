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



def convert_colon2dash(filepath: str) -> str:
    """
    Converts colons in the name to dashes. Primarily for iqtree, because iqtree uses colons for other stuff.

    Args:
        filepath: What is the name of the file you're converting names in?

    Returns:
        Converted filepath.

    """

    # STEP 1: Read in the file while replacing colons
    # lines = list(list(char if char != ":" else '-' for char in l) for l in open(filepath, "r").splitlines())
    lines = [l.replace(':', '-') for l in open(filepath, "r").read().splitlines()]

    # STEP 2: Write out to a fasta file and return the filepath
    output_filepath = f"{filepath.split('.')[0]}_colon2dash.{filepath.split('.')[1]}"
    
    with open(output_filepath, "w", encoding="utf-8") as f:
        for line in lines:
            f.write(line)
            f.write("\n")

    return output_filepath
