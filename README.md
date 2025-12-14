# AUTOMATED EVOLUTIONARY PRESSURE ANALYSIS (AutoEPA Pipeline) (Working title)

This is a little pipeline project with the aim of melding existing bioinformatics tools into an automatic pipeline. In particular, running CODEML positive selection analysis on a selected gene of interest.

# Requirements:

In no particular order, you will need the following CLI tools installed (until I can containerize this stuff):
- CODEML
- pal2nal
- IQTREE3
- MUSCLE

You will also need the following python libraries:
- Requests
- Biopython

# How it works:

Heavily simplified steps:
- You specify an NCBI gene id (eg. 173042) and taxa level id (eg. 6231).
- Orthologs are retrieved from OrthoDB
- 10 nucleotide sequences (9 Orthologs + your gene) are converted to protein sequences and aligned with MUSCLE
- The MSA is converted back to nucleotides and piped into IQTREE3 to infer the tree
- The MSA and tree are fed into CODEML for positive selection analysis.

# Quickstart:

Yeah this doesn't exist yet.

# Installation:

Good luck.




