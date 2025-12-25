import argparse


def parse() -> argparse.Namespace:
    """
    Parses input to determine which pipeline to run and how to run it.

    Args:
        None. What is this magic?

    Returns:
        Should return a list of arguments? I'm not sure how this works yet.

    """
    parser = argparse.ArgumentParser(
              prog="Phylogenetic Evolutionary Analysis Automation Logic (PEAAL)",
              description="A not so simple command line tool for running an optimized positive selection analysis on a gene of interest."
              )
    parser.add_argument('gene_id', nargs='?', default=None, help='The NCBI ID of the gene you want to analyze. Leave blank if you are using a pre-selected list of genes with -s')
    parser.add_argument('-s', '--specified', default=None, help='Use this if you already have a list of pre-selected ortholog NCBI IDs you would like to specify. Separate each ID with a comma, no spaces') ## for specifying your own list of orthologs 
    parser.add_argument('-t', '--taxa', help='Specify the taxanomic level ID at which you would like to conduct the analysis.')
    return parser.parse_args()
  

if __name__ == '__main__':
    print(parse())
