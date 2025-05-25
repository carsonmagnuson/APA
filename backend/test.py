from app import *
from Bio import AlignIO

def get_model_orthos():
    gene = "WBGene00004963" # choose a gene
    orthologs = pull_model_organism_orthologs(gene) # pulling the model organisms that have ortholog genes.
    records = fasta_to_seqrecord(pull_OG_fasta('430340at2759')) # this seems to be pulling the ortholog group organisms in general, and putting their sequences in a dict for future ref.
    model_organism_genes = get_model_organism_genes(orthologs) # this uses the model organism name from the model organism OG to get a list of all the potential ortholog genes available in model organisms.
    model_sequences = select_genes_from_seqrecord(model_organism_genes, records) # Uses records created earlier to get all the required model ortholog sequences. 

    taxa_sequences = select_genes_from_seqrecord(select_big_taxa_seq(model_sequences), model_sequences) # So we technically want only one gene to look at per taxa, right? Take the biggest.

    top_10 = select_biggest_k_seq(9, taxa_sequences) # Great now we take the 9 biggest to combine with our original


def test():
    gene = "WBGene00004963" # choose a gene
    gene_id = "6239_0:000672"
    output_phylip_file = 'output.phylip'
    output_newick_file = 'tree.nwk'

    orthologs = pull_model_organism_orthologs(gene) # pulling the model organisms that have ortholog genes.
    records = fasta_to_seqrecord(pull_OG_fasta('430340at2759')) # this seems to be pulling the ortholog group organisms in general, and putting their sequences in a dict for future ref. 
    model_organism_genes = get_model_organism_genes(orthologs)
    model_sequences = select_genes_from_seqrecord(model_organism_genes, records)
    taxa_sequences = select_genes_from_seqrecord(select_big_taxa_seq(model_sequences), model_sequences)
    top_10 = select_biggest_k_seq(10, taxa_sequences, gene_id, records)

    padded_sequences = pad_sequence_lengths(top_10)
    padded_and_removed_sequences = remove_stop_codons_in_multiple(padded_sequences)
    ## align = MultipleSeqAlignment(list(padded_and_removed_sequences.values()))

    align = format_ids_and_create_alignment(padded_and_removed_sequences)

    write_phylip_manual(align, output_phylip_file)
    
    tree = construct_phylo_tree(align)      
    write_newick_tree_with_header(tree, output_newick_file)

    results = run_codeml_positive_selection(output_newick_file, output_phylip_file)

    print(results)

    # Write in strict sequential PHYLIP format
    # with open('output.phylip', 'w') as f: # Use 'w' to overwrite during testing
    #     AlignIO.write(align, f, "phylip-sequential")
test()



