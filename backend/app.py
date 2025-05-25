import requests
import json
import os
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.Phylo.PAML import codeml

def cache_data(data, name):
    path = f"backend/cached/{name}"
    if type(data) == dict:
        print("Dumping json...")
        with open(path, "w") as file:
            json.dump(data, file)
    else:
        with open(path, 'a') as f:
            print(data, file=f)
    return


def pull_OG_fasta(ortholog_group):
    if os.path.exists(f'backend/cached/{ortholog_group}.fasta'):
        print("Data found in cache...")
        return open(f"backend/cached/{ortholog_group}.fasta").read()
    print("Pulling ortholog group data...")
    url = f"https://data.orthodb.org/v12/fasta"
    params = {
            'id': ortholog_group,
            'seqtype': 'cds'
            }
    response_text = requests.get(url, params=params).text
    cache_data(response_text,f"{ortholog_group}.fasta") #cache the response so I stop maybe accidentally crashing the db
    print("Ortholog group data pulled...")
    return response_text


def fasta_to_dict(response_text):
    response_text = response_text.split_lines()
    gene_index = {response_text[index].split(' ', 1)[0]: response_text[index+1] for index in range(0, len(response_text), 3)}
    return gene_index


def pull_model_organism_orthologs(gene):
    if os.path.exists(f"backend/cached/{gene}_search.json"):
        print("Data found in cache...")
        with open(f'backend/cached/{gene}_search.json', 'r') as file:
            response_json = json.load(file)
    else:
        gene_search = f"https://data.orthodb.org/v12/genesearch?query={gene}"
        response_json = requests.get(gene_search).json()
        cache_data(response_json,f"{gene}_search.json") #cache the response so I stop maybe accidentally crashing the db

    orthologs = response_json["orthologs_in_model_organisms"]
    ortholog_count = sum(sum(1 for found in organism["genes"]) for organism in orthologs)
    print(f"{ortholog_count} orthologs found in model organisms...")
    return orthologs

def fasta_to_seqrecord(response_text):
    fasta_handle = StringIO(response_text)
    records = {seq.id: seq for seq in SeqIO.parse(fasta_handle, "fasta")}

    print(f"Fasta data converted to seqrecord dict...")
    return records

def select_genes_from_seqrecord(genes, records):
    missing = []
    result = {}
    print("Paring down records for selected genes...")

    for gene in genes:
        if records.get(gene):
            result[gene] = records[gene]
        else:
            missing.append(gene)

    if missing:
        print(f"WARNING: Discrepancies found between records and selected genes...")
        for missed in missing:
            print(f"Gene id: {missed} not found in records...")

    # result = {gene: records[gene] for gene in genes}

    return result

def get_model_organism_genes(orthologs):
    return [gene_id["gene_id"]["param"] for organism in orthologs for gene_id in organism["genes"]]

def pad_sequence_lengths(records):
    big = max(len(record.seq) for record in records.values())
    print(f"Maximum sequence length is {big}...")

    for record in records.keys():
        records[record].seq = Seq(str(records[record].seq).ljust(big, '-'))
    return records

def select_big_taxa_seq(records, gene_id):
    result = {taxa: gene_id for taxa, gene_id in [key.split(":") for key in records.keys()]} ##initialize a dict to hold the biggest boyos with whatever default id
    for taxa_key in records.keys():
        taxa = taxa_key.split(":")[0]
        result[taxa] = max([result[taxa], taxa_key.split(":")[1]], key=lambda x: len(records[taxa + ":" + x].seq)) ##find the bigger of the two sequences and replace with that id belonging to that sequence and taxa
    print(f"{len(result)} unique taxa found...")
    return list(key + ":" + value for key, value in result.items())

def select_biggest_k_seq(k, records):
    #todo: add error checking to make sure k is smaller than or equal to records length
    top_k = sorted(list(records.values()), key=lambda x: len(x.seq))[-(k-1):] #don't forget to append the selected gene here later
    converted_dict = {seq.id: seq for seq in top_k}
    converted_dict[gene_id] = records[gene_id]
    return converted_dict

def remove_stop_codons_in_sequence(record, codon_stop_array = ["TAG", "TGA", "TAA"]): ## thanks to au_ndh at https://www.biostars.org/p/296261/ for this def
    tempRecordSeq = list(record.seq)
    for index in range(0, len(record.seq), 3):
            codon = record.seq[index:index+3]
            if codon in codon_stop_array:
                tempRecordSeq[index:index+3] = '-','-','-'
                print(f"removed a stop codon in {record.id}")
    record.seq = Seq("".join(tempRecordSeq))
    return record

def remove_stop_codons_in_multiple(records):
    result = {key:remove_stop_codons_in_sequence(record) for key, record in records.items()}
    return result

def construct_phylo_tree(alignment, rooted=False):
    print("Constructing phylogenetic tree...")
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix) if rooted else constructor.nj(distance_matrix)
    return tree

def convert_organism_id_to_names(records):
    for gene_id in records.keys():
        print(records[gene_id])
        organism_desc_str = records[gene_id].description
        organism_name = organism_desc_str[organism_desc_str.index("{"):].split(",")[4].split(":")[1][1:-1]
        print(organism_name)
        current_id = records[gene_id].id
        print(current_id)
        records[gene_id].id = organism_name
    return records

def run_CODEML_analysis(tree_path, seqfile_path):
    for record in records.values():
        if record[""]
    cml = codeml.Codeml()
    cml.alignment = seqfile_path
    cml.tree = tree_path
    cml.out_file = "results.out"
    cml.working_dir = "./scratch"
    cml.set_options(
        seqtype=1,
        verbose=0,
        noisy=0,
        RateAncestor=0,
        model=0,
        NSsites=[0, 1, 2],
        CodonFreq=2,
        cleandata=1,
        fix_alpha=1,
        kappa=4.54006,
    )
    results = cml.run()
    ns_sites = results.get("NSsites")
    m0 = ns_sites.get(0)
    m0_params = m0.get("parameters")
    print(m0_params.get("omega"))

def main():
    gene = "WBGene00004963"
    gene_id = "6239_0:000672"
    print("Starting up")
    orthologs = pull_model_organism_orthologs(gene)
    records = fasta_to_seqrecord(pull_OG_fasta('430340at2759'))
    print(list(records.items())[0][1])
    model_organism_genes = get_model_organism_genes(orthologs)
    model_sequences = select_genes_from_seqrecord(model_organism_genes, records)
    taxa_sequences = select_genes_from_seqrecord(select_big_taxa_seq(model_sequences), model_sequences)
    top_10 = select_biggest_k_seq(10, taxa_sequences, gene_id)
    padded_sequences = pad_sequence_lengths(top_10)
    padded_and_removed_sequences = remove_stop_codons_in_multiple(padded_sequences)
    padded_removed_and_stripped_sequences = convert_organism_id_to_names(padded_and_removed_sequences)
    align = MultipleSeqAlignment(list(padded_and_removed_sequences.values()))
    print(align)
    with open('output.phylip', 'a') as f:
        print(format(align, "phylip"), file=f)

    tree = construct_phylo_tree(align)

    Phylo.write(tree, 'tree.nex', "nexus")
    # run_CODEML_analysis('tree.nwk', 'output.phylip')
    # freq = {}
    # for gene in taxa_sequences.values():
    #     length = len(gene.seq)
    #     freq[length] = freq.get(length, 0) + 1
    #
    # print(len(records["6239_0:000672"].seq))
    # print((sorted(freq.items(), key=lambda x: x[1])[-1]))
