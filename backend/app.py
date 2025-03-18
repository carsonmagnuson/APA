import requests
import json
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
# from Bio.Align.Applications import MuscleCommandline

def pull_OG_fasta(ortholog_group):
    url = f"https://data.orthodb.org/v12/fasta"
    params = {
            'id': ortholog_group,
            'seqtype': 'cds'
            }
    response_text = requests.get(url, params=params).text
    print("Ortholog group data pulled...")
    return response_text


def fasta_to_dict(response_text):
    response_text = response_text.split_lines()
    gene_index = {response_text[index].split(' ', 1)[0]: response_text[index+1] for index in range(0, len(response_text), 3)}
    return gene_index


def pull_model_organism_orthologs(gene):
    gene_search = f"https://data.orthodb.org/v12/genesearch?query={gene}"
    response_json = requests.get(gene_search).json()
    orthologs = response_json["orthologs_in_model_organisms"]

    ortholog_count = sum(sum(1 for found in organism["genes"]) for organism in orthologs)
    print(f"{ortholog_count} orthologs found in model organisms...")
    return orthologs

def fasta_to_seqrecord(response_text):
    fasta_handle = StringIO(response_text)
    records = {seq.id: seq for seq in SeqIO.parse(fasta_handle, "fasta")}

    print(f"Fasta data converted to seqrecord dict...")
    return records

def gene_list_selection_from_seqrecord(genes, records):
    missing = []
    result = {}
    print("Paring down records for seleccted genes...")

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

def format_sequence_lengths(records):
    big = max(len(record.seq) for record in records.values())
    print(f"Maximum sequence length is {big}...")

    for record in records.keys():
        records[record].seq = Seq(str(records[record].seq).ljust(big, '-'))
    return records



gene = "WBGene00004963"
orthologs = pull_model_organism_orthologs(gene)
records = fasta_to_seqrecord(pull_OG_fasta('430340at2759'))
print(list(records.items())[0][1])
genes = gene_list_selection_from_seqrecord((get_model_organism_genes(orthologs)), records)

formatted_genes = format_sequence_lengths(genes)

print(list(genes.items())[0][1])


freq = {}
for gene in formatted_genes.values():
    length = len(gene.seq)
    freq[length] = freq.get(length, 0) + 1

print((sorted(freq.items(), key=lambda x: x[1])[-1]))


align = MultipleSeqAlignment(list(formatted_genes.values())[0:4])
print(align)

# print(format(align, "phylip"))

