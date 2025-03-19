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

def select_big_taxa_seq(records):
    result = {taxa: gene_id for taxa, gene_id in [key.split(":") for key in records.keys()]} ##initialize a dict to hold the biggest boyos with whatever default id
    for taxa_key in records.keys():
        taxa = taxa_key.split(":")[0]
        result[taxa] = max([result[taxa], taxa_key.split(":")[1]], key=lambda x: len(records[taxa + ":" + x].seq)) ##find the bigger of the two sequences and replace with that id belonging to that sequence and taxa
    print(f"{len(result)} unique taxa found...")
    return list(key + ":" + value for key, value in result.items())

def select_biggest_k_seq(k, records):
    #todo: add error checking to make sure k is smaller than or equal to records length
    top_k = sorted(list(records.values()), key=lambda x: len(x.seq))[-k:] #don't forget to append the selected gene here later
    converted_dict = {seq.id: seq for seq in top_k}
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

gene = "WBGene00004963"
orthologs = pull_model_organism_orthologs(gene)
records = fasta_to_seqrecord(pull_OG_fasta('430340at2759'))
print(list(records.items())[0][1])
model_organism_genes = get_model_organism_genes(orthologs)
model_sequences = select_genes_from_seqrecord(model_organism_genes, records)
taxa_sequences = select_genes_from_seqrecord(select_big_taxa_seq(model_sequences), model_sequences)
top_10 = select_biggest_k_seq(9, taxa_sequences)
padded_sequences = pad_sequence_lengths(top_10)
padded_and_removed_sequences = remove_stop_codons_in_multiple(padded_sequences)
align = MultipleSeqAlignment(list(padded_and_removed_sequences.values()))
print(align)
with open('output.txt', 'a') as f:
    print(format(align, "phylip"), file=f)
# freq = {}
# for gene in taxa_sequences.values():
#     length = len(gene.seq)
#     freq[length] = freq.get(length, 0) + 1
#
# print(len(records["6239_0:000672"].seq))
# print((sorted(freq.items(), key=lambda x: x[1])[-1]))




#
#
# for seqrecord in align:
#     seqrecord.id = seqrecord.id[-10:]
#     
# print(align)
#
#
