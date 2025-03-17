import itertools
import requests
import json

# gene = "WBGene00004963"
#
# search = f"https://data.orthodb.org/v12/search?query={gene}"
#
# gene_search = f"https://data.orthodb.org/v12/genesearch?query={gene}"
#
# response_json = requests.get(gene_search).json()
# response_json_generic = requests.get(search).json()
# ortholog_groups = response_json_generic["data"]

def pull_OG_fasta(ortholog_group):
    url = f"https://data.orthodb.org/v12/fasta"
    params = {
            'id': ortholog_group,
            'seqtype': 'cds'
            }
    response_text = requests.get(url, params=params).text.splitlines()
    print("data pulled...processing...")
    print(response_text[0][0:10])
    print(response_text[1][0:10])
    print(response_text[2][0:10])
    print(response_text[3][0:10])
    gene_index = {response_text[index].split(' ', 1)[0]: response_text[index+1] for index in range(0, len(response_text), 3)}
    return gene_index
    
    
# fasta_data = list(pull_OG_fasta(group) for group in ortholog_groups)

# print(len(fasta_data[0]))

test = pull_OG_fasta('430340at2759')

print(len(test.keys()))
print(list(test.items())[1])






# orthologs = response_json["orthologs_in_model_organisms"]
#
# ortholog_count = sum(sum(1 for found in organism["genes"]) for organism in orthologs)
#
# ortholog_groups = list(set(ortholog["lca"]["lca_cluster_id"] for ortholog in orthologs))
#
# print(ortholog_count)
# print(response_json["gene"].keys())
# print(orthologs[0].keys())
# print(orthologs[0]["lca"])
# print(ortholog_groups)

