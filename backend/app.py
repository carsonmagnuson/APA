import requests
import json

gene = "WBGene00004963"

gene_search = f"https://data.orthodb.org/v12/genesearch?query={gene}"

response = requests.get(gene_search).json()

print(json.dumps(response, indent=4))



