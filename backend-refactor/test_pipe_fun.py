
from pipe_fun import ortholog_functions, external_tools

test_ncbi_ids = ["173042", "173402", "3039"] # C. elegans spe-39, C. elegans lin-39, Human HBA1
test_levels = ["6231", "6231", "40674"]    # Nematoda, Nematoda, Mammalia
test_ortholog_group_ids = ["29865at6231"]

def test_select_ortholog_group():
    assert ortholog_functions.select_ortholog_group(test_ncbi_ids[0], test_levels[0]) == test_ortholog_group_ids[0]

def test_compile_orthologs():
    ortholog_compilation = ortholog_functions.compile_orthologs(test_ortholog_group_ids[0])
    print(len(ortholog_compilation))
    assert len(ortholog_compilation) > 0

def test_select_orthologs(k = 8):
    selection_path = ortholog_functions.select_orthologs(ortholog_functions.compile_orthologs(test_ortholog_group_ids[0]), k, test_ortholog_group_ids[0])
    assert int(selection_path.split('/')[1][0]) == k # Check that the file has the correct number of selected orthologs

def test_run_muscle(): 
    assert external_tools.run_muscle("29865at6231_run/8_selected_orthologs.fasta") == "29865at6231_run/8_selected_orthologs.phylip"

def test_run_iqtree(): 
    assert external_tools.run_iqtree("29865at6231_run/8_selected_orthologs.phylip") == "29865at6231_run/8_selected_orthologs.treefile"
