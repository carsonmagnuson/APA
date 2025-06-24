from pipe_fun import ortholog_functions


test_ncbi_ids = ["173042", "173402", "3039"] # C. elegans spe-39, C. elegans lin-39, Human HBA1
test_levels = ["6231", "6231", "40674"]    # Nematoda, Nematoda, Mammalia

def test_compile_orthologs():
    ortholog_compilation = ortholog_functions.compile_orthologs(test_ncbi_ids[0], test_levels[0])
    print(len(ortholog_compilation))
    assert len(ortholog_compilation) > 10
