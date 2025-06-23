from pipe_fun import ortholog_functions

def test_compile_orthologs():
    ortholog_compilation = ortholog_functions.compile_orthologs()
    assert ortholog_compilation
