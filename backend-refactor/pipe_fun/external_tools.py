import subprocess, os, shutil


def run_muscle(selection_path: str) -> str:
    """
    Runs MUSCLE tool on a fasta file of orthologs.

    Args:
        selection_path: What is the path of the fasta file with the selected orthologs?

    Returns:
        A path to a phylip file.

    """

    # STEP 1: Determine output path and design a MUSCLE command to execute.
    output_path = f"{selection_path.split('.')[0]}.phylip"
    muscle_command = [
        "muscle",
        "-in", selection_path,
        "-out", output_path,
        "-phyi" # Output in PHYLIP format (phyi for interleaved, phys for sequential)
    ]

    # STEP 2: Execute MUSCLE command and return path to phylip file.
    result = subprocess.run(muscle_command, check=True, capture_output=True, text=True)
    print(result)
    return output_path

def run_iqtree(phylip_path: str) -> str:
    """
    Runs IQTREE analysis on a phylip file of selected, aligned orthologs.

    Args:
        phylip_path: What is the path of the phylip file to be analyzed?

    Returns:
        A path to a tree.newick file.

    """
    # STEP 1: Determine output path and design a IQTREE command to execute.
    output_path = f"{phylip_path.split('.')[0]}"
    iqtree_command = [
        "iqtree3",
        "-s", phylip_path,
        "-m", "GTR+G", # Specifying to use General Time Reversible model with Gamma-dist rates
        "-nt", "auto", # Use all available CPU cores
        "-pre", output_path
    ]

    # STEP 2: Execute IQTREE command and return path to phylip file.
    result = subprocess.run(iqtree_command, check=True, capture_output=True, text=True)
    return f"{output_path}.treefile"


def run_codeml(
        phylip_path: str,
        tree_path: str,
        ) -> str:
    """
    Runs CODEML analysis on a phylip and tree files of selected, aligned orthologs.

    Args:
        phylip_path: What is the path of the phylip file to be analyzed?
        tree_path: What is the path of the tree file to be analyzed?

    Returns:
        A path to a results.out file.

    """

    # STEP 1: Check for/create a work directory for CODEML, copy phylip and tree files inside.
    run_path_name = phylip_path.split('.')[0] 
    output_directory = f"{run_path_name}"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    shutil.copy(phylip_path, f"{output_directory}/in.phylip")
    shutil.copy(tree_path, f"{output_directory}/in.treefile")

    # STEP 2: Create a control file for CODEML.
    results_name = f"results_M0.mlc"
    control_content = f"""
        seqfile = in.phylip
        treefile = in.treefile
        outfile = {results_name}

          noisy = 3
        verbose = 1
        runmode = 0
        seqtype = 1
      CodonFreq = 2
          model = 0
        NSsites = 0
          icode = 0
      fix_kappa = 0
          kappa = 2
      fix_omega = 0
          omega = 0.4
    """

    with open(f"{output_directory}/codeml.ctl", "w") as ctl_file:
        ctl_file.write(control_content)

    # STEP 3: Run CODEML analysis and return analysis file path.
    # result = subprocess.run("codeml", check=True, capture_output=True, text=True, cwd=output_directory)
    return f"{output_directory}/{results_name}"


if __name__ == "__main__":
    # print(run_muscle("29865at6231_run/8_selected_orthologs.fasta"))
    # print(run_iqtree("29865at6231_run/8_selected_orthologs.phylip"))
    print(run_codeml("29865at6231_run/8_selected_orthologs.phylip", "29865at6231_run/8_selected_orthologs.treefile"))

