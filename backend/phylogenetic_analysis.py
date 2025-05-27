import os
import shutil
import logging
from typing import Optional, Dict, Any
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.Phylo.PAML import codeml # Ensure PAML is installed and codeml executable is in PATH
# Removed problematic PamlError and _paml imports

# --- Configure Logging ---
if __name__ == '__main__': # Basic configuration for standalone testing
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(name)s - %(funcName)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
logger = logging.getLogger(__name__)

# --- Phylogenetic Analysis Functions ---

def construct_phylo_tree(alignment: MultipleSeqAlignment, rooted: bool = False) -> Optional[Phylo.BaseTree.Tree]:
    """
    Constructs a phylogenetic tree from a multiple sequence alignment.

    Uses the DistanceCalculator with the 'identity' model and then either
    UPGMA (if rooted=True) or Neighbor Joining (NJ) (if rooted=False)
    from DistanceTreeConstructor to build the tree.

    Args:
        alignment: A Bio.Align.MultipleSeqAlignment object.
        rooted: If True, constructs a UPGMA tree (rooted). 
                If False (default), constructs an NJ tree (unrooted).

    Returns:
        A Bio.Phylo.BaseTree.Tree object if successful, or None if an error occurs
        (e.g., alignment is too small or unsuitable).
    """
    if not alignment or len(alignment) < 2: 
        logger.warning(f"Alignment is too small (size {len(alignment)}) to construct a tree. Need at least 2 sequences.")
        return None

    logger.info(f"Constructing phylogenetic tree using {'UPGMA (rooted)' if rooted else 'NJ (unrooted)'} method...")
    try:
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)
    except ValueError as ve: 
        logger.error(f"Error calculating distance matrix, possibly due to alignment issues: {ve}", exc_info=True)
        return None
    except Exception as e:
        logger.error(f"Unexpected error calculating distance matrix: {e}", exc_info=True)
        return None

    constructor = DistanceTreeConstructor(calculator)
    
    try:
        if rooted:
            tree = constructor.upgma(distance_matrix)
        else:
            tree = constructor.nj(distance_matrix)
        logger.info("Phylogenetic tree constructed successfully.")
        return tree
    except Exception as e: 
        logger.error(f"Error constructing phylogenetic tree: {e}", exc_info=True)
        return None


def run_codeml_positive_selection(
    tree_filepath: str, 
    seqfile_filepath: str, 
    output_dir: str = "paml_output", 
    results_filename: str = "results.out",
    cleanup_working_dir: bool = True
) -> Optional[Dict[str, Any]]:
    """
    Runs CODEML from the PAML package to test for positive selection using site models.

    Args:
        tree_filepath: Path to the Newick tree file (must have PAML-style header).
        seqfile_filepath: Path to the PHYLIP sequence alignment file.
        output_dir: Directory where CODEML will run and store its output files.
        results_filename: Name of the main CODEML output file.
        cleanup_working_dir: If True, removes the `output_dir` after successful parsing.

    Returns:
        A dictionary containing parsed results or error information.
    """
    cml = codeml.Codeml()

    abs_seqfile_filepath = os.path.abspath(seqfile_filepath)
    abs_tree_filepath = os.path.abspath(tree_filepath)
    abs_output_dir = os.path.abspath(output_dir)
    results_out_full_path = os.path.join(abs_output_dir, results_filename)

    if not os.path.exists(abs_seqfile_filepath):
        logger.error(f"Sequence alignment file not found: {abs_seqfile_filepath}")
        return {"status": "SETUP_ERROR", "message": f"Sequence file not found: {abs_seqfile_filepath}"}
    if not os.path.exists(abs_tree_filepath):
        logger.error(f"Tree file not found: {abs_tree_filepath}")
        return {"status": "SETUP_ERROR", "message": f"Tree file not found: {abs_tree_filepath}"}

    cml.alignment = abs_seqfile_filepath
    cml.tree = abs_tree_filepath
    cml.out_file = results_filename 
    cml.working_dir = abs_output_dir

    if not os.path.exists(abs_output_dir):
        try:
            os.makedirs(abs_output_dir)
            logger.info(f"Created CODEML working directory: {abs_output_dir}")
        except OSError as e:
            logger.error(f"Error creating CODEML output directory {abs_output_dir}: {e}", exc_info=True)
            return {"status": "SETUP_ERROR", "message": f"Error creating output directory: {str(e)}"}

    logger.info(f"Running CODEML in directory: {abs_output_dir}")
    logger.info(f"Alignment file: {cml.alignment}")
    logger.info(f"Tree file: {cml.tree}")

    try:
        cml.set_options(
            noisy=0, verbose=1, runmode=0, seqtype=1, CodonFreq=2,
            model=0, NSsites=[0, 1, 2], icode=0, fix_kappa=0, kappa=2.0,
            fix_omega=0, omega=0.5, fix_alpha=1, alpha=0.0,
            cleandata=1, Small_Diff=1e-7, getSE=0
        )
    except ValueError as e:
        logger.error(f"Error setting PAML options: {e}", exc_info=True)
        return {"status": "SETUP_ERROR", "message": f"Error setting PAML options: {str(e)}"}

    parsed_results: Dict[str, Any] = {}
    
    try:
        logger.info("Starting CODEML run...")
        results_from_biopython = cml.run(verbose=True)
        logger.info("CODEML run command finished. Biopython parsing attempted.")

        if results_from_biopython and "NSsites" in results_from_biopython:
            for model_key, model_data in results_from_biopython["NSsites"].items():
                parsed_results[f"M{model_key}"] = {
                    "lnL": model_data.get("lnL"),
                    "parameters": model_data.get("parameters", {})
                }
            logger.info("Successfully parsed results from CODEML output.")
            if cleanup_working_dir:
                try:
                    shutil.rmtree(abs_output_dir)
                    logger.info(f"Cleaned up CODEML working directory: {abs_output_dir}")
                except OSError as e_clean:
                    logger.error(f"Error cleaning up working directory {abs_output_dir}: {e_clean}", exc_info=True)
            return parsed_results
        else: # CODEML ran, Biopython parsing ran, but 'NSsites' key is missing
            logger.warning("'NSsites' key not found in PAML results from Biopython. Parsing might be incomplete or CODEML output was not as expected.")
            logger.debug(f"Full PAML results object from Biopython: {results_from_biopython}")
            return {
                "status": "CODEML_RAN_PARSING_INCOMPLETE",
                "message": "CODEML completed, but the 'NSsites' section was not found in the parsed output. Check results.out.",
                "results_file_path": results_out_full_path,
                "raw_parsed_output": results_from_biopython if isinstance(results_from_biopython, dict) else str(results_from_biopython)
            }
    except FileNotFoundError as fnfe: # For codeml executable not found
        logger.error(f"CODEML executable not found or path issue: {fnfe}", exc_info=True)
        return {
            "status": "CODEML_SETUP_ERROR",
            "message": f"CODEML executable error: {str(fnfe)}. Ensure PAML is installed and codeml is in your PATH.",
            "error_details": str(fnfe)
        }
    except UnboundLocalError as ule: # Catch the specific parsing error
        logger.error(f"UnboundLocalError during Biopython's PAML result parsing (likely within cml.run): {ule}", exc_info=True)
        if os.path.exists(results_out_full_path):
            logger.info(f"results.out file exists at: {results_out_full_path}. Biopython parsing failed with UnboundLocalError.")
            return {
                "status": "CODEML_RAN_PARSING_FAILED",
                "message": "CODEML completed, but Biopython failed to parse results.out due to an UnboundLocalError. The raw output file is available.",
                "results_file_path": results_out_full_path,
                "error_details": str(ule)
            }
        else:
            logger.error(f"CODEML run might have failed before producing {results_filename}. UnboundLocalError occurred, and no output file found.")
            return {
                "status": "CODEML_RUN_OR_PARSING_FAILED_UNBOUNDLOCAL",
                "message": "CODEML run or parsing failed with UnboundLocalError, and results.out was not found.",
                "error_details": str(ule)
            }
    except Exception as e: # Catch other errors from cml.run() or our logic
        logger.error(f"An unexpected error occurred during CODEML execution or result handling: {e}", exc_info=True)
        error_status = {
            "status": "CODEML_UNEXPECTED_ERROR",
            "message": f"An unexpected error occurred: {str(e)}",
            "error_details": str(e)
        }
        if os.path.exists(results_out_full_path):
             error_status["results_file_path"] = results_out_full_path
             logger.info(f"results.out file exists at: {results_out_full_path}, despite the error.")
        else:
             logger.warning(f"results.out file not found at {results_out_full_path} after error.")
        return error_status

if __name__ == '__main__':
    logger.info("--- Phylogenetic Analysis Module Tests ---")
    test_output_dir = "paml_output_test"
    if not os.path.exists(test_output_dir):
        os.makedirs(test_output_dir)

    dummy_aln_file = os.path.join(test_output_dir, "dummy_aln.phylip")
    dummy_tree_file = os.path.join(test_output_dir, "dummy_tree.nwk")

    with open(dummy_aln_file, "w") as f:
        f.write(" 2 30\n") 
        f.write("seq1      ATGATGATGATGATGATGATGATGATGTGA\n")
        f.write("seq2      ATGATGATGATGATGATGATGATGATGTGA\n")

    with open(dummy_tree_file, "w") as f:
        f.write(" 2 1\n") 
        f.write("(seq1:0.1,seq2:0.1);\n")

    logger.info(f"Attempting CODEML run with dummy files in {test_output_dir}...")
    results = run_codeml_positive_selection(
        tree_filepath=dummy_tree_file,
        seqfile_filepath=dummy_aln_file,
        output_dir=test_output_dir,
        results_filename="test_results.out",
        cleanup_working_dir=False 
    )
    logger.info(f"Test CODEML run results: {results}")
    logger.info("--- Phylogenetic Analysis Module Tests Finished ---")

