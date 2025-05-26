import os
import shutil
import logging
from typing import Optional, Dict, Any
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.Phylo.PAML import codeml # Ensure PAML is installed and codeml executable is in PATH

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
    if not alignment or len(alignment) < 2: # Need at least 2 sequences for a tree
        logger.warning(f"Alignment is too small (size {len(alignment)}) to construct a tree. Need at least 2 sequences.")
        return None

    logger.info(f"Constructing phylogenetic tree using {'UPGMA (rooted)' if rooted else 'NJ (unrooted)'} method...")
    try:
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)
    except Exception as e:
        logger.error(f"Error calculating distance matrix: {e}", exc_info=True)
        return None

    constructor = DistanceTreeConstructor(calculator) # Can pass calculator directly
    
    try:
        if rooted:
            tree = constructor.upgma(distance_matrix)
        else:
            tree = constructor.nj(distance_matrix)
        logger.info("Phylogenetic tree constructed successfully.")
        return tree
    except Exception as e: # Catch errors during tree construction (e.g., invalid distance matrix)
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

    Specifically, it configures and runs CODEML for models M0 (one-ratio), 
    M1a (neutral), and M2a (selection) which are commonly used for
    Likelihood Ratio Tests (LRTs) to detect positive selection at specific codon sites.

    Args:
        tree_filepath: Path to the Newick tree file (must have PAML-style header).
        seqfile_filepath: Path to the PHYLIP sequence alignment file.
        output_dir: Directory where CODEML will run and store its output files.
                    This directory will be created if it doesn't exist.
                    Defaults to "paml_output".
        results_filename: Name of the main CODEML output file (e.g., "results.out").
                          Defaults to "results.out".
        cleanup_working_dir: If True, removes the `output_dir` after successful parsing
                             of results. Defaults to True.

    Returns:
        A dictionary containing parsed results for models M0, M1a, and M2a,
        including lnL values and key parameter estimates if the run is successful.
        Returns None if the CODEML run or result parsing fails.
        The structure is: {'M0': {'lnL': ..., 'parameters': ...}, ...}
    """
    cml = codeml.Codeml()

    # Ensure absolute paths for PAML, as it can be sensitive
    abs_seqfile_filepath = os.path.abspath(seqfile_filepath)
    abs_tree_filepath = os.path.abspath(tree_filepath)
    abs_output_dir = os.path.abspath(output_dir)

    if not os.path.exists(abs_seqfile_filepath):
        logger.error(f"Sequence alignment file not found: {abs_seqfile_filepath}")
        return None
    if not os.path.exists(abs_tree_filepath):
        logger.error(f"Tree file not found: {abs_tree_filepath}")
        return None

    cml.alignment = abs_seqfile_filepath
    cml.tree = abs_tree_filepath
    cml.out_file = results_filename  # This will be inside the working_dir

    if not os.path.exists(abs_output_dir):
        try:
            os.makedirs(abs_output_dir)
            logger.info(f"Created CODEML working directory: {abs_output_dir}")
        except OSError as e:
            logger.error(f"Error creating CODEML output directory {abs_output_dir}: {e}", exc_info=True)
            return None
    cml.working_dir = abs_output_dir

    logger.info(f"Running CODEML in directory: {abs_output_dir}")
    logger.info(f"Alignment file: {cml.alignment}")
    logger.info(f"Tree file: {cml.tree}")

    try:
        cml.set_options(
            noisy=0,          # 0-9: Amount of screen output. 0 or 1 is fine.
            verbose=1,        # 0: concise; 1: detailed output in results.out.
            runmode=0,        # 0: user tree. We provide the tree.
            seqtype=1,        # 1: codon sequences. Essential.
            CodonFreq=2,      # 2: F3x4 model (codon freqs from nt freqs at 3 positions).
            model=0,          # General setting used with NSsites.
            NSsites=[0, 1, 2],# M0 (one-ratio), M1a (neutral), M2a (selection).
            icode=0,          # 0: universal genetic code.
            fix_kappa=0,      # 0: estimate kappa (ts/tv ratio).
            kappa=2.0,        # Initial guess for kappa.
            fix_omega=0,      # 0: estimate omega. Essential for these models.
            omega=0.5,        # Initial guess for omega.
            fix_alpha=1,      # For M0, M1a, M2a, alpha is typically fixed at 0.
            alpha=0.0,        # No continuous gamma heterogeneity.
            cleandata=1,      # 1: Remove sites with alignment gaps/ambiguities.
            Small_Diff=1e-7,
            getSE=0           # 0: Don't calculate Standard Errors (faster).
        )
    except ValueError as e:
        logger.error(f"Error setting PAML options: {e}", exc_info=True)
        return None

    parsed_results: Dict[str, Dict[str, Any]] = {}
    try:
        logger.info("Starting CODEML run...")
        # Ensure 'codeml' executable is in system PATH or PAML_PATH env variable is set.
        results = cml.run(verbose=True) # verbose=True logs Biopython's call to codeml
        logger.info("CODEML run finished.")
        
        if "NSsites" in results:
            for model_key, model_results_dict in results["NSsites"].items():
                # model_results_dict is already a dictionary from Biopython's parsing
                parsed_results[f"M{model_key}"] = {
                    "lnL": model_results_dict.get("lnL"),
                    "parameters": model_results_dict.get("parameters", {}) # Ensure parameters is a dict
                }
            logger.info("Successfully parsed results from CODEML output.")
        else:
            logger.warning("'NSsites' key not found in PAML results. Parsing might be incomplete.")
            logger.debug(f"Full PAML results object: {results}")

    except Exception as e:
        logger.error(f"An error occurred during CODEML execution or parsing: {e}", exc_info=True)
        logger.error(f"Check files in the working directory for details: {abs_output_dir}")
        # It's often useful to inspect results.out manually if an error occurs.
        return None # Return None or partial results depending on desired behavior
    finally: # Ensure cleanup happens if specified, even if parsing failed after a successful run
        if cleanup_working_dir and os.path.exists(abs_output_dir):
            try:
                shutil.rmtree(abs_output_dir)
                logger.info(f"Cleaned up CODEML working directory: {abs_output_dir}")
            except OSError as e:
                logger.error(f"Error cleaning up working directory {abs_output_dir}: {e}", exc_info=True)
        
    return parsed_results if parsed_results else None


if __name__ == '__main__':
    logger.info("--- Phylogenetic Analysis Module Tests ---")
    # Note: Testing construct_phylo_tree and run_codeml_positive_selection effectively
    # requires valid alignment files, tree files, and potentially a PAML installation.
    # These are more suited for integration tests or require mock objects/files.
    
    # Example (conceptual) of how you might test construct_phylo_tree if you had a mock alignment
    # from Bio.Seq import Seq
    # from Bio.SeqRecord import SeqRecord
    # from Bio.Align import MultipleSeqAlignment
    # rec1 = SeqRecord(Seq("ATGCGT"), id="seq1")
    # rec2 = SeqRecord(Seq("ATGCCT"), id="seq2")
    # rec3 = SeqRecord(Seq("ATGCGA"), id="seq3")
    # mock_alignment = MultipleSeqAlignment([rec1, rec2, rec3])
    # if mock_alignment:
    #     tree = construct_phylo_tree(mock_alignment)
    #     if tree:
    #         logger.info("Mock tree constructed:")
    #         Phylo.draw_ascii(tree) # Simple way to visualize
    #     else:
    #         logger.error("Mock tree construction failed.")
    
    logger.info("Minimal tests for phylogenetic_analysis.py. Full tests require data files and PAML.")
    logger.info("--- Phylogenetic Analysis Module Tests Finished ---")

