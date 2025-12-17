import json
import os
import logging
from typing import Any, Dict, Union, Optional
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment

# --- Configure Logging ---
# Basic configuration for standalone testing of this module.
# In a larger application, configure logging in the main entry point.
if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(name)s - %(funcName)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

logger = logging.getLogger(__name__)

# --- Constants ---
CACHE_DIR = "cached" # Default cache directory

def cache_data(data: Any, filename: str, directory: str = CACHE_DIR) -> None:
    """
    Caches data to a file in the specified directory.

    If the directory doesn't exist, it will be created.
    If the data is a dictionary, it's saved as a JSON file.
    Otherwise, it's appended to the file as a string.

    Args:
        data: The data to cache (can be a dictionary or any object convertible to string).
        filename: The name of the file to save the data to (e.g., "gene_data.json").
        directory: The directory where the cache file will be stored.
                   Defaults to CACHE_DIR.
    """
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
            logger.info(f"Created cache directory: {directory}")
        except OSError as e:
            logger.error(f"Error creating cache directory {directory}: {e}", exc_info=True)
            return

    path = os.path.join(directory, filename)

    try:
        if isinstance(data, dict):
            logger.info(f"Dumping JSON data to: {path}")
            with open(path, "w") as file:
                json.dump(data, file, indent=4)
        else:
            logger.info(f"Appending string data to: {path}") # Original was 'a', consider 'w'
            with open(path, 'w') as f: # Changed to 'w' for typical caching behavior
                f.write(str(data)) # Removed newline for direct string saving
                if not isinstance(data, str) or not data.endswith("\n"):
                     f.write("\n") # Add newline if not already present or not a string that ends with it
    except IOError as e:
        logger.error(f"Error writing to cache file {path}: {e}", exc_info=True)

def write_phylip_manual(alignment: Optional[MultipleSeqAlignment], filepath: str) -> None:
    """
    Manually writes the alignment to a file in PHYLIP sequential format.

    This format is often required for compatibility with phylogenetic software like PAML/CODEML.
    It ensures correct spacing for sequence IDs (padded to 10 characters).

    Args:
        alignment: The Bio.Align.MultipleSeqAlignment object to write.
                   It's assumed that record IDs are already cleaned/formatted.
        filepath: The path to the output PHYLIP file.
    """
    if not alignment or len(alignment) == 0:
        logger.warning(f"Alignment is empty or None. Cannot write PHYLIP file to {filepath}.")
        # Optionally, write an empty valid PHYLIP header if downstream tools expect it
        # try:
        #     with open(filepath, 'w') as f_empty:
        #         f_empty.write(" 0 0\n")
        # except IOError as e:
        #     logger.error(f"Error writing empty PHYLIP file to {filepath}: {e}", exc_info=True)
        return

    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()

    try:
        with open(filepath, 'w') as f:
            # Write the header line: e.g., " 10 12531"
            f.write(f" {num_sequences} {alignment_length}\n")

            for record in alignment:
                # IDs should ideally be cleaned by format_ids_and_create_alignment
                # to be <=10 chars and PAML-compatible.
                # Pad with spaces on the right to ensure it's at least 10 characters.
                # PAML might truncate IDs longer than 10 chars, so ensure they are prepared.
                name = record.id.ljust(10)[:10] # Ensure exactly 10 characters
                
                sequence_string = str(record.seq)

                # Write the 10-character name, two spaces, then the sequence.
                f.write(f"{name}  {sequence_string}\n")
            
        logger.info(f"PHYLIP file '{filepath}' manually written in sequential format.")
    except IOError as e:
        logger.error(f"Error writing PHYLIP file to {filepath}: {e}", exc_info=True)


def write_newick_tree_with_header(tree: Optional[Phylo.BaseTree.Tree], filepath: str) -> None:
    """
    Writes the phylogenetic tree to a file in Newick format,
    prefixed with a PAML-style header (num_taxa num_trees).

    Args:
        tree: The Bio.Phylo.BaseTree.Tree object to write.
        filepath: The path to the output Newick file.
    """
    if not tree:
        logger.warning(f"Tree object is None. Cannot write Newick file to {filepath}.")
        return

    try:
        num_taxa = len(tree.get_terminals())  # Get the number of leaves (taxa)
        if num_taxa == 0:
            logger.warning(f"Tree has no taxa (terminals). Cannot write valid Newick file with header to {filepath}.")
            return
    except Exception as e: # Catching potential issues with tree object
        logger.error(f"Could not determine number of taxa from tree: {e}. Cannot write Newick file to {filepath}.", exc_info=True)
        return

    try:
        with open(filepath, 'w') as f_tree:
            f_tree.write(f" {num_taxa} 1\n") # PAML often expects a leading space here
            Phylo.write(tree, f_tree, "newick")
            # Phylo.write does not add a final newline by default, so we add one.
            f_tree.write("\n") 
            
        logger.info(f"Newick tree with header written to '{filepath}'.")
    except IOError as e:
        logger.error(f"Error writing Newick tree file to {filepath}: {e}", exc_info=True)
    except Exception as e: # Catch-all for other Phylo.write issues
        logger.error(f"An unexpected error occurred while writing Newick tree: {e}", exc_info=True)

if __name__ == '__main__':
    logger.info("--- File Utils Module Tests ---")
    
    # Test cache_data
    logger.info("\n--- Testing cache_data ---")
    test_dict_data = {"key": "value", "number": 123}
    cache_data(test_dict_data, "test_cache.json")
    test_str_data = "This is a test string for caching."
    cache_data(test_str_data, "test_cache.txt")
    logger.info(f"Cache tests complete. Check the '{CACHE_DIR}' directory.")

    # Note: Testing write_phylip_manual and write_newick_tree_with_header
    # would require creating mock Bio.Align.MultipleSeqAlignment and 
    # Bio.Phylo.BaseTree.Tree objects, which is more involved for a simple test block.
    # These are typically tested as part of the integration in the main application.
    logger.info("--- File Utils Module Tests Finished ---")

