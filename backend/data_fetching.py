
import requests
import json
import os
import logging # Import the logging module
from typing import Dict, Any, List, Optional

# --- Configure Logging ---
# This basic configuration sets up logging to the console.
# You can customize this further (e.g., log to a file, different formats).
logging.basicConfig(
    level=logging.INFO,  # Set the default logging level
    format='%(asctime)s - %(levelname)s - %(name)s - %(funcName)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# Create a logger instance for this module (optional, but good practice)
# If you don't get a specific logger, logging.info(), logging.error() etc. use the root logger.
logger = logging.getLogger(__name__) # Using __name__ gives you a logger named 'data_fetching'

# --- Constants ---
ORTHODB_FASTA_URL = "https://data.orthodb.org/v12/fasta"
ORTHODB_GENE_SEARCH_URL = "https://data.orthodb.org/v12/genesearch"
CACHE_DIR = "backend/cached"

# --- Utility Function (Consider moving to a separate file_utils.py) ---
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
            logger.info(f"Appending string data to: {path}")
            with open(path, 'a') as f:
                f.write(str(data) + "\n")
    except IOError as e:
        logger.error(f"Error writing to cache file {path}: {e}", exc_info=True)


# --- Data Fetching Functions ---
def pull_OG_fasta(ortholog_group: str) -> Optional[str]:
    """
    Retrieves FASTA formatted sequences for a given ortholog group from OrthoDB.

    This function first checks if the data is available in a local cache.
    If not, it fetches the data from the OrthoDB API and then caches it.

    Args:
        ortholog_group: The identifier for the ortholog group (e.g., "430340at2759").

    Returns:
        A string containing the FASTA sequences if successful, or None if an error occurs
        or the ortholog group is not found.
    """
    cache_filename = f"{ortholog_group}.fasta"
    cache_filepath = os.path.join(CACHE_DIR, cache_filename)

    if os.path.exists(cache_filepath):
        logger.info(f"Data for ortholog group '{ortholog_group}' found in cache: {cache_filepath}")
        try:
            with open(cache_filepath, "r") as f:
                return f.read()
        except IOError as e:
            logger.error(f"Error reading from cache file {cache_filepath}: {e}", exc_info=True)

    logger.info(f"Pulling FASTA data for ortholog group '{ortholog_group}' from OrthoDB...")
    params = {
        'id': ortholog_group,
        'seqtype': 'cds'
    }
    try:
        response = requests.get(ORTHODB_FASTA_URL, params=params)
        response.raise_for_status()
        
        response_text = response.text
        if not response_text or response_text.isspace():
            logger.warning(f"Received empty response for ortholog group '{ortholog_group}'.")
            return None
            
        cache_data(response_text, cache_filename)
        logger.info(f"Ortholog group data for '{ortholog_group}' pulled and cached.")
        return response_text
        
    except requests.exceptions.HTTPError as http_err:
        logger.error(f"HTTP error occurred while fetching FASTA for '{ortholog_group}': {http_err} - Status code: {response.status_code}", exc_info=True)
        if response.status_code == 404:
            logger.warning(f"Ortholog group '{ortholog_group}' not found on OrthoDB.")
        return None
    except requests.exceptions.RequestException as req_err:
        logger.error(f"Request error occurred while fetching FASTA for '{ortholog_group}': {req_err}", exc_info=True)
        return None
    except Exception as e:
        logger.exception(f"An unexpected error occurred while fetching FASTA for '{ortholog_group}'.") # .exception includes exc_info=True
        return None


def pull_model_organism_orthologs(gene_query: str) -> Optional[Dict[str, Any]]:
    """
    Retrieves ortholog information for a given gene query from OrthoDB,
    specifically focusing on model organisms.

    Args:
        gene_query: The gene identifier or name to search for (e.g., "WBGene00004963").

    Returns:
        A dictionary containing the ortholog data if successful, or None if an error occurs.
    """
    cache_filename = f"{gene_query}_search.json"
    cache_filepath = os.path.join(CACHE_DIR, cache_filename)

    if os.path.exists(cache_filepath):
        logger.info(f"Data for gene query '{gene_query}' found in cache: {cache_filepath}")
        try:
            with open(cache_filepath, 'r') as file:
                response_json = json.load(file)
            if "orthologs_in_model_organisms" in response_json:
                orthologs = response_json["orthologs_in_model_organisms"]
                ortholog_count = sum(len(organism.get("genes", [])) for organism in orthologs)
                logger.info(f"{ortholog_count} orthologs found in model organisms (from cache).")
                return response_json
            else:
                logger.warning(f"Cached data for '{gene_query}' is missing 'orthologs_in_model_organisms' key.")
        except (IOError, json.JSONDecodeError) as e:
            logger.error(f"Error reading or parsing cache file {cache_filepath}: {e}", exc_info=True)

    logger.info(f"Pulling model organism orthologs for gene query '{gene_query}' from OrthoDB...")
    gene_search_url = f"{ORTHODB_GENE_SEARCH_URL}?query={gene_query}"
    
    try:
        response = requests.get(gene_search_url)
        response.raise_for_status()
        
        response_json = response.json()
        cache_data(response_json, cache_filename)

        if "orthologs_in_model_organisms" in response_json:
            orthologs = response_json["orthologs_in_model_organisms"]
            ortholog_count = sum(len(organism.get("genes", [])) for organism in orthologs)
            logger.info(f"{ortholog_count} orthologs found in model organisms and cached.")
            return response_json
        else:
            logger.warning(f"'orthologs_in_model_organisms' key not found in API response for '{gene_query}'.")
            return None

    except requests.exceptions.HTTPError as http_err:
        logger.error(f"HTTP error for '{gene_query}': {http_err} - Status: {response.status_code}", exc_info=True)
        if response.status_code == 404:
             logger.warning(f"Gene query '{gene_query}' not found on OrthoDB.")
        return None
    except requests.exceptions.RequestException as req_err:
        logger.error(f"Request error for '{gene_query}': {req_err}", exc_info=True)
        return None
    except json.JSONDecodeError as json_err:
        logger.error(f"Error decoding JSON for '{gene_query}': {json_err}", exc_info=True)
        return None
    except Exception as e:
        logger.exception(f"An unexpected error occurred for '{gene_query}'.")
        return None

if __name__ == '__main__':
    logger.info("--- Starting Data Fetching Module Tests ---")
    
    # Test pull_OG_fasta
    logger.info("\n--- Testing pull_OG_fasta ---")
    og_id = "430340at2759"
    fasta_data = pull_OG_fasta(og_id)
    if fasta_data:
        logger.info(f"Successfully fetched FASTA for {og_id}. First 100 chars: {fasta_data[:100]}...")
    else:
        logger.error(f"Failed to fetch FASTA for {og_id}.")

    # Test pull_model_organism_orthologs
    logger.info("\n--- Testing pull_model_organism_orthologs ---")
    gene_query = "WBGene00004963"
    ortholog_data = pull_model_organism_orthologs(gene_query)
    if ortholog_data and "orthologs_in_model_organisms" in ortholog_data:
        logger.info(f"Successfully fetched orthologs for {gene_query}.")
    else:
        logger.error(f"Failed to fetch orthologs for {gene_query} or data is malformed.")

    logger.info("\n--- Testing with a non-existent gene ---")
    non_existent_gene = "NonExistentGene123XYZ"
    ortholog_data_non_existent = pull_model_organism_orthologs(non_existent_gene)
    if not ortholog_data_non_existent:
        logger.info(f"Correctly handled non-existent gene: {non_existent_gene}")

    logger.info("\n--- Testing with a non-existent ortholog group ---")
    non_existent_og = "0at0"
    fasta_data_non_existent = pull_OG_fasta(non_existent_og)
    if not fasta_data_non_existent:
        logger.info(f"Correctly handled non-existent ortholog group: {non_existent_og}")

    logger.info("--- Data Fetching Module Tests Finished ---")

