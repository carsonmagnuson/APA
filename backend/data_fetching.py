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
ORTHODB_SEARCH_URL = "https://data.orthodb.org/v12/search" # New constant for the search endpoint
CACHE_DIR = "backend/cached"

# --- Utility Function (Consider moving to a separate file_utils.py if not already there) ---
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
            logger.info(f"Writing string data to: {path}") # Changed from 'appending' for clarity
            with open(path, 'w') as f: # Use 'w' to overwrite, typical for caching new fetch
                f.write(str(data))
                if not isinstance(data, str) or not data.endswith("\n"):
                     f.write("\n") # Add newline if not already present or not a string that ends with it
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
            # Proceed to fetch from API if cache read fails

    logger.info(f"Pulling FASTA data for ortholog group '{ortholog_group}' from OrthoDB...")
    params = {
        'id': ortholog_group,
        'seqtype': 'cds'
    }
    response = None  # Initialize response to None
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
        logger.error(f"HTTP error occurred while fetching FASTA for '{ortholog_group}': {http_err} - Status code: {response.status_code if response else 'N/A'}", exc_info=True)
        if response and response.status_code == 404:
            logger.warning(f"Ortholog group '{ortholog_group}' not found on OrthoDB.")
        return None
    except requests.exceptions.RequestException as req_err:
        logger.error(f"Request error occurred while fetching FASTA for '{ortholog_group}': {req_err}", exc_info=True)
        return None
    except Exception as e:
        logger.exception(f"An unexpected error occurred while fetching FASTA for '{ortholog_group}'.")
        return None


def pull_model_organism_orthologs(gene_query: str) -> Optional[Dict[str, Any]]:
    """
    Retrieves ortholog information for a given gene query from OrthoDB /genesearch,
    specifically focusing on model organisms.

    Args:
        gene_query: The gene identifier or name to search for (e.g., "WBGene00004963", or an NCBI ID like "173042").

    Returns:
        A dictionary containing the gene search data if successful, or None if an error occurs.
        This dictionary includes the primary gene's details and its orthologs in model organisms.
    """
    cache_filename = f"{gene_query}_search.json" # Cache based on the query itself
    cache_filepath = os.path.join(CACHE_DIR, cache_filename)

    if os.path.exists(cache_filepath):
        logger.info(f"Gene search data for query '{gene_query}' found in cache: {cache_filepath}")
        try:
            with open(cache_filepath, 'r') as file:
                response_json = json.load(file)
            # Basic validation of cached structure
            if "gene" in response_json and "orthologs_in_model_organisms" in response_json:
                orthologs = response_json["orthologs_in_model_organisms"]
                ortholog_count = sum(len(organism.get("genes", [])) for organism in orthologs)
                logger.info(f"{ortholog_count} orthologs in model organisms found (from cache for query '{gene_query}').")
                return response_json
            else:
                logger.warning(f"Cached data for '{gene_query}' is missing expected keys ('gene' or 'orthologs_in_model_organisms'). Refetching.")
        except (IOError, json.JSONDecodeError) as e:
            logger.error(f"Error reading or parsing cache file {cache_filepath}: {e}", exc_info=True)
            # Proceed to fetch from API if cache is corrupted

    logger.info(f"Pulling gene search data for query '{gene_query}' from OrthoDB /genesearch...")
    # The API documentation suggests 'query' for general search, 'ncbi' or 'gid' for NCBI gene IDs.
    # Let's assume gene_query could be an NCBI ID or other gene name.
    # If it's specifically an NCBI gene ID, using the 'ncbi' parameter might be more precise.
    # For simplicity and consistency with original behavior, we'll use 'query'.
    # If this needs to be more specific (e.g., distinguish between NCBI ID and gene name),
    # the calling function (main_pipeline) should handle that or pass a type.
    params = {'query': gene_query}
    response = None # Initialize response
    try:
        response = requests.get(ORTHODB_GENE_SEARCH_URL, params=params)
        response.raise_for_status()
        
        response_json = response.json()
        
        # Validate structure before caching
        if "gene" in response_json and "orthologs_in_model_organisms" in response_json:
            cache_data(response_json, cache_filename) # Cache the valid response
            orthologs = response_json["orthologs_in_model_organisms"]
            ortholog_count = sum(len(organism.get("genes", [])) for organism in orthologs)
            logger.info(f"{ortholog_count} orthologs in model organisms found for '{gene_query}' and cached.")
            return response_json
        else:
            logger.error(f"API response for '{gene_query}' is missing expected keys ('gene' or 'orthologs_in_model_organisms'). Response: {response_json}")
            # Cache the problematic response for debugging
            cache_data(response_json, cache_filename)
            return None

    except requests.exceptions.HTTPError as http_err:
        logger.error(f"HTTP error for gene query '{gene_query}': {http_err} - Status: {response.status_code if response else 'N/A'}", exc_info=True)
        if response and response.status_code == 404:
             logger.warning(f"Gene query '{gene_query}' not found on OrthoDB /genesearch.")
        return None
    except requests.exceptions.RequestException as req_err:
        logger.error(f"Request error for gene query '{gene_query}': {req_err}", exc_info=True)
        return None
    except json.JSONDecodeError as json_err:
        logger.error(f"Error decoding JSON for gene query '{gene_query}': {json_err}. Response text: {response.text if response else 'No response'}", exc_info=True)
        return None
    except Exception as e:
        logger.exception(f"An unexpected error occurred for gene query '{gene_query}'.")
        return None

def pull_ortholog_group_ids_by_gene_and_level(
    gene_identifier: str,
    level_tax_id: str,
    identifier_type: str = "ncbi"
) -> Optional[List[str]]:
    """
    Retrieves OrthoDB Ortholog Group (OG) IDs for a given gene identifier
    at a specific taxonomic level using the /v12/search endpoint.

    Args:
        gene_identifier: The gene identifier (e.g., NCBI ID, OrthoDB gene ID).
        level_tax_id: The NCBI taxonomy ID for the level of orthology.
        identifier_type: Type of gene_identifier. Can be "ncbi" for NCBI gene ID
                         or "query" for OrthoDB gene ID or other text. Defaults to "ncbi".

    Returns:
        A list of ortholog group ID strings if successful, or None if an error occurs
        or no OGs are found.
    """
    # Sanitize gene_identifier for filename if it contains special characters
    safe_gene_identifier = "".join(c if c.isalnum() else "_" for c in gene_identifier)
    cache_filename = f"search_{safe_gene_identifier}_level_{level_tax_id}.json"
    cache_filepath = os.path.join(CACHE_DIR, cache_filename)

    if os.path.exists(cache_filepath):
        logger.info(f"Data for gene search '{gene_identifier}' at level '{level_tax_id}' found in cache: {cache_filepath}")
        try:
            with open(cache_filepath, 'r') as file:
                response_json = json.load(file)
            if "data" in response_json and isinstance(response_json["data"], list):
                og_ids = response_json["data"]
                logger.info(f"Found {len(og_ids)} OG IDs in cache for '{gene_identifier}' at level '{level_tax_id}'.")
                return og_ids
            else:
                logger.warning(f"Cached data for '{gene_identifier}' at level '{level_tax_id}' is missing 'data' key or it's not a list. Refetching.")
        except (IOError, json.JSONDecodeError) as e:
            logger.error(f"Error reading or parsing cache file {cache_filepath}: {e}", exc_info=True)
            # Proceed to fetch from API if cache is corrupted

    logger.info(f"Pulling OG IDs for gene '{gene_identifier}' (type: {identifier_type}) at level '{level_tax_id}' from OrthoDB /search...")
    params = {
        identifier_type: gene_identifier,
        'level': level_tax_id,
        # 'seqtype': 'cds' # seqtype is not a valid parameter for /search endpoint
    }
    
    response = None # Initialize response
    try:
        response = requests.get(ORTHODB_SEARCH_URL, params=params)
        response.raise_for_status() # Raises an HTTPError for bad responses (4XX or 5XX)
        
        response_json = response.json()
        
        if "data" in response_json and isinstance(response_json["data"], list):
            og_ids = response_json["data"]
            # Cache the response regardless of whether OGs were found, to avoid re-fetching known empty/problematic results
            cache_data(response_json, cache_filename)
            if og_ids:
                logger.info(f"{len(og_ids)} OG IDs found for '{gene_identifier}' at level '{level_tax_id}' and cached.")
                return og_ids
            else:
                logger.warning(f"No OG IDs found in 'data' array for '{gene_identifier}' at level '{level_tax_id}'. Response: {response_json}")
                return [] # Return empty list if 'data' is empty, indicating no OGs found
        else:
            logger.error(f"'data' key missing or not a list in API response for '{gene_identifier}' at level '{level_tax_id}'. Response: {response_json}")
            cache_data(response_json, cache_filename) # Cache the problematic response
            return None

    except requests.exceptions.HTTPError as http_err:
        logger.error(f"HTTP error for search '{gene_identifier}', level '{level_tax_id}': {http_err} - Status: {response.status_code if response else 'N/A'}", exc_info=True)
        if response and response.status_code == 404:
             logger.warning(f"Search for '{gene_identifier}', level '{level_tax_id}' not found on OrthoDB /search.")
        # Cache the error response if possible, or an error marker
        if response is not None:
            try:
                cache_data(response.json(), cache_filename)
            except json.JSONDecodeError:
                cache_data({"error": "HTTPError", "status_code": response.status_code, "text": response.text}, cache_filename)
        return None
    except requests.exceptions.RequestException as req_err:
        logger.error(f"Request error for search '{gene_identifier}', level '{level_tax_id}': {req_err}", exc_info=True)
        return None
    except json.JSONDecodeError as json_err:
        logger.error(f"Error decoding JSON for search '{gene_identifier}', level '{level_tax_id}': {json_err}. Response text: {response.text if response else 'No response'}", exc_info=True)
        if response is not None: # Cache the non-JSON response
             cache_data({"error": "JSONDecodeError", "text": response.text}, cache_filename)
        return None
    except Exception as e: # Catch-all for other unexpected errors
        logger.exception(f"An unexpected error occurred during pull_ortholog_group_ids_by_gene_and_level for '{gene_identifier}', level '{level_tax_id}'.")
        return None


if __name__ == '__main__':
    # Ensure cache directory exists for standalone testing
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)
        logger.info(f"Created cache directory for testing: {CACHE_DIR}")

    logger.info("--- Starting Data Fetching Module Tests ---")
    
    # Test pull_OG_fasta
    logger.info("\n--- Testing pull_OG_fasta ---")
    og_id_fasta_test = "430340at2759" # Eukaryota level OG
    fasta_data = pull_OG_fasta(og_id_fasta_test)
    if fasta_data:
        logger.info(f"Successfully fetched FASTA for {og_id_fasta_test}. First 100 chars: {fasta_data[:100]}...")
    else:
        logger.error(f"Failed to fetch FASTA for {og_id_fasta_test}.")

    # Test pull_model_organism_orthologs (uses /genesearch)
    logger.info("\n--- Testing pull_model_organism_orthologs ---")
    gene_query_search = "173042" # NCBI Gene ID for spe-39 in C. elegans
    ortholog_data = pull_model_organism_orthologs(gene_query_search)
    if ortholog_data and "gene" in ortholog_data and "orthologs_in_model_organisms" in ortholog_data:
        logger.info(f"Successfully fetched /genesearch data for {gene_query_search}.")
        # Example of extracting target_gene_orthodb_id
        target_gene_orthodb_id_test = ortholog_data.get("gene", {}).get("gene_id", {}).get("param")
        logger.info(f"Target OrthoDB ID from /genesearch: {target_gene_orthodb_id_test}")
    else:
        logger.error(f"Failed to fetch /genesearch data for {gene_query_search} or data is malformed.")

    # Test new function: pull_ortholog_group_ids_by_gene_and_level (uses /search)
    logger.info("\n--- Testing pull_ortholog_group_ids_by_gene_and_level ---")
    # Using an NCBI gene ID (e.g., for human hemoglobin subunit alpha 1) and a relevant level (e.g., Mammalia)
    test_ncbi_gene_id = "3039" # NCBI Gene ID for HBA1 (Human)
    test_level_tax_id = "40674"  # NCBI Taxonomy ID for Mammalia
    
    og_ids = pull_ortholog_group_ids_by_gene_and_level(
        gene_identifier=test_ncbi_gene_id,
        level_tax_id=test_level_tax_id,
        identifier_type="ncbi"
    )
    if og_ids is not None:
        if og_ids:
            logger.info(f"Successfully fetched {len(og_ids)} OG IDs for NCBI gene {test_ncbi_gene_id} at level {test_level_tax_id}. First few: {og_ids[:5]}")
        else:
            logger.info(f"No OG IDs found for NCBI gene {test_ncbi_gene_id} at level {test_level_tax_id}, but the API call was successful.")
    else:
        logger.error(f"Failed to fetch OG IDs for NCBI gene {test_ncbi_gene_id} at level {test_level_tax_id}.")

    logger.info("\n--- Testing with a non-existent gene for /search ---")
    non_existent_gene_search = "NonExistentGene123XYZ"
    og_ids_non_existent = pull_ortholog_group_ids_by_gene_and_level(non_existent_gene_search, "2759", "query")
    if og_ids_non_existent is not None and not og_ids_non_existent: # Expecting empty list for "not found"
        logger.info(f"Correctly handled non-existent gene for /search: {non_existent_gene_search} (returned empty list).")
    elif og_ids_non_existent is None:
        logger.error(f"Error occurred for non-existent gene search: {non_existent_gene_search}.")
    else:
        logger.warning(f"Unexpectedly found OGs for non-existent gene: {og_ids_non_existent}")


    logger.info("\n--- Testing with a non-existent ortholog group for /fasta ---")
    non_existent_og_fasta = "0at0" # Bogus OG ID
    fasta_data_non_existent = pull_OG_fasta(non_existent_og_fasta)
    if not fasta_data_non_existent:
        logger.info(f"Correctly handled non-existent ortholog group for /fasta: {non_existent_og_fasta}")

    logger.info("--- Data Fetching Module Tests Finished ---")

