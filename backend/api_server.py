import logging
import os
from typing import Dict, Any, Optional

from fastapi import FastAPI, HTTPException, Body, status
from fastapi.responses import JSONResponse
from pydantic import BaseModel
import uvicorn
# Assuming your custom modules are in the same directory or accessible via PYTHONPATH
# Ensure these modules and their dependencies (like BioPython, requests, PAML) are installed
try:
    from main_pipeline import (
        run_phylogenetic_pipeline,
        CACHE_DIR, # Import for ensuring directory existence
        DEFAULT_PAML_OUTPUT_DIR,
        DEFAULT_OUTPUT_PHYLIP_FILE,
        DEFAULT_OUTPUT_NEWICK_FILE
    )
except ImportError as e:
    print(f"Error importing modules: {e}. Make sure main_pipeline.py and its dependencies are accessible.")
    print("Ensure you are in the 'backend' directory or have set up PYTHONPATH correctly.")
    raise

# --- Configure Logging for the FastAPI app ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(name)s - %(module)s.%(funcName)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# --- Pydantic Model for Request Body ---
class GeneAnalysisRequest(BaseModel):
    """
    Defines the expected request body for the gene analysis endpoint.
    """
    ncbi_id: str  # NCBI ID, e.g., a GenBank accession or NCBI Gene ID that OrthoDB can search.
    k_sequences: Optional[int] = 10
    cleanup_paml_dir: Optional[bool] = True
    rooted_tree: Optional[bool] = False

# --- FastAPI Application ---
app = FastAPI(
    title="Phylogenetic Analysis API",
    description="An API to run phylogenetic analysis and CODEML.",
    version="0.1.0"
)

@app.on_event("startup")
async def startup_event():
    """
    Actions to perform on application startup.
    Ensures necessary directories exist.
    """
    if not os.path.exists(CACHE_DIR):
        try:
            os.makedirs(CACHE_DIR)
            logger.info(f"Created cache directory: {CACHE_DIR}")
        except OSError as e:
            logger.error(f"Error creating cache directory {CACHE_DIR}: {e}")

    if not os.path.exists(DEFAULT_PAML_OUTPUT_DIR):
        try:
            os.makedirs(DEFAULT_PAML_OUTPUT_DIR)
            logger.info(f"Created default PAML output directory: {DEFAULT_PAML_OUTPUT_DIR}")
        except OSError as e:
            logger.error(f"Error creating PAML output directory {DEFAULT_PAML_OUTPUT_DIR}: {e}")


@app.post("/analyze/", response_model=Optional[Dict[str, Any]])
async def analyze_gene(request_data: GeneAnalysisRequest = Body(...)):
    """
    Endpoint to perform phylogenetic analysis on a given gene, identified by NCBI ID.
    """
    logger.info(f"Received analysis request for NCBI ID: {request_data.ncbi_id}")

    phylip_file = os.path.join(DEFAULT_PAML_OUTPUT_DIR, DEFAULT_OUTPUT_PHYLIP_FILE)
    newick_file = os.path.join(DEFAULT_PAML_OUTPUT_DIR, DEFAULT_OUTPUT_NEWICK_FILE)

    try:
        results = run_phylogenetic_pipeline(
            ncbi_id=request_data.ncbi_id,
            k_sequences=request_data.k_sequences if request_data.k_sequences is not None else 10,
            output_phylip_file=phylip_file,
            output_newick_file=newick_file,
            paml_output_dir=DEFAULT_PAML_OUTPUT_DIR,
            cleanup_paml_dir=request_data.cleanup_paml_dir if request_data.cleanup_paml_dir is not None else True,
            rooted_tree=request_data.rooted_tree if request_data.rooted_tree is not None else False
        )

        if results:
            status_code = status.HTTP_200_OK
            if results.get("status") and "FAILED" in results.get("status", "").upper():
                logger.warning(f"Pipeline for {request_data.ncbi_id} completed with issues: {results.get('message')}")
            elif results.get("status") and "ERROR" in results.get("status", "").upper():
                logger.error(f"Pipeline for {request_data.ncbi_id} encountered an error: {results.get('message')}")
                return JSONResponse(
                    status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                    content=results
                )
            else:
                 logger.info(f"Analysis successful for NCBI ID: {request_data.ncbi_id}")

            return results
        else:
            logger.error(f"Analysis for NCBI ID: {request_data.ncbi_id} returned None from pipeline, indicating a setup or critical error.")
            raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Analysis pipeline failed to initialize or encountered a critical error before CODEML execution.")

    except FileNotFoundError as fnf_error:
        logger.error(f"FileNotFoundError during analysis setup: {fnf_error}", exc_info=True)
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"A required file/directory was not found: {fnf_error}")
    # Removed requests.exceptions.RequestException as it's better handled within the pipeline
    except Exception as e:
        logger.error(f"An unexpected error occurred in API for {request_data.ncbi_id}: {e}", exc_info=True)
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"An internal server error occurred: {str(e)}")

@app.get("/")
async def read_root():
    """
    Root endpoint to check if the API is running.
    """
    return {"message": "Phylogenetic Analysis API is running. Use the /analyze/ endpoint to submit jobs."}

if __name__ == "__main__":
    logger.info("Starting FastAPI server for Phylogenetic Analysis...")
    uvicorn.run("api_server:app", host="0.0.0.0", port=8000, reload=True)
