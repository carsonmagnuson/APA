import logging
import os
import json
import urllib.parse
import re # Ensure re is imported for the download endpoint

from fastapi import FastAPI, HTTPException, Body, status, Path as FastAPIPath
from fastapi.responses import FileResponse, JSONResponse
# Import CORSMiddleware
from fastapi.middleware.cors import CORSMiddleware # <<<< ADD THIS IMPORT
from pydantic import BaseModel
from typing import Optional, Dict, Any
import uvicorn

try:
    from main_pipeline import (
        run_phylogenetic_pipeline,
        CACHE_DIR,
        DEFAULT_PAML_BASE_OUTPUT_DIR,
        DEFAULT_OUTPUT_PHYLIP_FILENAME,
        DEFAULT_OUTPUT_NEWICK_FILENAME
    )
except ImportError as e:
    print(f"Error importing modules: {e}.")
    raise

# --- Configure Logging ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(name)s - %(module)s.%(funcName)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# --- Pydantic Model ---
class GeneAnalysisRequest(BaseModel):
    ncbi_id: str
    k_sequences: Optional[int] = 10
    ortho_level_tax_id: Optional[str] = "2759"
    rooted_tree: Optional[bool] = False

# --- FastAPI Application ---
app = FastAPI(
    title="Phylogenetic Analysis API",
    description="An API to run phylogenetic analysis and CODEML.",
    version="0.1.0"
)

# --- CORS Configuration ---
# Define the origins that are allowed to make requests to your API.
# For development, this will be your Vite frontend's address.
origins = [
    "http://localhost:5173", # Vite's default dev server
    "http://127.0.0.1:5173",
    # Add other origins if needed, e.g., your deployed frontend URL
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,  # List of origins that are allowed to make requests
    allow_credentials=True, # Allows cookies to be included in requests
    allow_methods=["*"],    # Allows all methods (GET, POST, OPTIONS, etc.)
    allow_headers=["*"],    # Allows all headers
)
# --- End CORS Configuration ---


@app.on_event("startup")
async def startup_event():
    if not os.path.exists(CACHE_DIR):
        try:
            os.makedirs(CACHE_DIR)
            logger.info(f"Created cache directory: {CACHE_DIR}")
        except OSError as e:
            logger.error(f"Error creating cache directory {CACHE_DIR}: {e}")

    if not os.path.exists(DEFAULT_PAML_BASE_OUTPUT_DIR):
        try:
            os.makedirs(DEFAULT_PAML_BASE_OUTPUT_DIR)
            logger.info(f"Created default PAML base output directory: {DEFAULT_PAML_BASE_OUTPUT_DIR}")
        except OSError as e:
            logger.error(f"Error creating PAML base output directory {DEFAULT_PAML_BASE_OUTPUT_DIR}: {e}")

@app.post("/analyze/", response_model=Optional[Dict[str, Any]])
async def analyze_gene(request_data: GeneAnalysisRequest = Body(...)):
    # ... (your existing /analyze/ endpoint logic)
    logger.info(f"Received analysis request for NCBI ID: {request_data.ncbi_id}, Level: {request_data.ortho_level_tax_id}")
    try:
        results = run_phylogenetic_pipeline(
            ncbi_id=request_data.ncbi_id,
            k_sequences=request_data.k_sequences if request_data.k_sequences is not None else 10,
            ortho_level_tax_id=request_data.ortho_level_tax_id if request_data.ortho_level_tax_id else "2759",
            paml_base_output_dir=DEFAULT_PAML_BASE_OUTPUT_DIR,
            rooted_tree=request_data.rooted_tree if request_data.rooted_tree is not None else False
        )
        # ... (rest of your /analyze/ logic for handling results and errors)
        if results:
            codeml_summary = results.get("codeml_analysis_summary", {})
            if codeml_summary.get("status") and "ERROR" in codeml_summary.get("status", "").upper():
                logger.error(f"Pipeline for {request_data.ncbi_id} encountered a CODEML error: {codeml_summary.get('message')}")
                return JSONResponse(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, content=results)
            elif results.get("status") and "ERROR" in results.get("status", "").upper():
                logger.error(f"Pipeline for {request_data.ncbi_id} encountered an error: {results.get('message')}")
                return JSONResponse(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, content=results)
            else:
                logger.info(f"Analysis for NCBI ID: {request_data.ncbi_id} completed.")
            return results
        else:
            logger.error(f"Analysis for NCBI ID: {request_data.ncbi_id} returned None from pipeline.")
            raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Analysis pipeline failed to produce results.")
    except Exception as e:
        logger.error(f"An unexpected error occurred in API for {request_data.ncbi_id}: {e}", exc_info=True)
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"An internal server error occurred: {str(e)}")


# Allowed filenames for download to prevent arbitrary file access
ALLOWED_FILENAMES = {
    "results.out": "text/plain",
    DEFAULT_OUTPUT_PHYLIP_FILENAME: "text/plain",
    DEFAULT_OUTPUT_NEWICK_FILENAME: "text/plain",
    "rst": "text/plain", 
    "rst1": "text/plain",
    "rub": "text/plain",
    "lnf": "text/plain",
    "2NG.t": "text/plain",
    "4fold.nuc": "text/plain"
}

@app.get("/download_paml_output/{run_identifier}/{filename}")
async def download_paml_file(
    run_identifier: str = FastAPIPath(..., description="The unique identifier for the PAML run"),
    filename: str = FastAPIPath(..., description=f"The name of the file to download. Allowed: {', '.join(ALLOWED_FILENAMES.keys())}")
):
    logger.info(f"Download request for run: '{run_identifier}', file: '{filename}'")
    if filename not in ALLOWED_FILENAMES:
        logger.warning(f"Attempt to download disallowed file: {filename}")
        raise HTTPException(status_code=400, detail="Invalid or disallowed filename.")
    if not re.match(r"^[a-zA-Z0-9_.-]+$", run_identifier) or ".." in run_identifier or "/" in run_identifier:
        logger.warning(f"Invalid run_identifier format: {run_identifier}")
        raise HTTPException(status_code=400, detail="Invalid run identifier format.")
    
    base_dir = os.path.abspath(DEFAULT_PAML_BASE_OUTPUT_DIR)
    file_path = os.path.abspath(os.path.join(base_dir, run_identifier, filename))

    if not file_path.startswith(base_dir):
        logger.error(f"Path traversal attempt detected. Requested path: {file_path}, Base path: {base_dir}")
        raise HTTPException(status_code=403, detail="Access forbidden.")
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        logger.error(f"File not found for download: {file_path}")
        raise HTTPException(status_code=404, detail=f"File '{filename}' not found for run '{run_identifier}'.")

    media_type = ALLOWED_FILENAMES[filename]
    return FileResponse(path=file_path, filename=filename, media_type=media_type)

@app.get("/")
async def read_root():
    return {"message": "Phylogenetic Analysis API is running."}

if __name__ == "__main__":
    logger.info("Starting FastAPI server for Phylogenetic Analysis...")
    uvicorn.run("api_server:app", host="0.0.0.0", port=8000, reload=True)

