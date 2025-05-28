// src/App.tsx
import { useState } from 'react';
import axios from 'axios';
import './App.css';

// Define constants for filenames that match your backend's ALLOWED_FILENAMES
//const FILENAME_RESULTS_OUT = 'results.out';
const FILENAME_PHYLIP = 'output.phylip'; // From DEFAULT_OUTPUT_PHYLIP_FILENAME in backend
const FILENAME_NEWICK = 'tree.nwk';     // From DEFAULT_OUTPUT_NEWICK_FILENAME in backend

interface CodemlAnalysisSummary {
  status?: string;
  message?: string;
  results_file_path_on_server?: string;
  results_filename?: string; // This might be present if run_codeml_positive_selection adds it
  error_details?: string;
  M0?: { lnL?: number; parameters?: any };
  M1a?: { lnL?: number; parameters?: any };
  M2a?: { lnL?: number; parameters?: any };
}

interface AnalysisResult {
  input_ncbi_id?: string;
  target_gene_orthodb_id?: string;
  ortholog_group_id_used?: string;
  ortho_level_tax_id_used?: string;
  k_sequences_requested?: number;
  sequences_in_alignment?: number;
  paml_run_identifier?: string;
  paml_results_filename?: string; // This is the top-level filename for 'results.out' from main_pipeline.py
  codeml_analysis_summary?: CodemlAnalysisSummary;
}

function App() {
  const [ncbiId, setNcbiId] = useState<string>('');
  const [orthoLevel, setOrthoLevel] = useState<string>('2759');
  const [results, setResults] = useState<AnalysisResult | null>(null);
  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);

  const BACKEND_URL = 'http://localhost:8000';

  const handleAnalyze = async () => {
    if (!ncbiId.trim()) {
      setError('Please enter an NCBI Gene ID.');
      return;
    }
    if (!orthoLevel.trim()) {
        setError('Please enter an Orthology Level Tax ID.');
        return;
    }
    setIsLoading(true);
    setError(null);
    setResults(null);
    console.log(`Sending analysis request for NCBI ID: ${ncbiId}, Level: ${orthoLevel}`);

    try {
      const response = await axios.post(`${BACKEND_URL}/analyze/`, {
        ncbi_id: ncbiId,
        ortho_level_tax_id: orthoLevel,
      });
      console.log("Backend response received:", response.data);
      setResults(response.data);
    } catch (err: any) {
      let errorMsg = 'An unexpected error occurred.';
      if (axios.isAxiosError(err) && err.response) {
        errorMsg = err.response.data.detail || 'An error occurred during analysis.';
        console.error('Analysis error response:', err.response.data);
      } else {
        console.error('Unexpected error during analysis:', err);
      }
      setError(errorMsg);
      setResults(null);
    } finally {
      setIsLoading(false);
    }
  };

  const handleDownloadFile = (filenameToDownload: string | undefined) => {
    console.log("handleDownloadFile called with:", filenameToDownload);
    console.log("Current results state:", results);

    const runIdentifier = results?.paml_run_identifier;
    
    if (runIdentifier && filenameToDownload) {
      const downloadUrl = `${BACKEND_URL}/download_paml_output/${runIdentifier}/${filenameToDownload}`;
      console.log("Attempting to download from URL:", downloadUrl);
      try {
        const newWindow = window.open(downloadUrl, '_blank');
        if (!newWindow || newWindow.closed || typeof newWindow.closed === 'undefined') {
            console.error('window.open failed. Pop-up blocker? Malformed URL?');
            setError(`Failed to initiate download for ${filenameToDownload}. Check browser pop-up blocker or console for errors.`);
        }
      } catch (e) {
        console.error("Error during window.open for download:", e);
        setError(`Error trying to open download link for ${filenameToDownload}.`);
      }
    } else {
      console.error('Download conditions not met:', {
        hasResults: !!results,
        runIdentifier: runIdentifier,
        filename: filenameToDownload,
      });
      setError(`Cannot download ${filenameToDownload || 'file'}: missing run identifier or filename in results.`);
    }
  };

  return (
    <div className="App">
      <header className="App-header">
        <h1>Phylogenetic Analysis</h1>
      </header>
      <main>
        <div>
          <label htmlFor="ncbiId">NCBI Gene ID:</label>
          <input
            type="text"
            id="ncbiId"
            value={ncbiId}
            onChange={(e) => setNcbiId(e.target.value)}
            placeholder="e.g., 173042"
          />
        </div>
        <div>
          <label htmlFor="orthoLevel">Orthology Level (NCBI TaxID):</label>
          <input
            type="text"
            id="orthoLevel"
            value={orthoLevel}
            onChange={(e) => setOrthoLevel(e.target.value)}
            placeholder="e.g., 2759 (Eukaryota)"
          />
        </div>
        <button onClick={handleAnalyze} disabled={isLoading}>
          {isLoading ? 'Analyzing...' : 'Run Analysis'}
        </button>

        {error && <p className="error-message">Error: {error}</p>}

        {isLoading && <p>Loading analysis results...</p>}

        {results && !isLoading && (
          <div className="results">
            <h2>Analysis Results for NCBI ID: {results.input_ncbi_id}</h2>
            <p>Target OrthoDB ID: {results.target_gene_orthodb_id || 'N/A'}</p>
            <p>Ortholog Group Used: {results.ortholog_group_id_used || 'N/A'} (Level: {results.ortho_level_tax_id_used || 'N/A'})</p>
            <p>Sequences Requested: {results.k_sequences_requested !== undefined ? results.k_sequences_requested : 'N/A'}</p>
            <p>Sequences in Final Alignment: {results.sequences_in_alignment !== undefined ? results.sequences_in_alignment : 'N/A'}</p>
            <p>PAML Run Identifier: {results.paml_run_identifier || 'N/A'}</p>
            
            {results.codeml_analysis_summary && (
              <div>
                <h3>CODEML Analysis Status:</h3>
                <p>Status: {results.codeml_analysis_summary.status || 'N/A'}</p>
                <p>Message: {results.codeml_analysis_summary.message || 'N/A'}</p>
                {results.codeml_analysis_summary.error_details && (
                  <p>Error Details: {results.codeml_analysis_summary.error_details}</p>
                )}
                {results.codeml_analysis_summary.M0?.lnL !== undefined && (
                    <p>M0 lnL: {results.codeml_analysis_summary.M0.lnL}</p>
                )}
                {results.codeml_analysis_summary.M1a?.lnL !== undefined && (
                    <p>M1a lnL: {results.codeml_analysis_summary.M1a.lnL}</p>
                )}
                 {results.codeml_analysis_summary.M2a?.lnL !== undefined && (
                    <p>M2a lnL: {results.codeml_analysis_summary.M2a.lnL}</p>
                )}

                {/* Download button for results.out (uses top-level paml_results_filename) */}
                {results.paml_run_identifier && results.paml_results_filename && (
                  <button onClick={() => handleDownloadFile(results.paml_results_filename)}>
                    Download {results.paml_results_filename}
                  </button>
                )}
                {/* Download buttons for phylip and newick files using defined constants */}
                {results.paml_run_identifier && (
                    <>
                        <button onClick={() => handleDownloadFile(FILENAME_PHYLIP)}>Download PHYLIP</button>
                        <button onClick={() => handleDownloadFile(FILENAME_NEWICK)}>Download Newick Tree</button>
                    </>
                )}
              </div>
            )}
            <hr />
            <h4>Raw JSON Response from Backend:</h4>
            <pre>
              {JSON.stringify(results, null, 2)}
            </pre>
          </div>
        )}
      </main>
    </div>
  );
}

export default App;
