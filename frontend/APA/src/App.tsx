// src/App.tsx
import React, { useState } from 'react';
import axios from 'axios';
import './App.css';

interface CodemlAnalysisSummary {
  status?: string;
  message?: string;
  results_file_path_on_server?: string; // Path on the server (for logging/debug)
  results_filename?: string;           // e.g., "results.out"
  error_details?: string;
  M0?: { lnL?: number; parameters?: any };
  M1a?: { lnL?: number; parameters?: any };
  M2a?: { lnL?: number; parameters?: any };
  // Add other parsed fields if necessary
}

interface AnalysisResult {
  input_ncbi_id?: string;
  target_gene_orthodb_id?: string;
  ortholog_group_id_used?: string;
  ortho_level_tax_id_used?: string;
  k_sequences_requested?: number;
  sequences_in_alignment?: number;
  paml_run_identifier?: string; // NEW: Identifier for the run directory
  paml_results_filename?: string; // NEW: Filename like 'results.out'
  codeml_analysis_summary?: CodemlAnalysisSummary;
}

function App() {
  const [ncbiId, setNcbiId] = useState<string>('');
  const [orthoLevel, setOrthoLevel] = useState<string>('2759'); // Default to Eukaryota
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

    try {
      const response = await axios.post(`${BACKEND_URL}/analyze/`, {
        ncbi_id: ncbiId,
        ortho_level_tax_id: orthoLevel,
        // k_sequences: 10, // You can make this an input too
        // cleanup_paml_dir: false, // For this download strategy, backend controls cleanup or it's off
      });
      setResults(response.data);
    } catch (err: any) {
      if (axios.isAxiosError(err) && err.response) {
        setError(err.response.data.detail || 'An error occurred during analysis.');
        console.error('Analysis error:', err.response.data);
      } else {
        setError('An unexpected error occurred.');
        console.error('Unexpected error:', err);
      }
    } finally {
      setIsLoading(false);
    }
  };

  const handleDownloadFile = (filenameToDownload: string | undefined) => {
    if (results?.paml_run_identifier && filenameToDownload) {
      const downloadUrl = `${BACKEND_URL}/download_paml_output/${results.paml_run_identifier}/${filenameToDownload}`;
      window.open(downloadUrl, '_blank');
    } else {
      setError(`Cannot download ${filenameToDownload || 'file'}: missing run identifier or filename.`);
    }
  };

  return (
    <div className="App">
      <header className="App-header">
        <h1>Phylogenetic Analysis</h1>
      </header>
      <main>
        <div>
          <label htmlFor="ncbiId">NCBI Gene ID: </label>
          <input
            type="text"
            id="ncbiId"
            value={ncbiId}
            onChange={(e) => setNcbiId(e.target.value)}
            placeholder="e.g., 173042"
          />
        </div>
        <div>
          <label htmlFor="orthoLevel">Orthology Level (NCBI TaxID): </label>
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

        {error && <p style={{ color: 'red' }}>Error: {error}</p>}

        {results && (
          <div className="results">
            <h2>Analysis Results for NCBI ID: {results.input_ncbi_id}</h2>
            <p>Target OrthoDB ID: {results.target_gene_orthodb_id}</p>
            <p>Ortholog Group Used: {results.ortholog_group_id_used} (Level: {results.ortho_level_tax_id_used})</p>
            <p>Sequences Requested: {results.k_sequences_requested}</p>
            <p>Sequences in Final Alignment: {results.sequences_in_alignment}</p>
            <p>PAML Run Identifier: {results.paml_run_identifier}</p>
            
            {results.codeml_analysis_summary && (
              <div>
                <h3>CODEML Analysis Status:</h3>
                <p>Status: {results.codeml_analysis_summary.status || 'N/A'}</p>
                <p>Message: {results.codeml_analysis_summary.message || 'N/A'}</p>
                {results.codeml_analysis_summary.error_details && (
                  <p>Error Details: {results.codeml_analysis_summary.error_details}</p>
                )}

                {/* Display parsed CODEML results if available */}
                {results.codeml_analysis_summary.M0?.lnL !== undefined && (
                    <p>M0 lnL: {results.codeml_analysis_summary.M0.lnL}</p>
                )}
                {results.codeml_analysis_summary.M1a?.lnL !== undefined && (
                    <p>M1a lnL: {results.codeml_analysis_summary.M1a.lnL}</p>
                )}
                 {results.codeml_analysis_summary.M2a?.lnL !== undefined && (
                    <p>M2a lnL: {results.codeml_analysis_summary.M2a.lnL}</p>
                )}

                {results.paml_run_identifier && results.codeml_analysis_summary.results_filename && (
                  <button onClick={() => handleDownloadFile(results.codeml_analysis_summary?.results_filename)}>
                    Download {results.codeml_analysis_summary.results_filename}
                  </button>
                )}
                {/* Add buttons for other files if needed */}
                {results.paml_run_identifier && (
                    <>
                        <button onClick={() => handleDownloadFile(DEFAULT_OUTPUT_PHYLIP_FILENAME)}>Download PHYLIP</button>
                        <button onClick={() => handleDownloadFile(DEFAULT_OUTPUT_NEWICK_FILENAME)}>Download Newick Tree</button>
                    </>
                )}
              </div>
            )}
            <hr />
            <h4>Raw JSON Response from Backend:</h4>
            <pre style={{ textAlign: 'left', backgroundColor: '#f0f0f0', padding: '10px', overflowX: 'auto', maxHeight: '200px' }}>
              {JSON.stringify(results, null, 2)}
            </pre>
          </div>
        )}
      </main>
    </div>
  );
}

export default App;
