# INCEPTION

I was helping a professor use CODEML from the PAML package and I realized that a lot of ortholog data is available online, from publicly pullable APIs.

# PLAN
I want to make an application to automate this process such that you can enter a gene name (and potentially specify some factors) and it will run a positive selection analysis on some number of orthologs and return the omega value (non-synonymous to synonymous mutation ratio).

# DELIVERABLES
- [X] ~~A Flask app~~ fastAPI app
  - [X] Pulling orthologs from OrthoDB
  - [X] Processing orthologs into a seq and tree file
    - [X] Sequence selection/wrangling/aligning 
    - [X] Sequence file export
    - [X] Phylogenetic tree generation
    - [X] Tree file export
    - ~~[ ] Control file export~~
  - [X] Running CODEML positive selection analysis on data
  - [X] Serving results

- [ ] A React app 
  - [X] Taking in a gene ID in a specific species
  - [X] Calling backend API with gene sequence
  - [X] Allowing downloadable phylip, tree and results files
  - [ ] Presenting relevant data from analysis results

