# INCEPTION

I was helping a professor use CODEML from the PAML package and I realized that a lot of ortholog data is available online, from publicly pullable APIs.

# PLAN
I want to make an application to automate this process such that you can enter a gene name (and potentially specify some factors) and it will run a positive selection analysis on some number of orthologs and return the omega value (non-synonymous to synonymous mutation ratio).

# CONSTITUENTS
- [ ] A Flask app
  - [ ] Pulling orthologs from OrthoDB
  - [ ] Processing orthologs into a seq and tree file
  - [ ] running CODEML positive selection analysis on data
  - [ ] serving results

- [ ] A React app 
  - [ ] Taking in a gene sequence in a specific species
  - [ ] Calling Flask API with gene sequence
  - [ ] Displaying omega value

