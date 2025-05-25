import requests
import json
import os
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.Phylo.PAML import codeml
import re # For regular expressions

def cache_data(data, name):
    path = f"backend/cached/{name}"
    if type(data) == dict:
        print("Dumping json...")
        with open(path, "w") as file:
            json.dump(data, file)
    else:
        with open(path, 'a') as f:
            print(data, file=f)
    return


def pull_OG_fasta(ortholog_group):
    if os.path.exists(f'backend/cached/{ortholog_group}.fasta'):
        print("Data found in cache...")
        return open(f"backend/cached/{ortholog_group}.fasta").read()
    print("Pulling ortholog group data...")
    url = f"https://data.orthodb.org/v12/fasta"
    params = {
            'id': ortholog_group,
            'seqtype': 'cds'
            }
    response_text = requests.get(url, params=params).text
    cache_data(response_text,f"{ortholog_group}.fasta") #cache the response so I stop maybe accidentally crashing the db
    print("Ortholog group data pulled...")
    return response_text


def fasta_to_dict(response_text):
    response_text = response_text.split_lines()
    gene_index = {response_text[index].split(' ', 1)[0]: response_text[index+1] for index in range(0, len(response_text), 3)}
    return gene_index


def pull_model_organism_orthologs(gene):
    if os.path.exists(f"backend/cached/{gene}_search.json"):
        print("Data found in cache...")
        with open(f'backend/cached/{gene}_search.json', 'r') as file:
            response_json = json.load(file)
    else:
        gene_search = f"https://data.orthodb.org/v12/genesearch?query={gene}"
        response_json = requests.get(gene_search).json()
        cache_data(response_json,f"{gene}_search.json") #cache the response so I stop maybe accidentally crashing the db

    orthologs = response_json["orthologs_in_model_organisms"]
    ortholog_count = sum(sum(1 for found in organism["genes"]) for organism in orthologs)
    print(f"{ortholog_count} orthologs found in model organisms...")
    return orthologs

def fasta_to_seqrecord(response_text):
    fasta_handle = StringIO(response_text)
    records = {seq.id: seq for seq in SeqIO.parse(fasta_handle, "fasta")}

    print(f"Fasta data converted to seqrecord dict...")
    return records

def select_genes_from_seqrecord(genes, records):
    missing = []
    result = {}
    print("Paring down records for selected genes...")

    for gene in genes:
        if records.get(gene):
            result[gene] = records[gene]
        else:
            missing.append(gene)

    if missing:
        print(f"WARNING: Discrepancies found between records and selected genes...")
        for missed in missing:
            print(f"Gene id: {missed} not found in records...")

    # result = {gene: records[gene] for gene in genes}

    return result

def get_model_organism_genes(orthologs):
    return [gene_id["gene_id"]["param"] for organism in orthologs for gene_id in organism["genes"]]

def pad_sequence_lengths(records):
    big = max(len(record.seq) for record in records.values())
    print(f"Maximum sequence length is {big}...")

    for record in records.keys():
        records[record].seq = Seq(str(records[record].seq).ljust(big, '-'))
    return records

def select_big_taxa_seq(records):
    result = {taxa: gene_id for taxa, gene_id in [key.split(":") for key in records.keys()]} ##initialize a dict to hold the biggest boyos with whatever default id
    for taxa_key in records.keys():
        taxa = taxa_key.split(":")[0]
        result[taxa] = max([result[taxa], taxa_key.split(":")[1]], key=lambda x: len(records[taxa + ":" + x].seq)) ##find the bigger of the two sequences and replace with that id belonging to that sequence and taxa
    print(f"{len(result)} unique taxa found...")
    return list(key + ":" + value for key, value in result.items())

def select_biggest_k_seq(k, records, gene_id, base_records):
    #todo: add error checking to make sure k is smaller than or equal to records length
    top_k = sorted(list(records.values()), key=lambda x: len(x.seq))[-(k-1):] #don't forget to append the selected gene here later
    converted_dict = {seq.id: seq for seq in top_k}
    if gene_id in records.keys():
        print("FOUND GENE IN RECORDS")
    converted_dict[gene_id] = base_records[gene_id]
    return converted_dict

def remove_stop_codons_in_sequence(record, codon_stop_array = ["TAG", "TGA", "TAA"]): ## thanks to au_ndh at https://www.biostars.org/p/296261/ for this def
    tempRecordSeq = list(record.seq)
    for index in range(0, len(record.seq), 3):
            codon = record.seq[index:index+3]
            if codon in codon_stop_array:
                tempRecordSeq[index:index+3] = '-','-','-'
                print(f"removed a stop codon in {record.id}")
    record.seq = Seq("".join(tempRecordSeq))
    return record

def remove_stop_codons_in_multiple(records):
    result = {key:remove_stop_codons_in_sequence(record) for key, record in records.items()}
    return result

def construct_phylo_tree(alignment, rooted=False):
    print("Constructing phylogenetic tree...")
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix) if rooted else constructor.nj(distance_matrix)
    return tree


def write_newick_tree_with_header(tree, filepath):
    """
    Writes the phylogenetic tree to a file in Newick format,
    prefixed with a PAML-style header (num_taxa num_trees).

    Args:
        tree (Bio.Phylo.BaseTree.Tree): The tree object to write.
        filepath (str): The path to the output file.
    """
    if not tree:
        print(f"Warning: Tree object is None. Cannot write Newick file to {filepath}.")
        return

    try:
        num_taxa = len(tree.get_terminals())  # Get the number of leaves (taxa)
    except Exception as e:
        print(f"Warning: Could not determine number of taxa from tree: {e}. Cannot write Newick file to {filepath}.")
        return

    with open(filepath, 'w') as f_tree:
        # Write the header line (number of taxa, 1 for one tree in the file)
        # PAML examples sometimes show a leading space, sometimes not.
        # Let's go with no leading space for this line as it's typical.
        f_tree.write(f"{num_taxa} 1\n") 
        
        # Write the tree using Phylo.write to the same file handle
        # This will append the Newick string after the header
        Phylo.write(tree, f_tree, "newick")
            
    print(f"Newick tree with header written to '{filepath}'.")

def convert_organism_id_to_names(records):
    for gene_id in records.keys():
        print(records[gene_id])
        organism_desc_str = records[gene_id].description
        organism_name = organism_desc_str[organism_desc_str.index("{"):].split(",")[4].split(":")[1][1:-1]
        print(organism_name)
        current_id = records[gene_id].id
        print(current_id)
        records[gene_id].id = organism_name
    return records

# Ensure these imports are at the top of your app.py
import os
import shutil
from Bio.Phylo.PAML import codeml # Should already be there

def run_codeml_positive_selection(tree_path, seqfile_path, 
                                  output_dir="paml_output", 
                                  results_filename="results.out",
                                  cleanup_working_dir=True):
    """
    Runs CODEML to test for positive selection using site models M0, M1a, and M2a.

    Args:
        tree_path (str): Path to the Newick tree file.
        seqfile_path (str): Path to the PHYLIP sequence alignment file.
        output_dir (str): Directory where CODEML will run and store output.
                          This directory will be created if it doesn't exist.
        results_filename (str): Name of the main CODEML output file.
        cleanup_working_dir (bool): If True, removes the output_dir after successful parsing.

    Returns:
        dict: A dictionary containing parsed results for models M0, M1a, and M2a,
              including lnL values and key parameter estimates. Returns None if
              CODEML run or parsing fails.
    """
    cml = codeml.Codeml()

    # Ensure absolute paths are used for PAML, as it can be sensitive to relative paths
    # when its working directory is different from where the script is run.
    cml.alignment = os.path.abspath(seqfile_path)
    cml.tree = os.path.abspath(tree_path)
    cml.out_file = results_filename # This will be inside the working_dir

    # Manage working directory
    # PAML creates a lot of intermediate files. It's best to run it in a dedicated directory.
    abs_output_dir = os.path.abspath(output_dir)
    if not os.path.exists(abs_output_dir):
        os.makedirs(abs_output_dir)
    cml.working_dir = abs_output_dir

    print(f"Running CODEML in directory: {abs_output_dir}")
    print(f"Alignment file: {cml.alignment}")
    print(f"Tree file: {cml.tree}")

    # --- Set CODEML options ---
    # These options are tailored for comparing site models M0, M1a, and M2a
    # to test for positive selection.
    try:
        cml.set_options(
            # --- Input/Output & General ---
            noisy=0,          # 0-9: Amount of screen output. 0 or 1 is usually fine.
            verbose=1,        # 0: concise; 1: detailed output in results.out.
            runmode=0,        # 0: user tree. We provide the tree.
            
            # --- Data Type ---
            seqtype=1,        # 1: codon sequences. Essential.
            CodonFreq=2,      # 2: F3x4 model (codon frequencies from nt freqs at 3 positions). Good general choice.
            
            # --- Site Models for Positive Selection (NSsites) ---
            # We will run models M0 (one-ratio), M1a (neutral), and M2a (selection).
            # This combination allows for likelihood ratio tests (LRTs) for positive selection.
            model=0,          # This is a general setting often used with NSsites.
                              # When model=0, NSsites defines specific site models.
            NSsites=[0, 1, 2],# 0: M0 (one dN/dS ratio for all sites)
                              # 1: M1a (neutral: nearly neutral (0 < w < 1) and neutral (w = 1) sites)
                              # 2: M2a (selection: extends M1a with a class for w > 1)
            
            # --- Substitution Model Parameters ---
            icode=0,          # 0: universal genetic code. Change if needed.
            fix_kappa=0,      # 0: estimate kappa (transition/transversion ratio). Generally recommended.
            kappa=2.0,        # Initial guess for kappa if fix_kappa=0.
            fix_omega=0,      # 0: estimate omega. Essential for these models.
            omega=0.5,        # Initial guess for omega (for M0, and starting point for others).
            
            fix_alpha=1,      # 0: estimate alpha (gamma shape for site rates); 1: fix alpha.
                              # For M0, M1a, M2a, alpha is typically fixed at 0 (no continuous gamma heterogeneity
                              # beyond what's defined by the discrete site classes).
            alpha=0.0,        # Value of alpha if fix_alpha=1. (alpha=0 means no gamma)
            
            # --- Other Important Settings ---
            cleandata=1,      # 1: Remove sites with alignment gaps or ambiguous characters. Recommended.
            Small_Diff=1e-7,  # A small number for iteration convergence. Default is often fine.
            getSE=0           # 0: Don't calculate Standard Errors (faster). Set to 1 if SEs are needed.
                              # RateAncestor=0 is default & fine (no ancestral state reconstruction).
        )
    except ValueError as e:
        print(f"Error setting PAML options: {e}")
        return None

    # --- Run CODEML ---
    try:
        print("Starting CODEML run...")
        # The Biopython Codeml object calls the codeml executable.
        # Ensure 'codeml' is in your system PATH or PAML_PATH env variable is set.
        results = cml.run(verbose=True) # verbose=True here prints Biopython's call to codeml
        print("CODEML run finished.")
        
        # --- Parse Results ---
        # Biopython parses results["NSsites"] into a dictionary where keys are
        # the model numbers (0 for M0, 1 for M1a, 2 for M2a from NSsites list).
        parsed_results = {}
        if "NSsites" in results:
            for model_key, model_results in results["NSsites"].items():
                parsed_results[f"M{model_key}"] = {
                    "lnL": model_results.get("lnL"),
                    "parameters": model_results.get("parameters")
                }
            print("Successfully parsed results.")
        else:
            print("Warning: 'NSsites' key not found in PAML results. Parsing might be incomplete.")
            print("Full PAML results object:", results) # Print the whole dict for debugging

        # --- Cleanup (Optional) ---
        if cleanup_working_dir:
            try:
                shutil.rmtree(abs_output_dir)
                print(f"Cleaned up working directory: {abs_output_dir}")
            except OSError as e:
                print(f"Error cleaning up working directory {abs_output_dir}: {e}")
        
        return parsed_results

    except Exception as e: # Catching generic Exception from cml.run() or parsing
        print(f"An error occurred during CODEML execution or parsing: {e}")
        print(f"Check files in the working directory for more details: {abs_output_dir}")
        # It's useful to inspect results.out and other files in abs_output_dir manually
        # if an error occurs here.
        return None

def format_ids_and_create_alignment(padded_sequences_dict):
    """
    Formats sequence IDs for CODEML compatibility (strictly alphanumeric, max 10 chars)
    and creates a MultipleSeqAlignment object.

    Args:
        padded_sequences_dict (dict): A dictionary where keys are original sequence IDs
                                     and values are Bio.SeqRecord objects with padded sequences.

    Returns:
        Bio.Align.MultipleSeqAlignment: An alignment object with formatted IDs.
    """
    codeml_records = []
    # Using a list to store final IDs to check for uniqueness after modification/truncation
    # This is more robust than a set if order matters for tie-breaking, though for IDs it's usually about uniqueness.
    # For PAML, IDs just need to be unique and correctly formatted.
    generated_ids = [] 

    for original_id_key, seq_record_val in padded_sequences_dict.items():
        # 1. Convert original ID to string
        base_id = str(original_id_key)
        
        # 2. Keep only alphanumeric characters
        alnum_id = re.sub(r'[^a-zA-Z0-9]', '', base_id)
        
        # 3. Truncate to a maximum of 10 characters
        truncated_id = alnum_id[:10]
        
        # 4. Handle empty IDs after cleaning (e.g., if original was all symbols)
        if not truncated_id:
            truncated_id = "SEQ" # Default base

        # 5. Ensure uniqueness by appending a counter if duplicates arise from truncation
        final_id = truncated_id
        counter = 1
        while final_id in generated_ids:
            suffix = str(counter)
            # Adjust length of base to make space for suffix, ensuring total is not > 10
            if len(truncated_id) + len(suffix) > 10:
                base_len = 10 - len(suffix)
                final_id = truncated_id[:base_len] + suffix
            else:
                final_id = truncated_id + suffix
            counter += 1
            # Safety break if it somehow loops too much making unique IDs (highly unlikely for 10 seqs)
            if counter > len(padded_sequences_dict) + 5: 
                # This indicates a deeper issue or very non-diverse starting IDs
                # For now, we'll just use what we have and PAML might complain about non-unique IDs
                # Or, you could raise an error.
                print(f"Warning: Could not generate a unique ID for base {truncated_id} within reasonable attempts.")
                break 
        
        generated_ids.append(final_id)
        
        current_sequence = Seq(str(seq_record_val.seq))
        new_seq_rec = SeqRecord(current_sequence, id=final_id, description="")
        codeml_records.append(new_seq_rec)

    return MultipleSeqAlignment(codeml_records)

# automated-phylogenetic-analysis/backend/app.py
# Add this function to your app.py file.
# Ensure other necessary imports like MultipleSeqAlignment are already at the top of app.py
# from Bio.Align import MultipleSeqAlignment # Should already be there or imported as needed

# ... (other functions in your app.py like pull_OG_fasta, format_ids_and_create_alignment, etc.) ...

def write_phylip_manual(alignment, filepath):
    """
    Manually writes the alignment to a file in PHYLIP sequential format,
    ensuring correct spacing for PAML/CODEML.

    Args:
        alignment (Bio.Align.MultipleSeqAlignment): The alignment object to write.
                                                   It's assumed that the record IDs within
                                                   this alignment object have already been
                                                   cleaned/formatted (e.g., by
                                                   format_ids_and_create_alignment).
        filepath (str): The path to the output file.
    """
    with open(filepath, 'w') as f:
        if not alignment or len(alignment) == 0: # Check if alignment is None or empty
            print(f"Warning: Alignment is empty or None. Cannot write PHYLIP file to {filepath}.")
            # Optionally write an empty valid PHYLIP header if that's preferred for downstream tools
            # f.write(" 0 0\n")
            return

        num_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()

        # Write the header line: e.g., " 10 12531"
        f.write(f" {num_sequences} {alignment_length}\n")

        for record in alignment:
            # IDs should already be cleaned by format_ids_and_create_alignment to be <=10 chars
            # Pad with spaces on the right to ensure it's exactly 10 characters.
            name = record.id.ljust(10) 
            
            sequence_string = str(record.seq)

            # Write the 10-character name, two spaces, then the sequence.
            f.write(f"{name}  {sequence_string}\n")
            
    print(f"PHYLIP file '{filepath}' manually written in sequential format.")


def main():
    gene = "WBGene00004963"
    gene_id = "6239_0:000672"
    print("Starting up")
    orthologs = pull_model_organism_orthologs(gene)
    records = fasta_to_seqrecord(pull_OG_fasta('430340at2759'))
    print(list(records.items())[0][1])
    model_organism_genes = get_model_organism_genes(orthologs)
    model_sequences = select_genes_from_seqrecord(model_organism_genes, records)
    taxa_sequences = select_genes_from_seqrecord(select_big_taxa_seq(model_sequences), model_sequences)
    top_10 = select_biggest_k_seq(10, taxa_sequences, gene_id, records)
    padded_sequences = pad_sequence_lengths(top_10)
    padded_and_removed_sequences = remove_stop_codons_in_multiple(padded_sequences)
    padded_removed_and_stripped_sequences = convert_organism_id_to_names(padded_and_removed_sequences)
    align = MultipleSeqAlignment(list(padded_and_removed_sequences.values()))
    print(align)
    with open('output.phylip', 'a') as f:
        print(format(align, "phylip"), file=f)

    tree = construct_phylo_tree(align)

    Phylo.write(tree, 'tree.nwk', "newick")
    # run_CODEML_analysis('tree.nwk', 'output.phylip')
    # freq = {}
    # for gene in taxa_sequences.values():
    #     length = len(gene.seq)
    #     freq[length] = freq.get(length, 0) + 1
    #
    # print(len(records["6239_0:000672"].seq))
    # print((sorted(freq.items(), key=lambda x: x[1])[-1]))
