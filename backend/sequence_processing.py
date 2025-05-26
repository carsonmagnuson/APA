import re
import logging
from io import StringIO
from typing import List, Dict, Set, Tuple, Optional, Any
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# --- Configure Logging ---
if __name__ == '__main__': # Basic configuration for standalone testing
    logging.basicConfig(
        level=logging.DEBUG, # Set to DEBUG to see more detailed logs during testing
        format='%(asctime)s - %(levelname)s - %(name)s - %(funcName)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
logger = logging.getLogger(__name__)

# --- Sequence Processing Functions ---

def fasta_to_seqrecord(fasta_string: str) -> Dict[str, SeqRecord]:
    """
    Converts a FASTA formatted string into a dictionary of SeqRecord objects.

    The dictionary keys are sequence IDs and values are the corresponding
    SeqRecord objects from BioPython.

    Args:
        fasta_string: A string containing FASTA formatted sequence data.

    Returns:
        A dictionary mapping sequence IDs (str) to SeqRecord objects.
        Returns an empty dictionary if the input string is empty or malformed.
    """
    if not fasta_string or fasta_string.isspace():
        logger.warning("fasta_to_seqrecord received an empty input string.")
        return {}
    
    fasta_handle = StringIO(fasta_string)
    records: Dict[str, SeqRecord] = {}
    try:
        for seq_record in SeqIO.parse(fasta_handle, "fasta"):
            if seq_record.id: # Ensure there's an ID
                 records[seq_record.id] = seq_record
            else:
                logger.warning(f"Sequence found with no ID: {seq_record.description[:50]}...")
    except Exception as e:
        logger.error(f"Error parsing FASTA string: {e}", exc_info=True)
        return {} # Return empty dict on parsing error
        
    logger.info(f"Converted FASTA data to {len(records)} SeqRecord objects.")
    return records

def fasta_to_dict(fasta_string: str) -> Dict[str, str]:
    """
    Parses a FASTA formatted string into a dictionary of sequence IDs to sequences.
    
    This function specifically extracts the first part of the FASTA header as the key
    and the sequence itself as the value. It assumes a specific FASTA structure where
    the header line is immediately followed by the sequence line, and then an optional
    quality/comment line (which is skipped).

    Note: This function seems to parse FASTA differently than `fasta_to_seqrecord`.
    It was present in the original app.py.

    Args:
        fasta_string: A string containing FASTA formatted sequence data.
                      Expected format is >id description\\nSEQUENCE\\n(optional other line).

    Returns:
        A dictionary where keys are sequence IDs (str) and values are sequence strings (str).
        Returns an empty dictionary if the input string is empty or malformed.
    """
    if not fasta_string or fasta_string.isspace():
        logger.warning("fasta_to_dict received an empty input string.")
        return {}

    lines = fasta_string.splitlines()
    gene_index: Dict[str, str] = {}
    
    if not lines:
        return gene_index

    try:
        # Original code iterates by 3, implying header, sequence, and an ignored line.
        # This is a bit fragile if the FASTA format varies.
        for i in range(0, len(lines), 3): # Original code used 3
            if i + 1 < len(lines): # Ensure sequence line exists
                header = lines[i]
                sequence = lines[i+1]
                if header.startswith(">"):
                    # Takes only the first part of the header before any space
                    seq_id = header[1:].split(' ', 1)[0] 
                    gene_index[seq_id] = sequence
                else:
                    logger.warning(f"Skipping line not starting with '>': {header[:50]}...")
            else: # Header line without a sequence line
                logger.warning(f"Header line found without a subsequent sequence line: {lines[i][:50]}...")
    except IndexError:
        logger.error("IndexError during FASTA parsing in fasta_to_dict. Input format might be unexpected.", exc_info=True)
    except Exception as e:
        logger.error(f"An unexpected error occurred in fasta_to_dict: {e}", exc_info=True)
        
    logger.info(f"Parsed FASTA string into dictionary with {len(gene_index)} entries.")
    return gene_index

def get_model_organism_genes(ortholog_data: List[Dict[str, Any]]) -> List[str]: # Corrected type hint
    """
    Extracts gene IDs for model organisms from ortholog data.

    The input is expected to be a list of dictionaries, where each dictionary
    represents an organism group and contains a "genes" key with a list of gene information.
    This structure is typical for OrthoDB's 'orthologs_in_model_organisms' field.

    Args:
        ortholog_data: A list of dictionaries containing ortholog information.

    Returns:
        A list of gene ID parameters (strings) for model organism orthologs.
        Returns an empty list if the input data is malformed or no genes are found.
    """
    if not isinstance(ortholog_data, list):
        logger.warning(f"get_model_organism_genes expected a list, but got {type(ortholog_data)}.")
        return []
        
    gene_ids: List[str] = []
    for organism_group in ortholog_data:
        if isinstance(organism_group, dict) and "genes" in organism_group and isinstance(organism_group["genes"], list):
            for gene_info in organism_group["genes"]:
                try:
                    gene_id_param = gene_info.get("gene_id", {}).get("param")
                    if gene_id_param:
                        gene_ids.append(gene_id_param)
                    else:
                        logger.debug(f"Missing 'param' in gene_id for gene_info: {gene_info}")
                except AttributeError: # Should not happen with .get() but good for safety
                    logger.warning(f"Unexpected structure in gene_info: {gene_info}", exc_info=True)
        else:
            logger.debug(f"Skipping malformed organism_group: {organism_group}")
            
    logger.info(f"Extracted {len(gene_ids)} model organism gene IDs.")
    return gene_ids

def select_genes_from_seqrecord(gene_ids: List[str], records: Dict[str, SeqRecord]) -> Dict[str, SeqRecord]:
    """
    Filters a dictionary of SeqRecord objects to include only those specified by a list of gene IDs.

    Args:
        gene_ids: A list of gene IDs (strings) to select.
        records: A dictionary of SeqRecord objects, where keys are gene IDs.

    Returns:
        A new dictionary containing only the selected SeqRecord objects.
        Logs warnings for any gene IDs not found in the records.
    """
    selected_records: Dict[str, SeqRecord] = {}
    missing_ids: List[str] = []

    for gene_id in gene_ids:
        record = records.get(gene_id)
        if record:
            selected_records[gene_id] = record
        else:
            missing_ids.append(gene_id)

    if missing_ids:
        logger.warning(f"Discrepancies found: {len(missing_ids)} gene ID(s) not found in records.")
        for missed_id in missing_ids:
            logger.debug(f"Gene ID not found: {missed_id}")
            
    logger.info(f"Selected {len(selected_records)} genes from records. {len(missing_ids)} ID(s) were missing.")
    return selected_records

def pad_sequence_lengths(records: Dict[str, SeqRecord]) -> Dict[str, SeqRecord]:
    """
    Pads all sequences in a dictionary of SeqRecord objects to the length of the longest sequence.

    Padding is done by appending '-' characters to the end of shorter sequences.
    This modifies the SeqRecord objects in place.

    Args:
        records: A dictionary of SeqRecord objects, where keys are sequence IDs.

    Returns:
        The input dictionary with sequences padded (modified in-place).
        Returns an empty dictionary if the input is empty or an error occurs.
    """
    if not records:
        logger.warning("pad_sequence_lengths received an empty records dictionary.")
        return {}

    try:
        # Ensure all values are SeqRecord objects and have a .seq attribute
        valid_records = [record for record in records.values() if isinstance(record, SeqRecord) and hasattr(record, 'seq')]
        if not valid_records:
            logger.warning("No valid SeqRecord objects with sequences found for padding.")
            return records
        max_len = max(len(record.seq) for record in valid_records)
        logger.info(f"Maximum sequence length for padding is {max_len}.")
    except ValueError: 
        logger.error("Could not determine maximum sequence length. Records might be empty or invalid.", exc_info=True)
        return records 

    for record_id in records:
        if isinstance(records[record_id], SeqRecord) and hasattr(records[record_id], 'seq'):
            current_seq_str = str(records[record_id].seq)
            records[record_id].seq = Seq(current_seq_str.ljust(max_len, '-'))
        else:
            logger.warning(f"Item with key '{record_id}' is not a valid SeqRecord or lacks a sequence; skipping padding.")
        
    logger.info(f"Padded sequences in {len(records)} records to length {max_len}.")
    return records

def select_big_taxa_seq(records: Dict[str, SeqRecord]) -> List[str]:
    """
    Selects one representative sequence ID per taxon based on the longest sequence.

    Assumes sequence IDs are in the format "TAXONID:GENESPECIFICID".
    For each unique TAXONID, it identifies the GENESPECIFICID corresponding
    to the longest sequence.

    Args:
        records: A dictionary of SeqRecord objects, where keys are sequence IDs
                 (e.g., "6239_0:000672").

    Returns:
        A list of combined sequence IDs (e.g., ["6239_0:000672"]) representing the
        longest sequence for each unique taxon found. Returns an empty list if
        no valid records are processed.
    """
    if not records:
        logger.warning("select_big_taxa_seq received an empty records dictionary.")
        return []

    taxon_to_longest_info: Dict[str, Tuple[str, int]] = {}

    for seq_id, record in records.items():
        if not isinstance(record, SeqRecord) or not hasattr(record, 'seq'):
            logger.warning(f"Skipping item '{seq_id}' as it's not a valid SeqRecord with a sequence.")
            continue
        try:
            taxon_part, gene_specific_part = seq_id.split(":", 1)
        except ValueError:
            logger.warning(f"Skipping SeqRecord with ID '{seq_id}' due to unexpected format (expected 'TAXONID:GENEID').")
            continue
        
        current_len = len(record.seq)

        if taxon_part not in taxon_to_longest_info or current_len > taxon_to_longest_info[taxon_part][1]:
            taxon_to_longest_info[taxon_part] = (gene_specific_part, current_len)
            
    result_ids = [f"{taxon}:{info[0]}" for taxon, info in taxon_to_longest_info.items()]
    logger.info(f"{len(result_ids)} unique taxa found, selected longest sequence for each.")
    return result_ids


def select_biggest_k_seq(
    k: int,
    records_to_filter: Dict[str, SeqRecord],
    target_gene_id: str,
    all_records: Dict[str, SeqRecord]
) -> Dict[str, SeqRecord]:
    """
    Selects a target gene and k-1 other longest unique sequences.

    The function aims to return a dictionary containing exactly 'k' sequences:
    1. The sequence corresponding to `target_gene_id` (fetched from `all_records`).
    2. The `k-1` longest sequences from `records_to_filter`, ensuring they are
       not the `target_gene_id` and are unique.

    If `target_gene_id` is not found in `all_records`, or if `k` is not positive,
    an empty dictionary is returned with an error logged.

    If fewer than `k-1` unique sequences (excluding the target) are available
    in `records_to_filter`, the function will return the target gene plus all
    available unique other longest sequences. In this case, the total number of
    returned sequences will be less than `k`, and a warning will be logged.

    Args:
        k: The total number of sequences desired in the final selection.
        records_to_filter: A dictionary of SeqRecord objects from which to select
                           the `k-1` longest *other* sequences.
        target_gene_id: The ID of the primary gene of interest. Its sequence will be
                        sourced from `all_records`.
        all_records: A comprehensive dictionary of SeqRecord objects, which must
                     contain the `target_gene_id`.

    Returns:
        A dictionary of SeqRecord objects. Ideally, it contains `k` sequences.
        Returns an empty dictionary on critical errors (e.g., target not found, invalid k).
    """
    if k <= 0:
        logger.warning(f"k must be a positive integer, but got {k}. Returning empty dict.")
        return {}

    target_gene_record = all_records.get(target_gene_id)
    if not target_gene_record:
        logger.error(f"Target gene ID '{target_gene_id}' not found in all_records. Cannot proceed.")
        return {}

    # Prepare the list of candidate "other" sequences
    other_candidates: List[SeqRecord] = []
    for rec_id, record in records_to_filter.items():
        if rec_id != target_gene_id:
            if isinstance(record, SeqRecord) and hasattr(record, 'seq'):
                other_candidates.append(record)
            else:
                logger.debug(f"Item '{rec_id}' in records_to_filter is not a valid SeqRecord, skipping.")
    
    # Sort these other candidates by sequence length in descending order
    sorted_other_candidates = sorted(other_candidates, key=lambda x: len(x.seq), reverse=True)

    # Initialize results with the target gene
    result_records: Dict[str, SeqRecord] = {target_gene_id: target_gene_record}
    
    # Add up to k-1 longest *other* unique sequences
    num_others_to_add = k - 1
    added_others_count = 0

    if num_others_to_add > 0: # Only add others if k > 1
        for record in sorted_other_candidates:
            if added_others_count < num_others_to_add:
                # The check `record.id not in result_records` ensures uniqueness if IDs in
                # `records_to_filter` might not be unique or if an ID could somehow match target_gene_id
                # (though the initial filter `rec_id != target_gene_id` should prevent the latter).
                if record.id not in result_records: 
                    result_records[record.id] = record
                    added_others_count += 1
            else:
                # We have added the required number of other sequences
                break
    
    final_count = len(result_records)
    if final_count < k:
        logger.warning(
            f"Requested k={k} sequences, but only {final_count} unique sequences could be formed. "
            f"Target gene '{target_gene_id}' was included. "
            f"Found {added_others_count} other unique sequences from records_to_filter (needed {num_others_to_add}). "
            f"This may be due to fewer than {num_others_to_add} unique, non-target sequences available in records_to_filter."
        )
    else:
        logger.info(
            f"Successfully selected {final_count} sequences: target '{target_gene_id}' and {added_others_count} longest others."
        )
            
    return result_records


def remove_stop_codons_in_sequence(record: SeqRecord, codon_stop_array: List[str] = ["TAG", "TGA", "TAA"]) -> SeqRecord:
    """
    Removes stop codons from a single DNA sequence by replacing them with '---'.

    Iterates through the sequence in codon steps (3 bases). If a codon matches
    any in `codon_stop_array`, it's replaced. The modification happens on a copy
    of the sequence, and the SeqRecord's seq attribute is updated.

    Args:
        record: The Bio.SeqRecord object containing the sequence.
        codon_stop_array: A list of stop codon strings. Defaults to ["TAG", "TGA", "TAA"].

    Returns:
        The modified Bio.SeqRecord object.
    """
    if not isinstance(record, SeqRecord) or not hasattr(record, 'seq'):
        logger.warning(f"Invalid record passed to remove_stop_codons_in_sequence: {record}")
        return record

    # Work on a mutable version of the sequence (list of characters)
    temp_record_seq_list = list(str(record.seq))
    original_length = len(temp_record_seq_list)
    codons_removed_count = 0
    # Convert stop codons to uppercase for case-insensitive comparison
    upper_codon_stop_array = [s.upper() for s in codon_stop_array]


    for i in range(0, original_length - (original_length % 3), 3): # Ensure we only process full codons
        codon = "".join(temp_record_seq_list[i:i+3]).upper()
        if codon in upper_codon_stop_array:
            temp_record_seq_list[i:i+3] = ['-', '-', '-']
            codons_removed_count += 1
            
    if codons_removed_count > 0:
        logger.debug(f"Removed {codons_removed_count} stop codon(s) in sequence {record.id}.")

    record.seq = Seq("".join(temp_record_seq_list))
    return record

def remove_stop_codons_in_multiple(records: Dict[str, SeqRecord]) -> Dict[str, SeqRecord]:
    """
    Applies stop codon removal to all sequences in a dictionary of SeqRecord objects.

    Args:
        records: A dictionary of SeqRecord objects, where keys are sequence IDs.

    Returns:
        The dictionary with SeqRecord objects having their stop codons removed (modified in-place).
    """
    if not records:
        logger.warning("remove_stop_codons_in_multiple received an empty records dictionary.")
        return {}
        
    for record_id in records:
        if isinstance(records[record_id], SeqRecord):
            records[record_id] = remove_stop_codons_in_sequence(records[record_id])
        else:
            logger.warning(f"Item with ID '{record_id}' is not a SeqRecord; skipping stop codon removal.")
    logger.info(f"Processed {len(records)} sequences for stop codon removal.")
    return records


def convert_organism_id_to_names(records: Dict[str, SeqRecord]) -> Dict[str, SeqRecord]:
    """
    Updates the ID of each SeqRecord in a dictionary to a parsed organism name.

    The organism name is extracted from the SeqRecord's description string.
    This function attempts to parse the organism name based on common OrthoDB FASTA header formats.
    If parsing fails, it defaults to using the original ID or parts of it.
    Handles potential ID collisions by appending a counter.

    Args:
        records: A dictionary of SeqRecord objects. IDs will be modified.

    Returns:
        A new dictionary with SeqRecord objects having updated IDs. Original SeqRecords are modified.
    """
    updated_records: Dict[str, SeqRecord] = {}
    id_collision_counts: Dict[str, int] = {} # To handle cases where different original IDs map to the same new ID

    if not records:
        logger.warning("convert_organism_id_to_names received an empty records dictionary.")
        return {}

    for original_id, record in records.items():
        if not isinstance(record, SeqRecord):
            logger.warning(f"Item '{original_id}' is not a SeqRecord, skipping name conversion.")
            updated_records[original_id] = record # Or skip adding it
            continue

        organism_name = original_id # Default to original ID
        description = record.description
        
        # Attempt to parse organism name from description
        # Example OrthoDB format: ">6239_0:000672 gene=CELE_C17D12.6, organism={Caenorhabditis elegans, tax_id=6239, N2}, ..."
        # Or sometimes: ">gene_id description organism=[Genus species]"
        match_curly = re.search(r"organism=\{([^,]+)", description)
        match_square = re.search(r"organism=\[([^\]]+)\]", description)

        parsed_name_candidate = None
        if match_curly:
            parsed_name_candidate = match_curly.group(1).strip()
        elif match_square:
            parsed_name_candidate = match_square.group(1).strip()
        
        if parsed_name_candidate:
            organism_name = parsed_name_candidate
            logger.debug(f"Parsed organism name '{organism_name}' for original ID '{original_id}'.")
        else:
            # Fallback if "organism=" pattern is not found, try to use the taxon part of original_id
            if ":" in original_id:
                taxon_part = original_id.split(":")[0]
                # Further clean up taxon_part if it's like "XXXX_0"
                taxon_part_cleaned = re.sub(r'_\d+$', '', taxon_part) 
                if taxon_part_cleaned:
                    organism_name = taxon_part_cleaned
            logger.debug(f"Could not parse organism name from description for ID '{original_id}'. Using '{organism_name}'. Description: {description[:100]}...")
        
        # Sanitize the name slightly (e.g., replace spaces for better ID handling later, though PAML formatting handles more)
        # This is a light sanitation; format_ids_and_create_alignment does the strict PAML formatting.
        # organism_name = organism_name.replace(" ", "_") # Example, might not be necessary here

        # Handle potential ID collisions if multiple original IDs map to the same organism name
        final_id_candidate = organism_name
        if final_id_candidate in updated_records:
            id_collision_counts[final_id_candidate] = id_collision_counts.get(final_id_candidate, 0) + 1
            final_id_candidate = f"{organism_name}_{id_collision_counts[final_id_candidate]}"
            logger.warning(f"ID collision for '{organism_name}'. Renaming to '{final_id_candidate}' for original ID '{original_id}'.")
        
        record.id = final_id_candidate
        # Store original ID in description for traceability if it's not already there in a clear way
        if f"original_id={original_id}" not in record.description:
            record.description = f"original_id={original_id} | {description}"
        
        updated_records[final_id_candidate] = record # Add to new dict with potentially new ID
            
    logger.info(f"Converted IDs for {len(records)} sequences to organism names (or derived IDs).")
    return updated_records


def format_ids_and_create_alignment(sequences_dict: Dict[str, SeqRecord]) -> Optional[MultipleSeqAlignment]:
    """
    Formats sequence IDs for CODEML compatibility and creates a MultipleSeqAlignment object.

    CODEML typically requires sequence IDs to be strictly alphanumeric and often
    has a length limit (e.g., 10-30 characters, though 10 is safest for older versions).
    This function:
    1. Takes original IDs (which might be organism names after `convert_organism_id_to_names`).
    2. Cleans them: keeps only alphanumeric characters.
    3. Truncates to a maximum of 10 characters.
    4. Ensures uniqueness by appending a counter if duplicates arise from truncation.
    5. Creates SeqRecord objects with these new IDs and the original sequences.
    6. Returns a MultipleSeqAlignment object.

    Args:
        sequences_dict: A dictionary where keys are original sequence IDs (e.g., organism names)
                        and values are Bio.SeqRecord objects with (potentially padded) sequences.

    Returns:
        A Bio.Align.MultipleSeqAlignment object with CODEML-compatible IDs,
        or None if the input is empty or an error occurs.
    """
    if not sequences_dict:
        logger.warning("format_ids_and_create_alignment received an empty dictionary.")
        return None

    codeml_records: List[SeqRecord] = []
    generated_ids: Set[str] = set() 

    for original_id_key, seq_record_val in sequences_dict.items():
        if not isinstance(seq_record_val, SeqRecord) or not hasattr(seq_record_val, 'seq'):
            logger.warning(f"Item '{original_id_key}' is not a valid SeqRecord with a sequence. Skipping.")
            continue

        base_id = str(original_id_key) # Use the current ID of the SeqRecord
        
        alnum_id = re.sub(r'[^a-zA-Z0-9]', '', base_id)
        truncated_id = alnum_id[:10] 
        
        if not truncated_id: 
            truncated_id = "SEQ" 
            logger.debug(f"Original ID '{base_id}' resulted in empty alphanumeric ID, using '{truncated_id}'.")

        final_id = truncated_id
        counter = 1
        # Ensure unique ID, even after truncation
        while final_id in generated_ids:
            suffix = str(counter)
            base_len = 10 - len(suffix)
            if base_len <= 0: # If suffix itself is too long or base_id was empty
                 # This case should be rare if truncated_id defaults to "SEQ"
                 final_id = f"SEQ{counter}" # Ensure it's somewhat unique and short
                 if len(final_id) > 10: final_id = final_id[:10] # Re-truncate if necessary
                 # If still not unique after this, there's a deeper problem or too many identical short IDs.
                 # For extreme cases, could add more random characters, but PAML has limits.
                 if final_id in generated_ids: # Highly unlikely now
                     logger.error(f"Could not generate a unique 10-char ID for base '{base_id}'. Using potentially non-unique '{final_id}'.")
                     break # Break to avoid infinite loop, PAML might error later
            else:
                final_id = truncated_id[:base_len] + suffix
            counter += 1
        
        generated_ids.add(final_id)
        
        new_seq_rec = SeqRecord(
            Seq(str(seq_record_val.seq)), 
            id=final_id, 
            description=f"original_id={original_id_key} | {seq_record_val.description}" # Preserve original mapping
        )
        codeml_records.append(new_seq_rec)

    if not codeml_records:
        logger.warning("No records were processed for alignment after ID formatting.")
        return None
        
    logger.info(f"Formatted IDs for {len(codeml_records)} sequences for CODEML.")
    try:
        alignment = MultipleSeqAlignment(codeml_records)
        return alignment
    except Exception as e:
        logger.error(f"Error creating MultipleSeqAlignment: {e}", exc_info=True)
        return None

if __name__ == '__main__':
    logger.info("--- Sequence Processing Module Tests ---")

    fasta_example = """>seq1 description1
ATGCGTAGCATCGATCGATCG
>seq2 description2 an organism={Homo sapiens, tax_id=9606}
CGATCGATCGATCGATCGATCGA
>37653_0:000a8b gene=LOC106873942, organism={Octopus bimaculoides, tax_id=37653}
ATGCGTAGCATCGATCGATCGATGCGTAGCATCGATCGATCG
>6239_0:000672 gene=CELE_C17D12.6, organism={Caenorhabditis elegans, tax_id=6239, N2}, ogs=[{Eukaryota, 430340at2759}]
ATGTAGTGACCC---
"""
    seq_records_dict = fasta_to_seqrecord(fasta_example)
    logger.info(f"fasta_to_seqrecord output keys: {list(seq_records_dict.keys())}")

    # Test get_model_organism_genes
    mock_orthologs_list = [
        {"genes": [{"gene_id": {"param": "geneA_param"}}, {"gene_id": {"param": "geneB_param"}}]},
        {"genes": [{"gene_id": {"param": "geneC_param"}}]}
    ]
    model_genes = get_model_organism_genes(mock_orthologs_list)
    logger.info(f"get_model_organism_genes output: {model_genes}")

    # Test select_biggest_k_seq
    # Create a more diverse set for testing select_biggest_k_seq
    all_srs = {
        "target_gene": SeqRecord(Seq("GATTACAGATTACA"), id="target_gene", description="organism={Targetus organismus}"),
        "other_long1": SeqRecord(Seq("AAAAAAAAAAAAAAAAAAAA"), id="other_long1", description="organism={Longus primus}"),
        "other_long2": SeqRecord(Seq("CCCCCCCCCCCCCCCCCCCC"), id="other_long2", description="organism={Longus secundus}"),
        "other_mid1": SeqRecord(Seq("GGGGGGGGGG"), id="other_mid1", description="organism={Mediumus unus}"),
        "other_short1": SeqRecord(Seq("TTTTT"), id="other_short1", description="organism={Brevis alpha}")
    }
    records_to_filter_k = {
        "other_long1": all_srs["other_long1"],
        "other_mid1": all_srs["other_mid1"],
        "other_short1": all_srs["other_short1"],
        "target_gene": all_srs["target_gene"] # Target might also be in records_to_filter
    }

    logger.info("\n--- Testing select_biggest_k_seq (k=3) ---")
    top_3 = select_biggest_k_seq(3, records_to_filter_k, "target_gene", all_srs)
    logger.info(f"select_biggest_k_seq (k=3) output keys: {list(top_3.keys())}")
    for sid, sr in top_3.items(): logger.info(f"  {sid}: len={len(sr.seq)}")
    
    logger.info("\n--- Testing select_biggest_k_seq (k=5, more than available others) ---")
    top_5 = select_biggest_k_seq(5, records_to_filter_k, "target_gene", all_srs) # k=5, but only 3 others
    logger.info(f"select_biggest_k_seq (k=5) output keys: {list(top_5.keys())}")
    for sid, sr in top_5.items(): logger.info(f"  {sid}: len={len(sr.seq)}")

    logger.info("\n--- Testing select_biggest_k_seq (k=1, only target) ---")
    top_1 = select_biggest_k_seq(1, records_to_filter_k, "target_gene", all_srs)
    logger.info(f"select_biggest_k_seq (k=1) output keys: {list(top_1.keys())}")
    for sid, sr in top_1.items(): logger.info(f"  {sid}: len={len(sr.seq)}")


    # Test stop codon removal
    logger.info("\n--- Testing stop codon removal ---")
    if "6239_0:000672" in seq_records_dict:
        rec_with_stop = seq_records_dict["6239_0:000672"]
        logger.debug(f"Original sequence for stop codon test ({rec_with_stop.id}): {rec_with_stop.seq}")
        remove_stop_codons_in_sequence(rec_with_stop)
        logger.info(f"remove_stop_codons_in_sequence output for {rec_with_stop.id}: {rec_with_stop.seq}")
    
    # Test convert_organism_id_to_names
    logger.info("\n--- Testing convert_organism_id_to_names ---")
    # Use a fresh copy for this test
    records_for_naming = fasta_to_seqrecord(fasta_example)
    if records_for_naming:
        named_records = convert_organism_id_to_names(records_for_naming)
        logger.info(f"convert_organism_id_to_names output keys: {list(named_records.keys())}")
        for nid, n_rec in named_records.items():
            logger.info(f"  New ID: {nid}, Original ID in desc: {n_rec.description.split('|')[0].strip()}")

        # Test format_ids_and_create_alignment
        logger.info("\n--- Testing format_ids_and_create_alignment ---")
        # Pad named_records before alignment for a more realistic test
        padded_named_records = pad_sequence_lengths(named_records)
        alignment = format_ids_and_create_alignment(padded_named_records)
        if alignment:
            logger.info("Alignment created successfully after ID formatting.")
            for record in alignment:
                logger.info(f"  Aligned ID: {record.id}, Length: {len(record.seq)}, Desc: {record.description}")
        else:
            logger.error("Alignment creation failed.")
            
    logger.info("--- Sequence Processing Module Tests Finished ---")

