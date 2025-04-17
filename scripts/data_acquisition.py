#!/usr/bin/env python3
"""
Data acquisition module for fetching protein sequences and domain information.
"""

import os
import time
import logging
import requests
from typing import List, Dict, Any, Optional, Tuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess
import tempfile
import io  # Add this import for StringIO

def parse_seed_proteins(seed_file_path: str) -> List[str]:
    """
    Parse seed protein UniProt IDs from a file.
    
    Args:
        seed_file_path: Path to the file containing seed UniProt IDs.
        
    Returns:
        List of UniProt IDs.
    """
    uniprot_ids = []
    
    with open(seed_file_path, 'r') as f:
        for line in f:
            # Skip comment lines and empty lines
            if line.strip() and not line.startswith('#'):
                # Take the first field, which should be the UniProt ID
                uniprot_id = line.strip().split('\t')[0]
                uniprot_ids.append(uniprot_id)
                
    logging.info(f"Parsed {len(uniprot_ids)} UniProt IDs from {seed_file_path}")
    return uniprot_ids

def fetch_uniprot_sequence(uniprot_id: str, retry_count: int = 3, retry_delay: int = 1) -> Optional[str]:
    """
    Fetch protein sequence from UniProt API.
    
    Args:
        uniprot_id: UniProt ID of the protein.
        retry_count: Number of times to retry the request if it fails.
        retry_delay: Delay in seconds between retries.
        
    Returns:
        Protein sequence as a string, or None if the request fails.
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    
    for attempt in range(retry_count):
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise exception for HTTP errors
            
            # Parse the FASTA response
            sequence_lines = response.text.strip().split('\n')
            # Skip the header line and concatenate the sequence lines
            sequence = ''.join(sequence_lines[1:])
            
            logging.info(f"Successfully fetched sequence for {uniprot_id}")
            return sequence
            
        except requests.exceptions.RequestException as e:
            logging.warning(f"Failed to fetch sequence for {uniprot_id}, attempt {attempt+1}/{retry_count}: {e}")
            if attempt < retry_count - 1:
                time.sleep(retry_delay)
    
    logging.error(f"Failed to fetch sequence for {uniprot_id} after {retry_count} attempts")
    return None

def fetch_uniprot_domains(uniprot_id: str, retry_count: int = 3, retry_delay: int = 1) -> Optional[List[Dict[str, Any]]]:
    """
    Fetch domain information for a protein from UniProt API.
    
    Args:
        uniprot_id: UniProt ID of the protein.
        retry_count: Number of times to retry the request if it fails.
        retry_delay: Delay in seconds between retries.
        
    Returns:
        List of domain information dictionaries, or None if the request fails.
    """
    url = f"https://www.ebi.ac.uk/proteins/api/features/{uniprot_id}"
    headers = {"Accept": "application/json"}
    
    for attempt in range(retry_count):
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()  # Raise exception for HTTP errors
            
            data = response.json()
            # Extract domain features (type: DOMAIN)
            domains = []
            if 'features' in data:
                for feature in data['features']:
                    if feature.get('type') == 'DOMAIN':
                        domain_info = {
                            'domain_id': feature.get('description', ''),
                            'start': feature['begin'],
                            'end': feature['end']
                        }
                        logging.debug(f"Found domain in {uniprot_id}: {domain_info}")
                        domains.append(domain_info)
            
            logging.info(f"Successfully fetched {len(domains)} domains for {uniprot_id}")
            return domains
            
        except requests.exceptions.RequestException as e:
            logging.warning(f"Failed to fetch domains for {uniprot_id}, attempt {attempt+1}/{retry_count}: {e}")
            if attempt < retry_count - 1:
                time.sleep(retry_delay)
    
    logging.error(f"Failed to fetch domains for {uniprot_id} after {retry_count} attempts")
    return None

def fetch_seed_proteins_data(uniprot_ids: List[str]) -> Tuple[Dict[str, str], Dict[str, List[Dict[str, Any]]]]:
    """
    Fetch sequences and domain information for seed proteins.
    
    Args:
        uniprot_ids: List of UniProt IDs.
        
    Returns:
        Tuple containing:
        - Dictionary mapping UniProt IDs to their sequences.
        - Dictionary mapping UniProt IDs to their domain information.
    """
    sequences = {}
    domains = {}
    
    for uniprot_id in uniprot_ids:
        # Fetch sequence
        sequence = fetch_uniprot_sequence(uniprot_id)
        if sequence:
            sequences[uniprot_id] = sequence
        
        # Fetch domains
        domain_info = fetch_uniprot_domains(uniprot_id)
        if domain_info:
            domains[uniprot_id] = domain_info
    
    logging.info(f"Fetched sequences for {len(sequences)}/{len(uniprot_ids)} proteins")
    logging.info(f"Fetched domain information for {len(domains)}/{len(uniprot_ids)} proteins")
    
    return sequences, domains

def save_sequences_to_fasta(sequences: Dict[str, str], output_file: str) -> bool:
    """
    Save protein sequences to a FASTA file.
    
    Args:
        sequences: Dictionary mapping sequence IDs to their sequences.
        output_file: Path to the output FASTA file.
        
    Returns:
        True if the sequences were saved successfully, False otherwise.
    """
    records = []
    
    for seq_id, sequence in sequences.items():
        record = SeqRecord(
            Seq(sequence),
            id=seq_id,
            description=""
        )
        records.append(record)
    
    try:
        SeqIO.write(records, output_file, "fasta")
        logging.info(f"Successfully saved {len(records)} sequences to {output_file}")
        return True
    except Exception as e:
        logging.error(f"Failed to save sequences to {output_file}: {e}")
        return False

def filter_cat_domains(domains: Dict[str, List[Dict[str, Any]]]) -> Dict[str, List[Dict[str, Any]]]:
    """
    Filter domain information to keep only C, A, and T domains.
    If only T domains are found, infer C and A domains based on NRPS module architecture.
    
    Args:
        domains: Dictionary mapping UniProt IDs to their domain information.
        
    Returns:
        Dictionary with filtered domain information.
    """
    cat_domains = {}
    
    # Define both Pfam IDs and keyword patterns for domain types
    domain_patterns = {
        "C": ["PF00668", "condensation", "condens", "c domain"],
        "A": ["PF00501", "adenyl", "amp-binding", "a domain", "amp bind"],
        "T": ["PF00550", "pp-binding", "carrier", "thiolation", "phosphopantetheine", "t domain", "acp", "pcp"]
    }
    
    for uniprot_id, domain_list in domains.items():
        cat_domains[uniprot_id] = []
        
        # First, identify all T domains (carriers)
        t_domains = []
        for domain in domain_list:
            domain_description = domain.get('domain_id', '').lower()
            
            # Try to find T domains
            for pattern in domain_patterns["T"]:
                if pattern.lower() in domain_description:
                    t_domain = {
                        'domain_type': "T",
                        'domain_id': pattern,
                        'start': int(domain['start']),
                        'end': int(domain['end'])
                    }
                    t_domains.append(t_domain)
                    logging.debug(f"Matched T domain in {uniprot_id} (pattern '{pattern}' matched '{domain_description}'): {domain['start']}-{domain['end']}")
                    break
        
        # Next, look for C and A domains
        c_domains = []
        a_domains = []
        for domain in domain_list:
            domain_description = domain.get('domain_id', '').lower()
            
            # Try to find C domains
            for pattern in domain_patterns["C"]:
                if pattern.lower() in domain_description:
                    c_domain = {
                        'domain_type': "C",
                        'domain_id': pattern,
                        'start': int(domain['start']),
                        'end': int(domain['end'])
                    }
                    c_domains.append(c_domain)
                    logging.debug(f"Matched C domain in {uniprot_id} (pattern '{pattern}' matched '{domain_description}'): {domain['start']}-{domain['end']}")
                    break
            
            # Try to find A domains
            for pattern in domain_patterns["A"]:
                if pattern.lower() in domain_description:
                    a_domain = {
                        'domain_type': "A",
                        'domain_id': pattern,
                        'start': int(domain['start']),
                        'end': int(domain['end'])
                    }
                    a_domains.append(a_domain)
                    logging.debug(f"Matched A domain in {uniprot_id} (pattern '{pattern}' matched '{domain_description}'): {domain['start']}-{domain['end']}")
                    break
        
        # If we have T domains but no or incomplete C/A domains, infer missing domains
        if t_domains and (len(c_domains) < len(t_domains) or len(a_domains) < len(t_domains)):
            logging.info(f"Inferring C and A domains for {uniprot_id} based on {len(t_domains)} T domain positions")
            
            # Sort T domains by position
            t_domains.sort(key=lambda x: x['start'])
            
            # Add all found domains to our filtered list
            cat_domains[uniprot_id].extend(t_domains)
            cat_domains[uniprot_id].extend(c_domains)
            cat_domains[uniprot_id].extend(a_domains)
            
            # For each T domain, create inferred C and A domains if they don't already exist
            for idx, t_domain in enumerate(t_domains):
                t_start = t_domain['start']
                t_end = t_domain['end']
                t_length = t_end - t_start + 1
                
                # Typical domain sizes and gaps in NRPS modules:
                # C domains: ~450 aa
                # A domains: ~550 aa
                # T domains: ~80-100 aa
                
                # Check if we need to infer an A domain for this T domain
                a_domain_exists = False
                for a_domain in a_domains:
                    # If an A domain is close to this T domain, consider it part of the same module
                    if abs(a_domain['end'] - t_start) < 100:
                        a_domain_exists = True
                        break
                
                if not a_domain_exists:
                    # Infer A domain position (just before T domain)
                    a_end = t_start - 1  # Directly before T domain
                    a_start = max(1, a_end - 550)  # Typical A domain size
                    
                    inferred_a_domain = {
                        'domain_type': "A",
                        'domain_id': "inferred_adenylation",
                        'start': a_start,
                        'end': a_end
                    }
                    cat_domains[uniprot_id].append(inferred_a_domain)
                    logging.debug(f"Inferred A domain for {uniprot_id} module {idx+1}: {a_start}-{a_end}")
                
                # Check if we need to infer a C domain for this module
                c_domain_exists = False
                for c_domain in c_domains:
                    # If a C domain is close to the A domain we found/inferred, consider it part of the same module
                    if a_domain_exists:
                        for a_domain in a_domains:
                            if abs(c_domain['end'] - a_domain['start']) < 100:
                                c_domain_exists = True
                                break
                    else:
                        # If we inferred an A domain, check if C domain is close to inferred A
                        if abs(c_domain['end'] - a_start) < 100:
                            c_domain_exists = True
                            break
                
                if not c_domain_exists:
                    # Infer C domain position (before A domain)
                    c_end = a_start - 1  # Directly before inferred A domain
                    c_start = max(1, c_end - 450)  # Typical C domain size
                    
                    # Special case: for first module, C domain might not be present
                    # Only infer C domain if this isn't the first module or if position makes sense
                    if idx > 0 or c_start > 50:  # Avoid creating C domains at the very beginning
                        inferred_c_domain = {
                            'domain_type': "C",
                            'domain_id': "inferred_condensation",
                            'start': c_start,
                            'end': c_end
                        }
                        cat_domains[uniprot_id].append(inferred_c_domain)
                        logging.debug(f"Inferred C domain for {uniprot_id} module {idx+1}: {c_start}-{c_end}")
        else:
            # If no T domains found or if we already have complete C-A-T sets, add all detected domains
            for domain in domain_list:
                domain_description = domain.get('domain_id', '').lower()
                
                # Try to find matching pattern for domain
                matched_type = None
                for domain_type, patterns in domain_patterns.items():
                    for pattern in patterns:
                        if pattern.lower() in domain_description:
                            matched_type = domain_type
                            # Create a new domain dict with domain_type field
                            cat_domain = {
                                'domain_type': domain_type,
                                'domain_id': pattern,
                                'start': int(domain['start']),
                                'end': int(domain['end'])
                            }
                            cat_domains[uniprot_id].append(cat_domain)
                            logging.debug(f"Added {domain_type} domain in {uniprot_id}: {domain['start']}-{domain['end']}")
                            break
                    if matched_type:
                        break
                
                if not matched_type:
                    logging.debug(f"Unknown domain type in {uniprot_id}: {domain_description}")
        
        # Sort domains by start position
        cat_domains[uniprot_id] = sorted(cat_domains[uniprot_id], key=lambda x: x['start'])
        
        # Log the number of C, A, T domains found
        c_domains = sum(1 for d in cat_domains[uniprot_id] if d['domain_type'] == 'C')
        a_domains = sum(1 for d in cat_domains[uniprot_id] if d['domain_type'] == 'A')
        t_domains = sum(1 for d in cat_domains[uniprot_id] if d['domain_type'] == 'T')
        logging.info(f"Final domains in {uniprot_id}: {c_domains} C domains, {a_domains} A domains, {t_domains} T domains")
    
    return cat_domains

def retrieve_sequences_from_fasta(protein_ids: List[str], fasta_file: str, 
                                 chunk_size: int = 1000, temp_dir: str = None, 
                                 num_processes: int = 4) -> Dict[str, str]:
    """
    Retrieve protein sequences from a large FASTA file (like UniRef90) efficiently.
    Uses multiple strategies for optimization:
    1. Index-based lookup with Biopython or SeqKit if available
    2. Parallel processing with multiple workers
    3. Optimized chunk size processing
    4. Deduplication of IDs to avoid retrieving the same sequence multiple times
    
    Args:
        protein_ids: List of protein IDs to retrieve.
        fasta_file: Path to the FASTA file.
        chunk_size: Number of IDs to process at once to avoid command line length limits.
        temp_dir: Directory for temporary files. If None, uses system temp directory.
        num_processes: Number of parallel processes to use for retrieval.
        
    Returns:
        Dictionary mapping protein IDs to their sequences.
    """
    import time
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor
    
    if not protein_ids:
        logging.warning("No protein IDs provided for sequence retrieval")
        return {}
        
    if not os.path.exists(fasta_file):
        logging.error(f"FASTA file not found: {fasta_file}")
        return {}
        
    start_time = time.time()
    sequences = {}
    
    # Deduplicate IDs (important when multiple hits per protein)
    unique_ids = list(set(protein_ids))
    total_ids = len(unique_ids)
    
    if len(unique_ids) < len(protein_ids):
        logging.info(f"Deduplicated input IDs from {len(protein_ids)} to {total_ids} unique IDs")
    
    # For very large datasets, consider sampling for debugging
    if total_ids > 1000000 and "DEBUG_SAMPLE" in os.environ:
        sample_size = int(os.environ.get("DEBUG_SAMPLE", 10000))
        logging.warning(f"DEBUG: Sampling {sample_size} IDs from {total_ids} for testing")
        import random
        random.seed(42)  # For reproducibility
        unique_ids = random.sample(unique_ids, min(sample_size, total_ids))
        total_ids = len(unique_ids)
    
    logging.info(f"Retrieving {total_ids} unique sequences from {fasta_file}")
    
    # Create temporary directory for IDs file if not provided
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp()
    else:
        os.makedirs(temp_dir, exist_ok=True)
        
    # First, try using the optimized UniRef90 retrieval function if it looks like a UniRef90 database
    # We determine this by checking if:
    # 1. The filename contains "uniref90" or "uniref_90" (case insensitive)
    # 2. Or if any of the IDs have a UniRef90_ prefix
    is_uniref90_db = False
    fasta_name_lower = fasta_file.lower()
    if "uniref90" in fasta_name_lower or "uniref_90" in fasta_name_lower:
        is_uniref90_db = True
        logging.info("Detected UniRef90 database from filename")
    else:
        # Check if any IDs have UniRef90_ prefix
        for pid in unique_ids[:min(100, len(unique_ids))]:  # Check first 100 IDs at most
            if pid.startswith('UniRef90_'):
                is_uniref90_db = True
                logging.info("Detected UniRef90 database from ID prefixes")
                break
    
    if is_uniref90_db:
        logging.info("Using optimized UniRef90 sequence retrieval")
        try:
            # Try the optimized UniRef90 retrieval function first
            sequences = retrieve_sequences_with_uniref90_prefix(
                unique_ids, 
                fasta_file, 
                temp_dir, 
                num_processes
            )
            
            if sequences:
                logging.info(f"Successfully retrieved {len(sequences)} sequences with optimized UniRef90 method")
                return sequences
            else:
                logging.warning("UniRef90 optimized retrieval returned no sequences, trying alternative methods")
        except Exception as e:
            logging.warning(f"UniRef90 optimized retrieval failed: {e}, trying alternative methods")
    
    # Check if we should use memory-mapped batch processing for extreme datasets
    use_mmap = total_ids > 100000 and os.path.getsize(fasta_file) > 10*1024*1024*1024  # > 10GB
    if use_mmap:
        logging.info("Using memory-mapped batch processing for large dataset")
        return _retrieve_with_mmap_batches(unique_ids, fasta_file, temp_dir, chunk_size, num_processes)
    
    # Check if we can use a faster method with SeqKit
    try:
        result = subprocess.run(["seqkit", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0:
            logging.info("Found seqkit - using fast index-based sequence retrieval")
            return _retrieve_with_seqkit(unique_ids, fasta_file, temp_dir, num_processes)
    except (FileNotFoundError, subprocess.SubprocessError):
        logging.info("SeqKit not found, trying other methods")
    
    # Expand protein IDs to include different variants that could be in the FASTA file
    expanded_ids = set()
    for pid in unique_ids:
        expanded_ids.add(pid)
        
        # Handle UniRef90 prefix
        if pid.startswith('UniRef90_'):
            base_id = pid.split('UniRef90_')[1]
            expanded_ids.add(base_id)
        else:
            expanded_ids.add(f"UniRef90_{pid}")
            
    logging.info(f"Expanded {len(unique_ids)} IDs to {len(expanded_ids)} search patterns")
    expanded_ids = list(expanded_ids)
    
    # Try to check if an index exists
    index_file = f"{fasta_file}.fai"
    use_index = False
    
    if os.path.exists(index_file):
        logging.info(f"Found FASTA index at {index_file}, using indexed retrieval")
        use_index = True
    else:
        # Try to create an index if it doesn't exist
        try:
            logging.info(f"Creating FASTA index for faster retrieval")
            from Bio import SeqIO
            SeqIO.index_db(f"{fasta_file}.idx", fasta_file, "fasta")
            use_index = True
            logging.info(f"Index created successfully")
        except Exception as e:
            logging.warning(f"Failed to create index: {e}. Using standard retrieval methods.")
    
    if use_index:
        try:
            from Bio import SeqIO
            record_dict = SeqIO.index(fasta_file, "fasta")
            
            for pid in expanded_ids:
                try:
                    record = record_dict.get(pid)
                    if record:
                        sequences[pid] = str(record.seq)
                        # Also store variants
                        if pid.startswith('UniRef90_'):
                            base_id = pid.split('UniRef90_')[1]
                            sequences[base_id] = str(record.seq)
                except Exception as e:
                    # Just continue if we can't find this ID
                    continue
                    
            if sequences:
                logging.info(f"Retrieved {len(sequences)} sequences using indexed lookup")
                return sequences
        except Exception as e:
            logging.warning(f"Index-based retrieval failed: {e}. Falling back to other methods.")
    
    # Try using seqtk which is faster and more efficient
    if _has_seqtk():
        sequences = _process_with_seqtk(expanded_ids, fasta_file, temp_dir, chunk_size)
        if sequences:
            logging.info(f"Successfully retrieved {len(sequences)} sequences using seqtk in {time.time() - start_time:.2f} seconds")
            return sequences
    
    # If we have multiple cores, use parallel processing for the fallback grep method
    if num_processes > 1 and total_ids > chunk_size:
        try:
            # Split the IDs into optimal chunks for parallel processing
            # For very large datasets, use larger chunks to reduce overhead
            if total_ids > 100000:
                optimal_chunk_size = max(chunk_size, total_ids // (num_processes * 2))
                logging.info(f"Adjusting chunk size to {optimal_chunk_size} for large dataset")
            else:
                optimal_chunk_size = chunk_size
                
            id_chunks = [expanded_ids[i:i+optimal_chunk_size] for i in range(0, len(expanded_ids), optimal_chunk_size)]
            logging.info(f"Using {num_processes} processes to retrieve sequences in {len(id_chunks)} chunks")
            
            # Define a function to process a single chunk
            def process_chunk(chunk_ids):
                return _grep_sequences(chunk_ids, fasta_file, temp_dir)
            
            # Process chunks in parallel with progress tracking
            completed = 0
            sequences_dict = {}
            with ProcessPoolExecutor(max_workers=num_processes) as executor:
                future_to_chunk = {executor.submit(process_chunk, chunk): i for i, chunk in enumerate(id_chunks)}
                
                for future in concurrent.futures.as_completed(future_to_chunk):
                    chunk_idx = future_to_chunk[future]
                    try:
                        result = future.result()
                        sequences_dict.update(result)
                        completed += 1
                        
                        # Log progress periodically
                        if completed % max(1, len(id_chunks) // 10) == 0:
                            logging.info(f"Progress: {completed}/{len(id_chunks)} chunks processed ({completed*100/len(id_chunks):.1f}%)")
                    except Exception as e:
                        logging.error(f"Error processing chunk {chunk_idx}: {e}")
                        
            logging.info(f"Retrieved {len(sequences_dict)} sequences using parallel grep in {time.time() - start_time:.2f} seconds")
            return sequences_dict
            
        except Exception as e:
            logging.warning(f"Parallel processing failed: {e}. Using sequential processing.")
    
    # Fallback to sequential grep method
    logging.info("Using sequential grep method for sequence retrieval")
    sequences = _grep_sequences(expanded_ids, fasta_file, temp_dir, chunk_size)
    
    logging.info(f"Retrieved {len(sequences)} sequences in {time.time() - start_time:.2f} seconds")
    return sequences

def _retrieve_with_mmap_batches(protein_ids, fasta_file, temp_dir, chunk_size=1000, num_processes=4):
    """Use memory mapping to process extremely large FASTA files in batches"""
    import mmap
    import os
    import time
    from concurrent.futures import ProcessPoolExecutor
    import concurrent.futures
    
    start_time = time.time()
    sequences = {}
    batch_size = chunk_size * 10  # Process larger batches
    
    # Create index file of sequence positions if it doesn't exist
    index_file = f"{fasta_file}.mmidx"
    if not os.path.exists(index_file):
        logging.info(f"Creating memory map index for {fasta_file}...")
        
        # Create a simple index - map header to file position
        headers_to_pos = {}
        with open(fasta_file, 'rb') as f:
            mmapped_file = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            pos = 0
            
            while True:
                header_pos = mmapped_file.find(b'>', pos)
                if header_pos == -1:
                    break
                    
                # Find end of header line
                eol_pos = mmapped_file.find(b'\n', header_pos)
                if eol_pos == -1:
                    break
                    
                # Extract header
                header = mmapped_file[header_pos+1:eol_pos].decode('ascii').split()[0]
                headers_to_pos[header] = header_pos
                
                # Move to next position
                pos = eol_pos + 1
                
                # Show progress for large files
                if len(headers_to_pos) % 1000000 == 0:
                    logging.info(f"Indexed {len(headers_to_pos)} sequences...")
                    
            mmapped_file.close()
            
        # Save index to file
        import pickle
        with open(index_file, 'wb') as f:
            pickle.dump(headers_to_pos, f)
            
        logging.info(f"Created index with {len(headers_to_pos)} entries in {time.time() - start_time:.2f} seconds")
    
    # Load index
    import pickle
    with open(index_file, 'rb') as f:
        headers_to_pos = pickle.load(f)
    
    logging.info(f"Loaded index with {len(headers_to_pos)} entries")
    
    # Function to extract sequences for a batch of IDs
    def extract_batch(id_batch):
        batch_sequences = {}
        
        # Find positions in file
        positions = []
        matched_ids = []
        
        for pid in id_batch:
            # Try different variants of ID
            if pid in headers_to_pos:
                positions.append(headers_to_pos[pid])
                matched_ids.append(pid)
            elif pid.startswith('UniRef90_') and pid[9:] in headers_to_pos:
                positions.append(headers_to_pos[pid[9:]])
                matched_ids.append(pid)
            elif f"UniRef90_{pid}" in headers_to_pos:
                positions.append(headers_to_pos[f"UniRef90_{pid}"])
                matched_ids.append(pid)
                
        if not positions:
            return {}
            
        # Extract sequences using memory mapping
        with open(fasta_file, 'rb') as f:
            mmapped_file = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            
            for i, pos in enumerate(positions):
                pid = matched_ids[i]
                
                # Find header end
                header_end = mmapped_file.find(b'\n', pos)
                if header_end == -1:
                    continue
                    
                # Find next header or end of file
                next_header = mmapped_file.find(b'>', header_end)
                if next_header == -1:
                    next_header = mmapped_file.size()
                    
                # Extract sequence lines
                seq_lines = mmapped_file[header_end+1:next_header].decode('ascii').replace('\n', '')
                batch_sequences[pid] = seq_lines
                
                # Also store without UniRef90_ prefix if present
                if pid.startswith('UniRef90_'):
                    batch_sequences[pid[9:]] = seq_lines
                    
            mmapped_file.close()
            
        return batch_sequences
    
    # Process in batches with multiple processes
    batches = [protein_ids[i:i+batch_size] for i in range(0, len(protein_ids), batch_size)]
    logging.info(f"Processing {len(batches)} batches with {num_processes} processes")
    
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        batch_futures = [executor.submit(extract_batch, batch) for batch in batches]
        
        completed = 0
        for future in concurrent.futures.as_completed(batch_futures):
            try:
                batch_result = future.result()
                sequences.update(batch_result)
                completed += 1
                
                if completed % max(1, len(batches) // 10) == 0:
                    logging.info(f"Progress: {completed}/{len(batches)} batches processed ({len(sequences)} sequences found)")
            except Exception as e:
                logging.error(f"Error processing batch: {e}")
    
    logging.info(f"Retrieved {len(sequences)} sequences using memory mapping in {time.time() - start_time:.2f} seconds")
    return sequences

def _has_seqtk():
    """Check if seqtk is available"""
    try:
        result = subprocess.run(["seqtk"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.returncode == 1  # seqtk without args returns 1 but prints help
    except (FileNotFoundError, subprocess.SubprocessError):
        return False

def _check_seqkit_version():
    """Check if seqkit is available and determine if it has faidx capability"""
    try:
        # Try to get version
        result = subprocess.run(["seqkit", "version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            version_text = result.stdout.strip()
            logging.info(f"Found seqkit: {version_text}")
            
            # Check if faidx command is available (this should be in most versions)
            faidx_check = subprocess.run(["seqkit", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            has_faidx = "faidx" in faidx_check.stdout or "faidx" in faidx_check.stderr
            
            return True, has_faidx
        else:
            # Try to check just if seqkit is available
            help_result = subprocess.run(["seqkit", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if help_result.returncode == 0 or help_result.returncode == 1:  # Some versions return 1 for help
                return True, False  # seqkit available but assume no faidx capability
        
        return False, False
    except (FileNotFoundError, subprocess.SubprocessError):
        return False, False

def _retrieve_with_seqkit(protein_ids, fasta_file, temp_dir, num_processes=4):
    """Use seqkit for fast indexed retrieval"""
    sequences = {}
    ids_file = os.path.join(temp_dir, "all_ids.txt")
    
    # Process IDs to ensure we can match both with and without UniRef90_ prefix
    processed_ids = set()
    for pid in protein_ids:
        processed_ids.add(pid)  # Original ID
        # Handle UniRef90 prefix
        if pid.startswith('UniRef90_'):
            processed_ids.add(pid[9:])  # Remove prefix
        else:
            processed_ids.add(f"UniRef90_{pid}")  # Add prefix
    
    # Write all processed IDs to a file
    with open(ids_file, 'w') as f:
        for pid in processed_ids:
            f.write(f"{pid}\n")
    
    try:
        # Check seqkit version and capabilities
        has_seqkit, has_faidx = _check_seqkit_version()
        if not has_seqkit:
            raise Exception("Seqkit not found")
            
        # Only try to build index if seqkit has faidx capability
        idx_file = f"{fasta_file}.fai"
        if has_faidx and not os.path.exists(idx_file):
            logging.info(f"Building FASTA index for {fasta_file}")
            try:
                subprocess.run(["seqkit", "faidx", fasta_file], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                logging.info(f"FASTA index built successfully")
            except subprocess.CalledProcessError as e:
                logging.warning(f"Failed to build FASTA index: {e}")
                # Continue anyway, seqkit grep can still work without index but will be slower
        elif not has_faidx:
            logging.info("This version of seqkit doesn't support indexing. Continuing without index.")
        else:
            logging.info(f"FASTA index already exists at {idx_file}")
        
        # Extract sequences using SeqKit with pattern match
        logging.info(f"Retrieving sequences for {len(processed_ids)} IDs using seqkit")
        
        # Build the command based on available options
        cmd = ["seqkit", "grep", "-f", ids_file, fasta_file]
        
        # Add other options if supported
        try:
            # Check if seqkit grep supports -n and -p
            grep_help = subprocess.run(["seqkit", "grep", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            help_text = grep_help.stdout + grep_help.stderr
            
            if "-n," in help_text or "--only-id" in help_text:
                cmd.append("-n")  # Only ID matching
            
            if "-p," in help_text or "--pattern" in help_text:
                cmd.append("-p")  # Pattern matching
                
            if "-j," in help_text or "--threads" in help_text:
                cmd.extend(["-j", str(num_processes)])  # Multi-threading
                
        except Exception:
            # If can't determine options, use basic command
            pass
            
        logging.info(f"Running seqkit command: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Parse the output FASTA
        from Bio import SeqIO
        for record in SeqIO.parse(io.StringIO(result.stdout.decode("utf-8")), "fasta"):
            record_id = record.id
            sequence = str(record.seq)
            
            # Store with original ID
            sequences[record_id] = sequence
            
            # Extract ID from UniRef90 header patterns
            # Handle various header formats:
            # 1. UniRef90_P12345 (standard UniRef)
            # 2. >UniRef90_P12345 other text
            # 3. >tr|P12345|... (UniProt TrEMBL)
            # 4. >sp|P12345|... (UniProt SwissProt)
            
            # Store the sequence under all relevant ID formats
            if record_id.startswith('UniRef90_'):
                base_id = record_id[9:]  # Remove UniRef90_ prefix
                sequences[base_id] = sequence
                
                # Also check for UniProt-style IDs (with |)
                if '|' in base_id:
                    parts = base_id.split('|')
                    for part in parts:
                        if part.strip():
                            sequences[part.strip()] = sequence
            else:
                # If it doesn't have UniRef90_ prefix, add it as an alternative
                sequences[f"UniRef90_{record_id}"] = sequence
                
                # Also check for UniProt-style IDs (with |)
                if '|' in record_id:
                    parts = record_id.split('|')
                    for part in parts:
                        if part.strip():
                            sequences[part.strip()] = sequence
                            sequences[f"UniRef90_{part.strip()}"] = sequence
        
        logging.info(f"Retrieved {len(sequences)} sequences using SeqKit")
        
    except Exception as e:
        logging.warning(f"SeqKit retrieval failed: {e}")
    
    finally:
        if os.path.exists(ids_file):
            os.remove(ids_file)
    
    return sequences

def retrieve_sequences_with_uniref90_prefix(protein_ids, fasta_file, temp_dir=None, num_processes=4):
    """
    Optimized function for retrieving sequences from UniRef90 database using seqkit.
    This function is designed to work directly with UniRef90 prefixes and automatically
    handle various ID formats.
    
    Args:
        protein_ids: List of protein IDs to retrieve (with or without UniRef90_ prefix)
        fasta_file: Path to the UniRef90 FASTA file
        temp_dir: Directory for temporary files. If None, uses system temp directory
        num_processes: Number of parallel processes to use
        
    Returns:
        Dictionary mapping protein IDs to their sequences
    """
    if not temp_dir:
        temp_dir = tempfile.mkdtemp()
    else:
        os.makedirs(temp_dir, exist_ok=True)
    
    logging.info(f"Retrieving {len(protein_ids)} sequences from UniRef90 database using seqkit")
    
    # Check if seqkit is available and get capabilities
    has_seqkit, has_faidx = _check_seqkit_version()
    if not has_seqkit:
        logging.warning("SeqKit not found. Please install seqkit for faster sequence retrieval.")
        return {}
    
    # Process IDs to ensure correct format for UniRef90 database
    uniref_ids = []
    id_mapping = {}  # Maps UniRef90 IDs to original IDs
    
    for pid in protein_ids:
        if pid.startswith('UniRef90_'):
            uniref_id = pid
        else:
            uniref_id = f"UniRef90_{pid}"
        
        uniref_ids.append(uniref_id)
        id_mapping[uniref_id] = pid
    
    # Write IDs to a temporary file
    ids_file = os.path.join(temp_dir, "uniref90_ids.txt")
    with open(ids_file, 'w') as f:
        for uid in uniref_ids:
            f.write(f"{uid}\n")
    
    sequences = {}
    
    try:
        # Try to build or verify faidx index
        idx_file = f"{fasta_file}.fai"
        if has_faidx and not os.path.exists(idx_file):
            logging.info(f"Building FASTA index for {fasta_file} using seqkit faidx")
            try:
                subprocess.run(["seqkit", "faidx", fasta_file], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                logging.info(f"FASTA index built successfully")
            except subprocess.CalledProcessError as e:
                logging.warning(f"Failed to build FASTA index: {e}. Continuing without index.")
        elif not has_faidx:
            logging.info("This version of seqkit doesn't support faidx. Continuing without index.")
        else:
            logging.info(f"FASTA index already exists at {idx_file}")
        
        # Extract sequences using pattern matching
        # Build the command based on available options
        cmd = ["seqkit", "grep", "-f", ids_file, fasta_file]
        
        # Check if seqkit grep supports additional options
        try:
            grep_help = subprocess.run(["seqkit", "grep", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            help_text = grep_help.stdout + grep_help.stderr
            
            if "-p," in help_text or "--pattern" in help_text:
                cmd.append("-p")  # Pattern matching
                
            if "-n," in help_text or "--only-id" in help_text:
                cmd.append("-n")  # ID matching
                
            if "-j," in help_text or "--threads" in help_text:
                cmd.extend(["-j", str(num_processes)])  # Multi-threading
        except Exception:
            # If can't determine options, use basic command
            pass
        
        logging.info(f"Running seqkit command: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Parse the output FASTA
        from Bio import SeqIO
        fasta_output = result.stdout.decode("utf-8")
        
        # Check if we got any output
        if not fasta_output.strip():
            logging.warning("SeqKit returned empty results. Check your ID format and database.")
            return {}
        
        for record in SeqIO.parse(io.StringIO(fasta_output), "fasta"):
            record_id = record.id
            sequence = str(record.seq)
            
            # Store with original record ID
            sequences[record_id] = sequence
            
            # Store with original requested ID (without UniRef90_ if that's how it was requested)
            if record_id in id_mapping:
                orig_id = id_mapping[record_id]
                sequences[orig_id] = sequence
            
            # Handle case where record_id might have additional information
            for uniref_id in uniref_ids:
                if record_id.startswith(uniref_id):
                    orig_id = id_mapping[uniref_id]
                    sequences[orig_id] = sequence
                    break
        
        logging.info(f"Retrieved {len(sequences)} sequences from UniRef90 database")
        
    except Exception as e:
        logging.error(f"Error retrieving sequences with seqkit: {e}")
    
    finally:
        # Clean up
        if os.path.exists(ids_file):
            os.remove(ids_file)
    
    return sequences

def _process_with_seqtk(protein_ids, fasta_file, temp_dir, chunk_size=1000):
    """Process sequence retrieval using seqtk in chunks"""
    sequences = {}
    total_ids = len(protein_ids)
    
    # Process in chunks to avoid command line length limits
    for i in range(0, total_ids, chunk_size):
        chunk = protein_ids[i:i+chunk_size]
        logging.info(f"Processing chunk {i//chunk_size + 1}/{(total_ids-1)//chunk_size + 1} ({len(chunk)} IDs)")
        
        # Create a temporary file with the IDs
        ids_file = os.path.join(temp_dir, f"ids_chunk_{i}.txt")
        with open(ids_file, 'w') as f:
            for pid in chunk:
                f.write(f"{pid}\n")
        
        # Use seqtk to extract sequences
        cmd = ["seqtk", "subseq", fasta_file, ids_file]
        try:
            result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # Parse the output FASTA
            fasta_data = result.stdout
            for record in SeqIO.parse(io.StringIO(fasta_data), "fasta"):
                # Store with the full ID
                sequences[record.id] = str(record.seq)
                
                # Also store with cleaned ID variants
                id_parts = record.id.split('|')
                for part in id_parts:
                    if part.strip():
                        sequences[part.strip()] = str(record.seq)
                
                # Store without UniRef90_ prefix if present
                if record.id.startswith('UniRef90_'):
                    base_id = record.id.split('UniRef90_')[1]
                    sequences[base_id] = str(record.seq)
            
            # Clean up temporary file
            os.remove(ids_file)
            
        except subprocess.CalledProcessError:
            logging.warning("seqtk failed, trying alternative method...")
            os.remove(ids_file)
            # Continue to alternative method
            break
    
    return sequences

def _grep_sequences(protein_ids, fasta_file, temp_dir, chunk_size=1000):
    """Extract sequences using grep, optionally in chunks"""
    sequences = {}
    total_ids = len(protein_ids)
    
    # Create patterns that can match different header formats for each ID
    # Process in chunks to avoid command line length limits
    for i in range(0, total_ids, chunk_size):
        chunk = protein_ids[i:i+chunk_size]
        logging.info(f"Processing chunk {i//chunk_size + 1}/{(total_ids-1)//chunk_size + 1} ({len(chunk)} IDs)")
        
        # Create patterns for each ID
        patterns = []
        for pid in chunk:
            # Match ID at the start of the header
            patterns.append(f"^>{pid}\\b")
            # Match ID in UniProt format (sp|ID|...)
            patterns.append(f"^>.*\\|{pid}\\|")
            # Match ID in TrEMBL format (tr|ID|...)
            patterns.append(f"^>.*\\|{pid}$")
            # Match ID with UniRef90 prefix
            if not pid.startswith('UniRef90_'):
                patterns.append(f"^>UniRef90_{pid}\\b")
        
        # Join patterns with OR
        pattern = '|'.join(patterns)
        
        # Create a temporary file for the output
        output_file = os.path.join(temp_dir, f"sequences_chunk_{i}.fasta")
        
        # Try using ripgrep if available (much faster than grep)
        ripgrep_cmd = f"rg -A 100 -m 1 --no-ignore-case '{pattern}' {fasta_file} > {output_file}"
        grep_cmd = f"grep -A 100 -m 1 -E '{pattern}' {fasta_file} > {output_file}"
        
        try:
            # Try ripgrep first
            result = subprocess.run("rg --version", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode == 0:
                cmd = ripgrep_cmd
            else:
                cmd = grep_cmd
                
            subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE)
            
            # Check if output file exists and has content
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                try:
                    # Parse the output in FASTA format
                    with open(output_file, 'r') as f:
                        lines = f.readlines()
                        
                    current_id = None
                    current_seq = []
                    
                    for line in lines:
                        line = line.strip()
                        if line.startswith('>'):
                            # Save previous record if there was one
                            if current_id and current_seq:
                                sequences[current_id] = ''.join(current_seq)
                                
                                # Store ID variants
                                id_parts = current_id.split('|')
                                for part in id_parts:
                                    if part.strip():
                                        sequences[part.strip()] = ''.join(current_seq)
                                
                                # Store without UniRef90_ prefix
                                if current_id.startswith('UniRef90_'):
                                    base_id = current_id.split('UniRef90_')[1]
                                    sequences[base_id] = ''.join(current_seq)
                            
                            # Start new record
                            header = line[1:].strip()  # Remove '>' and whitespace
                            current_id = header.split()[0]  # First word is the ID
                            current_seq = []
                        elif current_id and not line.startswith('--'):  # Skip grep separator
                            current_seq.append(line)
                    
                    # Save the last record
                    if current_id and current_seq:
                        sequences[current_id] = ''.join(current_seq)
                        
                        # Store ID variants
                        id_parts = current_id.split('|')
                        for part in id_parts:
                            if part.strip():
                                sequences[part.strip()] = ''.join(current_seq)
                        
                        # Store without UniRef90_ prefix
                        if current_id.startswith('UniRef90_'):
                            base_id = current_id.split('UniRef90_')[1]
                            sequences[base_id] = ''.join(current_seq)
                except Exception as e:
                    logging.error(f"Error parsing FASTA file {output_file}: {e}")
            
            # Clean up the temporary file
            if os.path.exists(output_file):
                os.remove(output_file)
                
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to run grep command: {e}")
            if os.path.exists(output_file):
                os.remove(output_file)
    
    return sequences 