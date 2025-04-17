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
                                 chunk_size: int = 500, temp_dir: str = None) -> Dict[str, str]:
    """
    Retrieve protein sequences from a large FASTA file (like UniRef90) efficiently.
    Uses blastdbcmd or seqtk if available, otherwise falls back to a chunked approach with grep.
    
    Args:
        protein_ids: List of protein IDs to retrieve.
        fasta_file: Path to the FASTA file.
        chunk_size: Number of IDs to process at once to avoid command line length limits.
        temp_dir: Directory for temporary files. If None, uses system temp directory.
        
    Returns:
        Dictionary mapping protein IDs to their sequences.
    """
    if not protein_ids:
        logging.warning("No protein IDs provided for sequence retrieval")
        return {}
        
    if not os.path.exists(fasta_file):
        logging.error(f"FASTA file not found: {fasta_file}")
        return {}
        
    sequences = {}
    total_ids = len(protein_ids)
    logging.info(f"Retrieving {total_ids} sequences from {fasta_file}")
    
    # Create temporary directory for IDs file if not provided
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp()
    else:
        os.makedirs(temp_dir, exist_ok=True)
    
    # Expand protein IDs to include different variants that could be in the FASTA file
    expanded_ids = set()
    for pid in protein_ids:
        expanded_ids.add(pid)
        
        # Handle UniRef90 prefix
        if pid.startswith('UniRef90_'):
            base_id = pid.split('UniRef90_')[1]
            expanded_ids.add(base_id)
        else:
            expanded_ids.add(f"UniRef90_{pid}")
            
    logging.info(f"Expanded {len(protein_ids)} IDs to {len(expanded_ids)} search patterns")
    protein_ids = list(expanded_ids)
    
    # Try using seqtk which is faster and more efficient
    try:
        # Process in chunks to avoid command line length limits
        for i in range(0, len(protein_ids), chunk_size):
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
                    
                    # Store first word of description (usually the ID)
                    first_word = record.id.split()[0]
                    sequences[first_word] = str(record.seq)
                
                # Clean up temporary file
                os.remove(ids_file)
                
            except subprocess.CalledProcessError:
                logging.warning("seqtk failed, trying alternative method...")
                os.remove(ids_file)
                # Continue to alternative method
                break
        
        if sequences:
            logging.info(f"Successfully retrieved {len(sequences)} sequences using seqtk")
            return sequences
    
    except (FileNotFoundError, subprocess.SubprocessError):
        logging.warning("seqtk not available, using fallback method")
    
    # Fallback method: Create expanded patterns that match different header formats
    try:
        # Create patterns that can match different header formats for each ID
        # Process in chunks to avoid command line length limits
        for i in range(0, len(protein_ids), chunk_size):
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
            
            # Get header lines and then some lines after to get the sequence
            # We'll use -A 30 to get a reasonable chunk of sequence, which should be enough for most headers
            cmd = f"grep -A 30 -E '{pattern}' {fasta_file} > {output_file}"
            try:
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
                        
                        # Try more direct parsing since BioPython might fail on partial records
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
            except subprocess.CalledProcessError as e:
                logging.error(f"Error running grep command: {e}")
                logging.error(f"Stderr: {e.stderr}")
                
            # Clean up temporary file
            if os.path.exists(output_file):
                os.remove(output_file)
    
    except Exception as e:
        logging.error(f"Error retrieving sequences from FASTA: {e}")
    
    # Log some diagnostics
    if sequences:
        logging.info(f"Successfully retrieved {len(sequences)} sequences")
        logging.debug(f"Retrieved sequence IDs: {list(sequences.keys())[:10]}...")
    else:
        logging.warning("No sequences retrieved from FASTA file")
        # Try a simple test to verify file access
        try:
            with open(fasta_file, 'r') as f:
                first_line = f.readline().strip()
            logging.info(f"First line of FASTA file: {first_line}")
        except Exception as e:
            logging.error(f"Error reading FASTA file: {e}")
    
    return sequences 