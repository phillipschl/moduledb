#!/usr/bin/env python3
"""
Validation module for NRPS C-A-T modules.
Contains functions for validating domain order and continuity.
"""

import logging
import os
import tempfile
from typing import List, Tuple, Dict, Any, Optional, Union
from scripts.alignment import run_hmmscan

def validate_domains_with_pfam(
    sequences: Dict[str, str], 
    pfam_hmm_path: str, 
    e_value_threshold: float = 1e-10,
    cpu: int = 4,
    temp_dir: Optional[str] = None
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Validate protein sequences against Pfam-A.hmm database to identify domains.
    
    Args:
        sequences: Dictionary mapping sequence IDs to their sequences.
        pfam_hmm_path: Path to the Pfam-A.hmm database file.
        e_value_threshold: E-value threshold for HMMER hits.
        cpu: Number of CPU threads to use.
        temp_dir: Directory for temporary files (created if None).
        
    Returns:
        Dictionary mapping sequence IDs to lists of domain hits.
    """
    # Check if Pfam-A.hmm exists
    if not os.path.exists(pfam_hmm_path):
        logging.error(f"Pfam-A.hmm file not found at {pfam_hmm_path}")
        return {}
    
    # Create a temp directory if not provided
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp()
    else:
        os.makedirs(temp_dir, exist_ok=True)
    
    # Create a temporary FASTA file with the sequences
    temp_fasta = os.path.join(temp_dir, "temp_sequences.fasta")
    
    try:
        # Write sequences to temporary FASTA file
        with open(temp_fasta, 'w') as f:
            for seq_id, sequence in sequences.items():
                f.write(f">{seq_id}\n{sequence}\n")
        
        # Run hmmscan against Pfam-A.hmm
        output_prefix = os.path.join(temp_dir, "pfam_scan")
        domtblout_file = run_hmmscan(pfam_hmm_path, temp_fasta, output_prefix, e_value_threshold, cpu)
        
        if not domtblout_file or not os.path.exists(domtblout_file):
            logging.error("hmmscan against Pfam-A.hmm failed")
            return {}
        
        # Parse the domtblout file
        hits_by_protein = parse_hmmscan_domtblout(domtblout_file, e_value_threshold)
        
        # Process domain hits to identify specific domains of interest
        processed_hits = {}
        
        # Define domain mappings for NRPS domains (Pfam domains known to be associated with NRPS)
        domain_mappings = {
            "Condensation": "C",
            "AMP-binding": "A",
            "PP-binding": "T",
            "Thioesterase": "TE",
            "Epimerase": "E",
            "Methyltransf_11": "MT", 
            "NMT": "MT",  # N-methyltransferase
            "Thioesterase": "TE",
            "Formyl_trans_N": "F",
            "Formyl_trans_C": "F",
            "Reductase_2": "R",
            "Reductase_3": "R",
            "PKS_AT": "AT"
        }
        
        # Process hits to identify domain types
        for protein_id, hits in hits_by_protein.items():
            processed_hits[protein_id] = []
            
            for hit in hits:
                domain_id = hit['domain_id']
                domain_type = None
                
                # Try direct mapping first
                for pfam_name, domain_code in domain_mappings.items():
                    if pfam_name.lower() in domain_id.lower():
                        domain_type = domain_code
                        break
                
                # Add to processed hits with identified domain type if found
                if domain_type:
                    processed_hits[protein_id].append({
                        'domain_type': domain_type,
                        'domain_id': domain_id,
                        'start': hit['start'],
                        'end': hit['end'],
                        'e_value': hit['e_value']
                    })
        
        logging.info(f"Validated {len(processed_hits)} sequences against Pfam-A.hmm")
        return processed_hits
        
    except Exception as e:
        logging.error(f"Error validating sequences against Pfam-A.hmm: {e}")
        return {}
    
    finally:
        # Clean up temporary files
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)

def validate_cat_order(domain_hits: List[Dict[str, Any]], max_gap_allowed: int) -> Tuple[bool, Optional[Tuple[int, int]]]:
    """
    Validate that domain hits form a valid C-A-T module with proper order and spacing.
    
    Args:
        domain_hits: List of domain hits, each a dictionary with keys:
                    'domain_type', 'start', 'end'.
        max_gap_allowed: Maximum allowed gap between consecutive domains.
        
    Returns:
        Tuple of (is_valid, coordinates) where:
        - is_valid: Boolean indicating if the module is valid
        - coordinates: Tuple of (start, end) if valid, None otherwise
    """
    if not domain_hits:
        return False, None
    
    # For combined C-A-T module hits, we just need to check if there's at least one hit
    # and return its coordinates
    if len(domain_hits) == 1 and domain_hits[0]['domain_type'] == 'CAT':
        return True, (domain_hits[0]['start'], domain_hits[0]['end'])
    
    # For individual domain hits, check the traditional C-A-T order
    # Sort domains by start position
    sorted_domains = sorted(domain_hits, key=lambda x: x['start'])
    
    # Check if we have at least one of each domain type
    domain_types = [d['domain_type'] for d in sorted_domains]
    if not ('C' in domain_types and 'A' in domain_types and 'T' in domain_types):
        return False, None
    
    # Find the first C domain
    c_domain = next((d for d in sorted_domains if d['domain_type'] == 'C'), None)
    if not c_domain:
        return False, None
    
    # Find the first A domain after C
    a_domain = next((d for d in sorted_domains if d['domain_type'] == 'A' and d['start'] > c_domain['start']), None)
    if not a_domain:
        return False, None
    
    # Find the first T domain after A
    t_domain = next((d for d in sorted_domains if d['domain_type'] == 'T' and d['start'] > a_domain['start']), None)
    if not t_domain:
        return False, None
    
    # Check gaps between domains
    if (a_domain['start'] - c_domain['end'] > max_gap_allowed or
        t_domain['start'] - a_domain['end'] > max_gap_allowed):
        return False, None
    
    return True, (c_domain['start'], t_domain['end'])

def parse_hmmscan_domtblout(domtblout_file: str, e_value_threshold: float) -> Dict[str, List[Dict[str, Any]]]:
    """
    Parse HMMER domtblout file and filter hits by E-value.
    
    Args:
        domtblout_file: Path to the HMMER domtblout file.
        e_value_threshold: E-value threshold for filtering hits.
        
    Returns:
        Dictionary mapping protein IDs to lists of domain hits.
        Each domain hit is a dictionary with keys:
        'domain_type', 'domain_id', 'start', 'end', 'e_value'.
    """
    hits_by_protein = {}
    
    # Ensure e_value_threshold is a float
    try:
        e_value_threshold = float(e_value_threshold)
    except (ValueError, TypeError):
        logging.error(f"Invalid e_value_threshold: {e_value_threshold}. Using default 1e-10.")
        e_value_threshold = 1e-10
    
    with open(domtblout_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split()
            
            if len(fields) < 22:
                logging.warning(f"Invalid line format in {domtblout_file}: {line.strip()}")
                continue
            
            raw_target_name = fields[0]  # Protein/sequence ID with potential UniRef90_ prefix
            
            # Extract just the ID part from UniRef90_<ID> format
            if raw_target_name.startswith('UniRef90_'):
                target_name = raw_target_name.split('UniRef90_')[1]
                logging.debug(f"Extracted ID {target_name} from {raw_target_name}")
            else:
                target_name = raw_target_name
            
            domain_id = fields[3]    # Domain ID (HMM name)
            i_evalue = float(fields[11])  # Domain i-evalue
            ali_from = int(fields[17])  # Alignment start in target sequence
            ali_to = int(fields[18])    # Alignment end in target sequence
            
            # Skip hits that don't meet the E-value threshold
            if i_evalue > e_value_threshold:
                continue
            
            # For custom HMM hits, we'll treat the entire hit as a C-A-T module
            # and let the validation function determine if it's valid
            domain_hit = {
                'domain_type': 'CAT',  # Combined C-A-T module
                'domain_id': domain_id,
                'start': ali_from,
                'end': ali_to,
                'e_value': i_evalue
            }
            
            if target_name not in hits_by_protein:
                hits_by_protein[target_name] = []
                
            hits_by_protein[target_name].append(domain_hit)
    
    # Sort domain hits by start position for each protein
    for protein_id, hits in hits_by_protein.items():
        hits_by_protein[protein_id] = sorted(hits, key=lambda x: x['start'])
    
    return hits_by_protein

def extract_cat_module_sequences(
    hits_by_protein: Dict[str, List[Dict[str, Any]]], 
    sequences_by_id: Dict[str, str],
    max_gap_allowed: int
) -> Dict[str, str]:
    """
    Extract C-A-T module sequences from proteins with valid modules.
    
    Args:
        hits_by_protein: Dictionary mapping protein IDs to lists of domain hits.
        sequences_by_id: Dictionary mapping protein IDs to their full amino acid sequences.
        max_gap_allowed: Maximum allowable amino acid gap between consecutive domains.
        
    Returns:
        Dictionary mapping module IDs to their amino acid sequences.
    """
    module_sequences = {}
    
    for protein_id, domain_hits in hits_by_protein.items():
        # Process protein_id to handle various ID formats
        possible_ids = [protein_id]
        
        # Handle UniRef90 prefix
        if protein_id.startswith('UniRef90_'):
            base_id = protein_id.split('UniRef90_')[1]
            possible_ids.extend([base_id, f"tr|{base_id}", f"sp|{base_id}"])
            logging.debug(f"Looking up {base_id} (extracted from {protein_id})")
        else:
            # Add possible variants of the ID
            possible_ids.extend([
                f"UniRef90_{protein_id}",
                f"tr|{protein_id}",
                f"sp|{protein_id}"
            ])
        
        # Try all possible ID formats
        full_sequence = None
        found_id = None
        
        for pid in possible_ids:
            if pid in sequences_by_id:
                full_sequence = sequences_by_id[pid]
                found_id = pid
                break
                
            # Also try checking if ID is part of a longer ID key
            for seq_id in sequences_by_id:
                if pid in seq_id.split('|') or pid in seq_id.split():
                    full_sequence = sequences_by_id[seq_id]
                    found_id = seq_id
                    break
            
            if full_sequence:
                break
        
        if not full_sequence:
            # If we've tried all variants and still couldn't find it, log a warning
            logging.warning(f"Sequence not found for protein {protein_id} (tried {', '.join(possible_ids)})")
            continue
        else:
            logging.debug(f"Found sequence for {protein_id} using ID {found_id}")
            
        # Validate C-A-T order and get module coordinates
        is_valid, coords = validate_cat_order(domain_hits, max_gap_allowed)
        
        if is_valid and coords:
            c_start, t_end = coords
            
            # Make sure the coordinates are within the sequence length
            if c_start < 1:
                c_start = 1
            if t_end > len(full_sequence):
                t_end = len(full_sequence)
                
            # Extract the module sequence
            # Note: Adjusting by -1 for 0-based indexing in Python strings
            try:
                module_seq = full_sequence[c_start-1:t_end]
                
                # Create a unique module ID
                module_id = f"{protein_id}_CAT_Module_{c_start}_{t_end}"
                
                module_sequences[module_id] = module_seq
                logging.info(f"Extracted C-A-T module: {module_id}, length: {len(module_seq)} aa")
            except Exception as e:
                logging.error(f"Error extracting module sequence for {protein_id}: {e}")
                logging.error(f"Coordinates: {c_start}-{t_end}, sequence length: {len(full_sequence)}")
                continue
    
    return module_sequences 