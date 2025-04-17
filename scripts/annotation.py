#!/usr/bin/env python3
"""
Annotation module for annotating C-A-T modules with biological information.
"""

import os
import time
import logging
import requests
import csv
import json
from typing import Dict, List, Any, Optional
from Bio import SeqIO

def fetch_uniprot_annotations(
    uniprot_id: str, 
    retry_count: int = 3, 
    retry_delay: int = 1
) -> Optional[Dict[str, Any]]:
    """
    Fetch annotations for a protein from UniProt API.
    
    Args:
        uniprot_id: UniProt ID of the protein.
        retry_count: Number of times to retry the request if it fails.
        retry_delay: Delay in seconds between retries.
        
    Returns:
        Dictionary of annotations, or None if the request fails.
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.json"
    
    for attempt in range(retry_count):
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise exception for HTTP errors
            
            data = response.json()
            
            # Extract relevant annotations
            annotations = {
                'uniprot_id': uniprot_id,
                'protein_name': data.get('protein', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                'organism': data.get('organism', {}).get('scientificName', ''),
                'gene_name': ', '.join([name.get('value', '') for name in data.get('gene', [])]) if 'gene' in data else '',
                'taxonomy': '; '.join([taxon.get('scientificName', '') for taxon in data.get('organism', {}).get('lineage', [])]),
                'sequence_length': data.get('sequence', {}).get('length', 0),
                'review_status': data.get('entryType', '')
            }
            
            logging.info(f"Successfully fetched annotations for {uniprot_id}")
            return annotations
            
        except requests.exceptions.RequestException as e:
            logging.warning(f"Failed to fetch annotations for {uniprot_id}, attempt {attempt+1}/{retry_count}: {e}")
            if attempt < retry_count - 1:
                time.sleep(retry_delay)
    
    logging.error(f"Failed to fetch annotations for {uniprot_id} after {retry_count} attempts")
    return None

def extract_uniprot_id_from_module_id(module_id: str) -> Optional[str]:
    """
    Extract the UniProt ID from a module ID.
    
    Args:
        module_id: Module ID in the format 'UniProtID_CAT_Module_Cstart_Tend'.
        
    Returns:
        UniProt ID or None if it can't be extracted.
    """
    try:
        # The module ID format is expected to be 'UniProtID_CAT_Module_Cstart_Tend'
        parts = module_id.split('_')
        if len(parts) < 2:
            return None
            
        # The UniProt ID is the first part
        uniprot_id = parts[0]
        return uniprot_id
        
    except Exception as e:
        logging.error(f"Failed to extract UniProt ID from module ID {module_id}: {e}")
        return None

def extract_coordinates_from_module_id(module_id: str) -> Optional[Dict[str, int]]:
    """
    Extract the coordinates from a module ID.
    
    Args:
        module_id: Module ID in the format 'UniProtID_CAT_Module_Cstart_Tend'.
        
    Returns:
        Dictionary with coordinates or None if it can't be extracted.
    """
    try:
        # The module ID format is expected to be 'UniProtID_CAT_Module_Cstart_Tend'
        parts = module_id.split('_')
        if len(parts) < 5:
            return None
            
        # The coordinates are the last two parts
        c_start = int(parts[-2])
        t_end = int(parts[-1])
        
        coordinates = {
            'C_start': c_start,
            'T_end': t_end
        }
        
        return coordinates
        
    except Exception as e:
        logging.error(f"Failed to extract coordinates from module ID {module_id}: {e}")
        return None

def annotate_modules(
    module_sequences: Dict[str, str],
    domain_hits_by_protein: Dict[str, List[Dict[str, Any]]]
) -> Dict[str, Dict[str, Any]]:
    """
    Annotate C-A-T modules with domain coordinates and source protein information.
    
    Args:
        module_sequences: Dictionary mapping module IDs to sequences.
        domain_hits_by_protein: Dictionary mapping source protein IDs to domain hits.
        
    Returns:
        Dictionary mapping module IDs to annotation dictionaries.
    """
    annotations = {}
    
    for module_id, sequence in module_sequences.items():
        # Extract the source protein ID and module coordinates
        uniprot_id = extract_uniprot_id_from_module_id(module_id)
        coordinates = extract_coordinates_from_module_id(module_id)
        
        if not uniprot_id or not coordinates:
            logging.warning(f"Could not parse module ID {module_id}, skipping annotation")
            continue
        
        # Fetch UniProt annotations
        uniprot_annotations = fetch_uniprot_annotations(uniprot_id)
        
        if not uniprot_annotations:
            logging.warning(f"Failed to fetch annotations for UniProt ID {uniprot_id}, using minimal annotations")
            uniprot_annotations = {'uniprot_id': uniprot_id}
        
        # Find domain boundaries within the module
        module_annotation = {
            'module_id': module_id,
            'source_uniprot_id': uniprot_id,
            'module_start': coordinates['C_start'],
            'module_end': coordinates['T_end'],
            'module_length': len(sequence),
            **uniprot_annotations
        }
        
        # If we have domain hits for the source protein, add domain coordinates
        if uniprot_id in domain_hits_by_protein:
            domain_hits = domain_hits_by_protein[uniprot_id]
            
            # Find C, A, T domains within the module coordinates
            c_domains = [hit for hit in domain_hits if hit['domain_type'] == 'C' 
                        and coordinates['C_start'] <= hit['start'] <= coordinates['T_end']]
            a_domains = [hit for hit in domain_hits if hit['domain_type'] == 'A' 
                        and coordinates['C_start'] <= hit['start'] <= coordinates['T_end']]
            t_domains = [hit for hit in domain_hits if hit['domain_type'] == 'T' 
                        and coordinates['C_start'] <= hit['start'] <= coordinates['T_end']]
            
            # Use the first domain of each type that falls within the module
            if c_domains:
                module_annotation['C_domain_start'] = c_domains[0]['start']
                module_annotation['C_domain_end'] = c_domains[0]['end']
            
            if a_domains:
                module_annotation['A_domain_start'] = a_domains[0]['start']
                module_annotation['A_domain_end'] = a_domains[0]['end']
            
            if t_domains:
                module_annotation['T_domain_start'] = t_domains[0]['start']
                module_annotation['T_domain_end'] = t_domains[0]['end']
        
        annotations[module_id] = module_annotation
    
    return annotations

def save_annotations_to_tsv(
    annotations: Dict[str, Dict[str, Any]], 
    output_file: str
) -> bool:
    """
    Save annotations to a TSV file.
    
    Args:
        annotations: Dictionary mapping module IDs to annotation dictionaries.
        output_file: Path to the output TSV file.
        
    Returns:
        True if the annotations were saved successfully, False otherwise.
    """
    if not annotations:
        logging.error("No annotations to save")
        return False
    
    try:
        # Get all unique keys from all annotation dictionaries
        fieldnames = set()
        for annotation in annotations.values():
            fieldnames.update(annotation.keys())
        
        # Sort the fieldnames for consistent output
        fieldnames = sorted(list(fieldnames))
        
        with open(output_file, 'w', newline='') as tsvfile:
            writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            for module_id, annotation in annotations.items():
                writer.writerow(annotation)
        
        logging.info(f"Successfully saved {len(annotations)} annotations to {output_file}")
        return True
        
    except Exception as e:
        logging.error(f"Failed to save annotations to {output_file}: {e}")
        return False

def save_fasta_with_annotations(
    sequences: Dict[str, str],
    annotations: Dict[str, Dict[str, Any]],
    output_file: str
) -> bool:
    """
    Save sequences to a FASTA file with annotations in the headers.
    
    Args:
        sequences: Dictionary mapping sequence IDs to their sequences.
        annotations: Dictionary mapping sequence IDs to annotation dictionaries.
        output_file: Path to the output FASTA file.
        
    Returns:
        True if the sequences were saved successfully, False otherwise.
    """
    try:
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        
        records = []
        
        for seq_id, sequence in sequences.items():
            # Create description with annotations if available
            description = ""
            if seq_id in annotations:
                annotation_items = []
                for key, value in annotations[seq_id].items():
                    if key != 'module_id':  # Skip the ID since it's already in the ID field
                        annotation_items.append(f"{key}={value}")
                description = " ".join(annotation_items)
            
            record = SeqRecord(
                Seq(sequence),
                id=seq_id,
                description=description
            )
            records.append(record)
        
        SeqIO.write(records, output_file, "fasta")
        logging.info(f"Successfully saved {len(records)} annotated sequences to {output_file}")
        return True
        
    except Exception as e:
        logging.error(f"Failed to save annotated sequences to {output_file}: {e}")
        return False 