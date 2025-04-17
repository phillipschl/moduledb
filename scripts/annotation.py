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
    Fetch annotations for a protein from UniProt API or UniRef90 API.
    
    Args:
        uniprot_id: UniProt ID or UniRef90 ID of the protein.
        retry_count: Number of times to retry the request if it fails.
        retry_delay: Delay in seconds between retries.
        
    Returns:
        Dictionary of annotations, or None if the request fails.
    """
    # Check if the ID is a UniRef90 ID
    is_uniref90 = uniprot_id.startswith('UniRef90_')
    
    # Extract the base UniProt ID if it's a UniRef90 ID
    if is_uniref90:
        base_id = uniprot_id[9:]  # Remove 'UniRef90_' prefix
        # Handle case where base ID might be a UniProt entry with pipes
        if '|' in base_id:
            parts = base_id.split('|')
            base_id = parts[1] if len(parts) > 1 else base_id
    else:
        base_id = uniprot_id
    
    # Try the UniRef90 API first for UniRef90 IDs
    if is_uniref90:
        url = f"https://rest.uniprot.org/uniref/UniRef90_{base_id}"
    else:
        # For standard UniProt IDs, use the UniProt API
        url = f"https://rest.uniprot.org/uniprotkb/{base_id}"
    
    for attempt in range(retry_count):
        try:
            headers = {"Accept": "application/json"}
            response = requests.get(url, headers=headers)
            response.raise_for_status()  # Raise exception for HTTP errors
            
            data = response.json()
            
            # Process differently based on API type
            if is_uniref90:
                # Process UniRef90 API response
                representative = data.get('representativeMember', {})
                member_data = representative.get('proteinAttributes', {})
                organism_data = representative.get('organism', {})
                
                annotations = {
                    'uniprot_id': uniprot_id,
                    #'cluster_id': data.get('id', ''),
                    #'cluster_name': data.get('name', ''),
                    'protein_name': member_data.get('submittedName', '') or member_data.get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    'organism': organism_data.get('scientificName', ''),
                    'taxonomy': organism_data.get('taxonId', ''),
                    'sequence_length': data.get('memberCount', 0),
                    'identity': data.get('identity', ''),
                    #'cluster_size': data.get('memberCount', 0),
                    'seed_id': base_id
                }
            else:
                # Process UniProtKB API response
                annotations = {
                    'uniprot_id': uniprot_id,
                    'protein_name': data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    'organism': data.get('organism', {}).get('scientificName', ''),
                    'gene_name': ', '.join([name.get('value', '') for name in data.get('genes', [])]) if 'genes' in data else '',
                    'taxonomy': '; '.join([taxon.get('scientificName', '') for taxon in data.get('organismHosts', [])]),
                    'sequence_length': data.get('sequence', {}).get('length', 0),
                    'review_status': data.get('entryType', '')
                }
            
            logging.info(f"Successfully fetched annotations for {uniprot_id} using {'UniRef90' if is_uniref90 else 'UniProt'} API")
            return annotations
            
        except requests.exceptions.RequestException as e:
            logging.warning(f"Failed to fetch annotations for {uniprot_id}, attempt {attempt+1}/{retry_count}: {e}")
            
            if attempt < retry_count - 1:
                time.sleep(retry_delay)
                
                # If UniRef90 API failed, try UniProt API as fallback
                if is_uniref90 and attempt == 0:
                    logging.info(f"Trying standard UniProt API as fallback for {uniprot_id}")
                    url = f"https://rest.uniprot.org/uniprotkb/{base_id}"
                    is_uniref90 = False  # Switch to standard UniProt processing
    
    logging.error(f"Failed to fetch annotations for {uniprot_id} after {retry_count} attempts")
    return None

def extract_uniprot_id_from_module_id(module_id: str) -> Optional[str]:
    """
    Extract the UniProt ID or UniRef90 ID from a module ID.
    
    Args:
        module_id: Module ID in the format 'UniProtID_CAT_Module_Cstart_Tend'
                  or 'UniRef90_UniProtID_CAT_Module_Cstart_Tend'.
        
    Returns:
        UniProt ID or UniRef90 ID, or None if it can't be extracted.
    """
    try:
        # Split the module ID on underscores
        parts = module_id.split('_')
        if len(parts) < 2:
            return None
        
        # Check if first part is "UniRef90"
        if parts[0] == "UniRef90":
            # This is a UniRef90 ID, which consists of 'UniRef90' and the second part
            if len(parts) < 3:
                return None
            uniprot_id = f"{parts[0]}_{parts[1]}"
            
            # Handle case where UniRef90 might have underscores in its ID
            # Look for a subpattern like "CAT_Module" which should follow the UniRef90 ID
            for i in range(2, len(parts)):
                if parts[i] == "CAT" and i+1 < len(parts) and parts[i+1] == "Module":
                    break
                uniprot_id += f"_{parts[i]}"
                
            return uniprot_id
        else:
            # This is a standard UniProt ID (the first part)
            # Check if it contains pipe (|) characters which need to be handled
            if '|' in parts[0]:
                # For UniProt format like sp|P12345|NAME, extract the accession (P12345)
                pipe_parts = parts[0].split('|')
                if len(pipe_parts) >= 2:
                    return pipe_parts[1]
            
            return parts[0]
            
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
            logging.warning(f"Failed to fetch annotations for ID {uniprot_id}, using minimal annotations")
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
        # Check several possible ID formats in the domain_hits_by_protein dictionary
        domain_hits = None
        
        # Check for the ID directly
        if uniprot_id in domain_hits_by_protein:
            domain_hits = domain_hits_by_protein[uniprot_id]
        
        # If it's a UniRef90 ID, try the base ID
        elif uniprot_id.startswith('UniRef90_'):
            base_id = uniprot_id[9:]  # Remove 'UniRef90_' prefix
            if base_id in domain_hits_by_protein:
                domain_hits = domain_hits_by_protein[base_id]
            elif '|' in base_id:
                # Handle format like UniRef90_sp|P12345|NAME
                parts = base_id.split('|')
                if len(parts) > 1 and parts[1] in domain_hits_by_protein:
                    domain_hits = domain_hits_by_protein[parts[1]]
        
        # If it's a base ID, try the UniRef90 version
        elif f"UniRef90_{uniprot_id}" in domain_hits_by_protein:
            domain_hits = domain_hits_by_protein[f"UniRef90_{uniprot_id}"]
            
        # Process domain hits if found
        if domain_hits:
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

def fetch_uniref90_cluster_info(
    cluster_id: str,
    retry_count: int = 3,
    retry_delay: int = 1
) -> Optional[Dict[str, Any]]:
    """
    Fetch information about a UniRef90 cluster.
    
    Args:
        cluster_id: UniRef90 cluster ID (with or without UniRef90_ prefix).
        retry_count: Number of times to retry the request if it fails.
        retry_delay: Delay in seconds between retries.
        
    Returns:
        Dictionary of cluster information, or None if the request fails.
    """
    # Add UniRef90_ prefix if not present
    if not cluster_id.startswith('UniRef90_'):
        cluster_id = f"UniRef90_{cluster_id}"
    
    url = f"https://rest.uniprot.org/uniref/{cluster_id}"
    headers = {"Accept": "application/json"}
    
    for attempt in range(retry_count):
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            
            data = response.json()
            
            # Extract cluster members
            members = []
            if 'members' in data:
                for member in data.get('members', []):
                    member_info = {
                        'id': member.get('memberId', ''),
                        'name': member.get('memberIdType', ''),
                        'organism': member.get('organism', {}).get('scientificName', ''),
                        'taxon_id': member.get('organism', {}).get('taxonId', '')
                    }
                    members.append(member_info)
            
            # Get representative member information
            representative = data.get('representativeMember', {})
            seed_id = representative.get('memberId', '')
            
            # Build the cluster info
            cluster_info = {
                'cluster_id': data.get('id', ''),
                'name': data.get('name', ''),
                'updated': data.get('updated', ''),
                'member_count': data.get('memberCount', 0),
                'identity': data.get('identity', ''),
                'seed_id': seed_id,
                'members': members,
                'organism_groups': [],
                'taxonomy_distribution': {}
            }
            
            # Calculate taxonomy distribution
            taxonomy_counts = {}
            for member in members:
                organism = member.get('organism', '')
                if organism:
                    taxonomy_counts[organism] = taxonomy_counts.get(organism, 0) + 1
            
            # Sort and add to the cluster info
            taxonomy_distribution = [{'organism': org, 'count': count} 
                                   for org, count in sorted(taxonomy_counts.items(), 
                                                          key=lambda x: x[1], reverse=True)]
            cluster_info['taxonomy_distribution'] = taxonomy_distribution
            
            logging.info(f"Successfully fetched UniRef90 cluster info for {cluster_id}")
            return cluster_info
            
        except requests.exceptions.RequestException as e:
            logging.warning(f"Failed to fetch UniRef90 cluster info for {cluster_id}, attempt {attempt+1}/{retry_count}: {e}")
            if attempt < retry_count - 1:
                time.sleep(retry_delay)
    
    logging.error(f"Failed to fetch UniRef90 cluster info for {cluster_id} after {retry_count} attempts")
    return None 