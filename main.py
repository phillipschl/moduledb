#!/usr/bin/env python3
"""
Main orchestrator script for the NRPS C-A-T Module Sequence Database Generation Pipeline.
"""

import os
import sys
import argparse
import logging
import yaml
import time
from pathlib import Path
import subprocess

# Import pipeline modules
from scripts.setup import setup_environment
from scripts.data_acquisition import (
    parse_seed_proteins, 
    fetch_seed_proteins_data, 
    filter_cat_domains, 
    save_sequences_to_fasta,
    retrieve_sequences_from_fasta
)
from scripts.validation import (
    validate_cat_order, 
    parse_hmmscan_domtblout, 
    extract_cat_module_sequences
)
from scripts.alignment import (
    run_mafft_alignment, 
    build_hmm, 
    run_hmmsearch, 
    run_hmmscan
)
from scripts.clustering import (
    run_mmseqs_cluster, 
    run_mmseqs_search, 
    remove_redundancy, 
    parse_cluster_file, 
    parse_network_file, 
    create_network_graph, 
    plot_network
)
from scripts.annotation import (
    annotate_modules, 
    save_annotations_to_tsv, 
    save_fasta_with_annotations
)

def setup_logging(log_dir: str, log_level: str = "INFO") -> None:
    """
    Set up logging configuration.
    
    Args:
        log_dir: Directory to store log files.
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
    """
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    log_file = os.path.join(log_dir, f"pipeline_{timestamp}.log")
    
    # Configure root logger
    log_level_num = getattr(logging, log_level.upper())
    
    # Configure both file and console logging
    logging.basicConfig(
        level=log_level_num,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    logging.info(f"Logging initialized. Log file: {log_file}")

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    
    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="NRPS C-A-T Module Sequence Database Generation Pipeline")
    
    parser.add_argument("--uniref90", required=True, help="Path to the UniRef90 database FASTA file")
    parser.add_argument("--config", default="config.yaml", help="Path to the configuration file")
    parser.add_argument("--seed_proteins", help="Path to a file containing seed UniProt IDs (overrides config)")
    parser.add_argument("--output_dir", help="Directory to store results (overrides config)")
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], 
                         help="Logging level")
    
    # Add restart options
    parser.add_argument("--start_phase", type=int, default=1, choices=[1, 2, 3, 4, 5],
                        help="Phase to start execution from (1-5)")
    parser.add_argument("--end_phase", type=int, default=5, choices=[1, 2, 3, 4, 5],
                        help="Phase to end execution at (1-5)")
    parser.add_argument("--seed_cat_fasta", help="Path to seed C-A-T modules FASTA file (required if start_phase > 1)")
    parser.add_argument("--use_custom_seed", action="store_true", help="Use a custom seed C-A-T modules FASTA file directly for phase 2 without requiring other phase 1 outputs")
    parser.add_argument("--seed_cat_msa", help="Path to seed alignment file (required if start_phase > 1)")
    parser.add_argument("--cat_hmm", help="Path to HMM file (required if start_phase > 1)")
    parser.add_argument("--domtblout_file", help="Path to HMMER domtblout file (required if start_phase > 2)")
    parser.add_argument("--filtered_validated_fasta", help="Path to filtered validated FASTA file (required if start_phase > 3)")
    parser.add_argument("--final_fasta", help="Path to final FASTA file (required if start_phase > 3)")
    parser.add_argument("--annotations_tsv", help="Path to annotations TSV file (required if start_phase > 4)")
    
    return parser.parse_args()

def load_config(config_path: str) -> dict:
    """
    Load configuration from a YAML file.
    
    Args:
        config_path: Path to the configuration file.
        
    Returns:
        Configuration dictionary.
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config

def update_config_with_args(config: dict, args: argparse.Namespace) -> dict:
    """
    Update configuration with command-line arguments.
    
    Args:
        config: Configuration dictionary.
        args: Parsed command-line arguments.
        
    Returns:
        Updated configuration dictionary.
    """
    # Update UniRef90 path
    config['paths']['uniref90'] = args.uniref90
    
    # Update seed proteins path if specified
    if args.seed_proteins:
        config['paths']['seed_proteins'] = args.seed_proteins
    
    # Update output directory if specified
    if args.output_dir:
        config['paths']['output_dir'] = args.output_dir
    
    return config

def run_phase1(config: dict) -> tuple:
    """
    Run Phase 1: Setup, Seed Data Acquisition, and Seed Alignment.
    
    Args:
        config: Configuration dictionary.
        
    Returns:
        Tuple containing:
        - Path to the seed C-A-T modules FASTA file.
        - Path to the seed alignment file.
        - Path to the HMM file.
    """
    logging.info("Starting Phase 1: Setup, Seed Data Acquisition, and Seed Alignment")
    
    # Parse seed protein IDs
    uniprot_ids = parse_seed_proteins(config['paths']['seed_proteins'])
    
    # Fetch seed protein sequences and domain information
    sequences, domains = fetch_seed_proteins_data(uniprot_ids)
    
    # Filter domains to keep only C, A, T domains
    cat_domains = filter_cat_domains(domains)
    
    # Extract C-A-T modules
    max_gap_allowed = config['parameters']['max_gap_allowed']
    seed_cat_modules = extract_cat_module_sequences(cat_domains, sequences, max_gap_allowed)
    
    # Save seed C-A-T modules to FASTA
    seed_cat_fasta = os.path.join(config['paths']['output_dir'], "seed_cat_modules.fasta")
    save_sequences_to_fasta(seed_cat_modules, seed_cat_fasta)
    
    assert os.path.exists(seed_cat_fasta), "Seed C-A-T modules FASTA file not created"
    
    # Run multiple sequence alignment
    seed_cat_msa = os.path.join(config['paths']['output_dir'], "seed_cat_msa.sto")
    success = run_mafft_alignment(seed_cat_fasta, seed_cat_msa, config['parameters']['threads'])
    
    assert success, "MAFFT alignment failed"
    assert os.path.exists(seed_cat_msa), "MSA file not created"
    
    # Build HMM
    cat_hmm = os.path.join(config['paths']['output_dir'], "cat_module.hmm")
    success = build_hmm(seed_cat_msa, cat_hmm)
    
    assert success, "HMM building failed"
    assert os.path.exists(cat_hmm), "HMM file not created"
    
    logging.info("Phase 1 completed successfully")
    
    return seed_cat_fasta, seed_cat_msa, cat_hmm

def run_phase2(config: dict, cat_hmm: str) -> str:
    """
    Run Phase 2: pHMM Construction and Database Search.
    
    Args:
        config: Configuration dictionary.
        cat_hmm: Path to the CAT module HMM file.
        
    Returns:
        Path to the HMMER domtblout file.
    """
    logging.info("Starting Phase 2: pHMM Construction and Database Search")
    
    # Run HMMER search against UniRef90
    uniref90_path = config['paths']['uniref90']
    output_prefix = os.path.join(config['paths']['output_dir'], "uniref90_hits")
    e_value = config['parameters']['hmmer_evalue']
    threads = config['parameters']['threads']
    
    domtblout_file = run_hmmsearch(cat_hmm, uniref90_path, output_prefix, e_value, threads)
    
    assert domtblout_file, "HMMER search failed"
    assert os.path.exists(domtblout_file), "HMMER domtblout file not created"
    
    logging.info(f"Phase 2 completed successfully, results in {domtblout_file}")
    
    return domtblout_file

def run_phase3(config: dict, domtblout_file: str) -> tuple:
    """
    Run Phase 3: Hit Filtering, Validation, and Sequence Retrieval.
    
    Args:
        config: Configuration dictionary.
        domtblout_file: Path to the HMMER domtblout file.
        
    Returns:
        Tuple containing:
        - Path to the filtered validated C-A-T modules FASTA file.
        - Path to the final C-A-T modules FASTA file.
    """
    logging.info("Starting Phase 3: Hit Filtering, Validation, and Sequence Retrieval")
    
    # Parse HMMER results
    e_value_threshold = float(config['parameters']['hmmer_evalue'])
    hits_by_protein = parse_hmmscan_domtblout(domtblout_file, e_value_threshold)
    
    # Count total hits and proteins with multiple hits
    total_proteins = len(hits_by_protein)
    total_hits = sum(len(hits) for hits in hits_by_protein.values())
    multi_hit_proteins = sum(1 for hits in hits_by_protein.values() if len(hits) > 1)
    
    logging.info(f"Parsed {total_hits} hits across {total_proteins} proteins from HMMER results")
    logging.info(f"Found {multi_hit_proteins} proteins with multiple hits ({multi_hit_proteins/total_proteins*100:.1f}%)")
    
    # Create a directory for temporary files
    temp_dir = os.path.join(config['paths']['temp_dir'], f"phase3_{int(time.time())}")
    os.makedirs(temp_dir, exist_ok=True)
    
    # Initialize variables
    sequences_by_id = {}
    max_gap_allowed = config['parameters']['max_gap_allowed']
    validated_modules = {}
    filtered_validated_fasta = None
    final_fasta = None
    
    try:
        # Get all protein IDs that need sequences
        protein_ids = list(hits_by_protein.keys())
        logging.info(f"Retrieving sequences for {len(protein_ids)} proteins from UniRef90")
        
        # Set batch size based on available memory
        available_mem_gb = config['parameters'].get('available_memory_gb', 8)
        batch_size = min(10000, max(1000, int(available_mem_gb * 1000)))
        if config['parameters'].get('batch_size_override'):
            batch_size = config['parameters']['batch_size_override']
            logging.info(f"Using manually overridden batch size of {batch_size}")
        else:
            logging.info(f"Using batch size of {batch_size} sequences based on available memory")
        
        # Optimize thread count based on system resources
        thread_count = config['parameters'].get('threads', 4)
        if thread_count > 8 and total_proteins > 10000:
            # For very large datasets with many cores, limit threads to avoid I/O bottlenecks
            effective_threads = min(thread_count, 12)
            if effective_threads != thread_count:
                logging.info(f"Limiting thread count from {thread_count} to {effective_threads} for optimal I/O performance")
                thread_count = effective_threads
        
        # Retrieve sequences from the UniRef90 database
        uniref90_file = config['paths']['uniref90']
        
        # Check if UniRef90 file exists
        if not os.path.exists(uniref90_file):
            logging.error(f"UniRef90 FASTA file not found: {uniref90_file}")
            raise FileNotFoundError(f"UniRef90 FASTA file not found: {uniref90_file}")
        
        # Log the file size before retrieval to estimate memory requirements
        file_size_gb = os.path.getsize(uniref90_file) / (1024**3)
        logging.info(f"UniRef90 file size: {file_size_gb:.2f} GB")
        
        # Choose retrieval method based on config
        retrieval_method = config['parameters'].get('retrieval_method', 'auto')
        
        if retrieval_method == 'direct' or (retrieval_method == 'auto' and file_size_gb < 10):
            # For smaller files, just use Bio.SeqIO to parse directly - simpler but uses more memory
            logging.info(f"Using direct Bio.SeqIO parsing for the FASTA file (size: {file_size_gb:.2f}GB)")
            try:
                from Bio import SeqIO
                import io
                
                # Create pattern for faster ID matching
                logging.info("Creating ID lookup patterns for direct matching")
                id_patterns = {}
                for pid in protein_ids:
                    id_patterns[pid] = True
                    if pid.startswith('UniRef90_'):
                        id_patterns[pid[9:]] = True
                    else:
                        id_patterns[f"UniRef90_{pid}"] = True
                
                # Parse FASTA file in chunks to avoid memory issues
                logging.info("Starting direct FASTA parsing")
                seq_count = 0
                with open(uniref90_file, 'r') as fasta_handle:
                    batch = []
                    for record in SeqIO.parse(fasta_handle, "fasta"):
                        # Check if this ID matches any of our target patterns
                        record_id = record.id
                        base_id = record_id.split('|')[0] if '|' in record_id else record_id
                        
                        if record_id in id_patterns or base_id in id_patterns:
                            # Store with original ID
                            sequences_by_id[record_id] = str(record.seq)
                            
                            # Also store with base ID if different
                            if base_id != record_id:
                                sequences_by_id[base_id] = str(record.seq)
                            
                            # Handle UniRef90 prefix variants
                            if record_id.startswith('UniRef90_'):
                                sequences_by_id[record_id[9:]] = str(record.seq)
                            else:
                                sequences_by_id[f"UniRef90_{record_id}"] = str(record.seq)
                            
                            seq_count += 1
                            if seq_count % 1000 == 0:
                                logging.info(f"Parsed {seq_count} matching sequences so far")
                
                logging.info(f"Direct parsing completed - found {len(sequences_by_id)} sequences")
                
            except Exception as e:
                logging.error(f"Direct parsing failed: {e}")
                logging.info("Falling back to chunked retrieval")
                # Clear any partial results
                sequences_by_id = {}
                
                # Process in optimally sized batches if the dataset is very large
                if total_proteins > 50000:
                    logging.info(f"Large protein set ({total_proteins} IDs) - processing in batches")
                    
                    for i in range(0, len(protein_ids), batch_size):
                        batch_ids = protein_ids[i:i+batch_size]
                        logging.info(f"Processing batch {i//batch_size + 1}/{(len(protein_ids)-1)//batch_size + 1} ({len(batch_ids)} IDs)")
                        
                        batch_sequences = retrieve_sequences_from_fasta(
                            batch_ids, 
                            uniref90_file,
                            chunk_size=1000,
                            temp_dir=temp_dir,
                            num_processes=thread_count
                        )
                        
                        # Update sequence dictionary and log stats
                        prev_size = len(sequences_by_id)
                        sequences_by_id.update(batch_sequences)
                        logging.info(f"Batch added {len(sequences_by_id) - prev_size} new sequences (total: {len(sequences_by_id)})")
                        
                        # Free memory after each batch
                        batch_sequences = None
                        import gc
                        gc.collect()
                else:
                    # Get sequences for all protein IDs at once if the dataset is smaller
                    sequences_by_id = retrieve_sequences_from_fasta(
                        protein_ids, 
                        uniref90_file,
                        chunk_size=1000,
                        temp_dir=temp_dir,
                        num_processes=thread_count
                    )
        
        else:
            # Process in optimally sized batches if the dataset is very large
            if total_proteins > 50000:
                logging.info(f"Large protein set ({total_proteins} IDs) - processing in batches")
                
                for i in range(0, len(protein_ids), batch_size):
                    batch_ids = protein_ids[i:i+batch_size]
                    logging.info(f"Processing batch {i//batch_size + 1}/{(len(protein_ids)-1)//batch_size + 1} ({len(batch_ids)} IDs)")
                    
                    batch_sequences = retrieve_sequences_from_fasta(
                        batch_ids, 
                        uniref90_file,
                        chunk_size=1000,
                        temp_dir=temp_dir,
                        num_processes=thread_count
                    )
                    
                    # Update sequence dictionary and log stats
                    prev_size = len(sequences_by_id)
                    sequences_by_id.update(batch_sequences)
                    logging.info(f"Batch added {len(sequences_by_id) - prev_size} new sequences (total: {len(sequences_by_id)})")
                    
                    # Free memory after each batch
                    batch_sequences = None
                    import gc
                    gc.collect()
            else:
                # Get sequences for all protein IDs at once if the dataset is smaller
                sequences_by_id = retrieve_sequences_from_fasta(
                    protein_ids, 
                    uniref90_file,
                    chunk_size=1000,
                    temp_dir=temp_dir,
                    num_processes=thread_count
                )
        
        # If we still have no sequences, try fallback method one last time
        if len(sequences_by_id) == 0:
            logging.warning("All retrieval methods failed. Trying basic fallback method...")
            try:
                # Basic retrieval with grep as a last resort
                sequences_by_id = {}
                for i in range(0, len(protein_ids), 100):
                    chunk = protein_ids[i:i+100]
                    for pid in chunk:
                        pattern = f"grep -A 1 -m 1 '{pid}' {uniref90_file} || true"
                        try:
                            result = subprocess.run(pattern, shell=True, stdout=subprocess.PIPE, text=True)
                            lines = result.stdout.strip().split('\n')
                            if len(lines) >= 2 and lines[0].startswith('>'):
                                header = lines[0]
                                seq = lines[1]
                                sequences_by_id[pid] = seq
                        except Exception as e:
                            logging.error(f"Error retrieving sequence {pid}: {e}")
                    
                    if i % 1000 == 0:
                        logging.info(f"Fallback retrieved {len(sequences_by_id)} sequences so far")
                
                logging.info(f"Fallback retrieval completed - found {len(sequences_by_id)} sequences")
            except Exception as e:
                logging.error(f"Fallback retrieval failed: {e}")
        
        # Log retrieval statistics
        missing_sequences = [pid for pid in protein_ids if pid not in sequences_by_id]
        missing_count = len(missing_sequences)
        if missing_count > 0:
            missing_percent = missing_count / len(protein_ids) * 100
            logging.warning(f"Unable to retrieve {missing_count} sequences ({missing_percent:.1f}%)")
            if missing_count < 10:
                logging.warning(f"Missing IDs: {', '.join(missing_sequences)}")
            elif config.get('log_level', 'INFO').upper() == 'DEBUG':
                missing_file = os.path.join(temp_dir, "missing_sequences.txt")
                with open(missing_file, 'w') as f:
                    for pid in missing_sequences:
                        f.write(f"{pid}\n")
                logging.debug(f"List of missing sequences saved to {missing_file}")
        
        if len(sequences_by_id) == 0:
            logging.error("Failed to retrieve any sequences. Check UniRef90 file format and protein IDs.")
            raise RuntimeError("No sequences retrieved after multiple attempts. Cannot proceed.")
        
        # Process hits to extract valid modules
        logging.info(f"Extracting valid C-A-T modules from {len(sequences_by_id)} sequences")
        
        # Extract valid C-A-T modules
        validated_modules = extract_cat_module_sequences(hits_by_protein, sequences_by_id, max_gap_allowed)
        
        # Log extraction statistics
        module_count = len(validated_modules)
        hits_with_sequences = sum(1 for pid in hits_by_protein if pid in sequences_by_id)
        success_rate = module_count / hits_with_sequences * 100 if hits_with_sequences > 0 else 0
        
        logging.info(f"Extracted {module_count} validated C-A-T modules from {hits_with_sequences} proteins with sequences")
        logging.info(f"Extraction success rate: {success_rate:.1f}%")
        
        # Save validated modules to FASTA
        filtered_validated_fasta = os.path.join(config['paths']['output_dir'], "filtered_validated_cat_modules.fasta")
        save_sequences_to_fasta(validated_modules, filtered_validated_fasta)
        
        assert os.path.exists(filtered_validated_fasta), "Filtered validated modules FASTA file not created"
        
        # Remove redundancy with optimized settings
        final_fasta = os.path.join(config['paths']['output_dir'], "final_cat_modules.fasta")
        success = remove_redundancy(
            filtered_validated_fasta,
            final_fasta,
            temp_dir,
            min_seq_id=0.99,
            threads=thread_count
        )
        
        assert success, "Redundancy removal failed"
        assert os.path.exists(final_fasta), "Final modules FASTA file not created"
        
        # Count sequences in final FASTA
        from Bio import SeqIO
        final_count = sum(1 for _ in SeqIO.parse(final_fasta, "fasta"))
        redundant_count = module_count - final_count
        redundant_percent = redundant_count / module_count * 100 if module_count > 0 else 0
        
        logging.info(f"Removed {redundant_count} redundant sequences ({redundant_percent:.1f}%)")
        logging.info(f"Final database contains {final_count} unique C-A-T modules")
        
        logging.info(f"Phase 3 completed successfully, results in {final_fasta}")
        
    except Exception as e:
        logging.exception(f"Phase 3 failed: {e}")
        raise
        
    finally:
        # Clean up temporary directory if requested in config
        if config.get('parameters', {}).get('cleanup_temp', True):
            import shutil
            logging.info(f"Cleaning up temporary directory: {temp_dir}")
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                logging.warning(f"Failed to clean up temporary directory: {e}")
    
    return filtered_validated_fasta, final_fasta

def run_phase4(config: dict, final_fasta: str, hits_by_protein: dict) -> str:
    """
    Run Phase 4: Annotation.
    
    Args:
        config: Configuration dictionary.
        final_fasta: Path to the final C-A-T modules FASTA file.
        hits_by_protein: Dictionary mapping protein IDs to domain hits.
        
    Returns:
        Path to the annotations TSV file.
    """
    logging.info("Starting Phase 4: Annotation")
    
    # For this example, we'll use a simplified annotation process
    # In a real implementation, you would fetch annotations from UniProt
    
    # Read the final sequences
    from Bio import SeqIO
    sequences = {}
    for record in SeqIO.parse(final_fasta, "fasta"):
        sequences[record.id] = str(record.seq)
    
    # Annotate modules
    annotations = annotate_modules(sequences, hits_by_protein)
    
    # Save annotations
    annotations_tsv = os.path.join(config['paths']['output_dir'], "final_cat_modules_annotations.tsv")
    success = save_annotations_to_tsv(annotations, annotations_tsv)
    
    assert success, "Saving annotations failed"
    assert os.path.exists(annotations_tsv), "Annotations TSV file not created"
    
    # Save annotated FASTA
    annotated_fasta = os.path.join(config['paths']['output_dir'], "final_cat_modules_annotated.fasta")
    success = save_fasta_with_annotations(sequences, annotations, annotated_fasta)
    
    assert success, "Saving annotated FASTA failed"
    assert os.path.exists(annotated_fasta), "Annotated FASTA file not created"
    
    logging.info(f"Phase 4 completed successfully, annotations in {annotations_tsv}")
    
    return annotations_tsv

def run_phase5(config: dict, final_fasta: str, annotations_tsv: str) -> tuple:
    """
    Run Phase 5: Clustering and Network Visualization.
    
    Args:
        config: Configuration dictionary.
        final_fasta: Path to the final C-A-T modules FASTA file.
        annotations_tsv: Path to the annotations TSV file.
        
    Returns:
        Tuple containing:
        - Path to the cluster file.
        - Path to the network file.
        - Path to the network visualization file.
    """
    logging.info("Starting Phase 5: Clustering and Network Visualization")
    
    # Create a directory for temporary files
    temp_dir = config['paths']['temp_dir']
    os.makedirs(temp_dir, exist_ok=True)
    
    # Run sequence clustering
    output_prefix = os.path.join(config['paths']['output_dir'], "cluster_results")
    min_seq_id = config['parameters']['cluster_min_seq_id']
    coverage = config['parameters']['cluster_cov']
    cov_mode = config['parameters']['cluster_cov_mode']
    threads = config['parameters']['threads']
    
    cluster_file = run_mmseqs_cluster(
        final_fasta,
        output_prefix,
        temp_dir,
        min_seq_id,
        coverage,
        cov_mode,
        threads
    )
    
    assert cluster_file, "Clustering failed"
    assert os.path.exists(cluster_file), "Cluster file not created"
    
    # Run all-against-all search for network generation
    network_prefix = os.path.join(config['paths']['output_dir'], "network")
    network_file = run_mmseqs_search(
        final_fasta,
        network_prefix,
        temp_dir,
        sensitivity=7.0,
        threads=threads
    )
    
    assert network_file, "Network search failed"
    assert os.path.exists(network_file), "Network file not created"
    
    # Parse network file
    min_identity = config['parameters']['cluster_min_seq_id']
    edges = parse_network_file(network_file, min_identity, config['parameters']['hmmer_evalue'])
    
    # Read annotations
    import pandas as pd
    annotations = {}
    try:
        df = pd.read_csv(annotations_tsv, sep='\t')
        for _, row in df.iterrows():
            module_id = row['module_id']
            annotations[module_id] = row.to_dict()
    except Exception as e:
        logging.error(f"Failed to read annotations: {e}")
    
    # Create network graph
    G = create_network_graph(edges, annotations)
    
    # Plot network
    network_viz = os.path.join(config['paths']['output_dir'], "network_visualization.png")
    success = plot_network(
        G,
        network_viz,
        color_attribute='taxonomy',  # Assuming this is in the annotations
        size_attribute='module_length'  # Assuming this is in the annotations
    )
    
    logging.info(f"Phase 5 completed successfully, network visualization in {network_viz}")
    
    return cluster_file, network_file, network_viz

def main():
    """Main function."""
    # Parse command-line arguments
    args = parse_arguments()
    
    # Load configuration
    config = load_config(args.config)
    
    # Update configuration with command-line arguments
    config = update_config_with_args(config, args)
    
    # Set up logging
    setup_logging(config['paths']['log_dir'], args.log_level)
    
    # Log start time
    start_time = time.time()
    logging.info(f"Starting NRPS C-A-T Module Database Generation Pipeline (Phases {args.start_phase}-{args.end_phase})")
    
    try:
        # Set up the environment
        if not setup_environment(config):
            logging.error("Environment setup failed. Exiting.")
            sys.exit(1)
        
        # Initialize variables that may be set in different phases
        seed_cat_fasta = args.seed_cat_fasta
        seed_cat_msa = args.seed_cat_msa
        cat_hmm = args.cat_hmm
        domtblout_file = args.domtblout_file
        filtered_validated_fasta = args.filtered_validated_fasta
        final_fasta = args.final_fasta
        annotations_tsv = args.annotations_tsv
        
        # Validate required files based on start phase
        if args.start_phase > 1:
            # Check if custom seed FASTA option is used
            if args.use_custom_seed:
                if not seed_cat_fasta:
                    logging.error("When using custom seed FASTA, you must provide --seed_cat_fasta")
                    sys.exit(1)
                if not os.path.exists(seed_cat_fasta):
                    logging.error("Custom seed FASTA file does not exist")
                    sys.exit(1)
                
                logging.info(f"Using custom seed FASTA file: {seed_cat_fasta}")
                
                # Create the MSA and HMM from the custom seed FASTA
                os.makedirs(config['paths']['output_dir'], exist_ok=True)
                
                seed_cat_msa = os.path.join(config['paths']['output_dir'], "seed_cat_msa.sto")
                logging.info(f"Creating MSA from custom seed FASTA: {seed_cat_msa}")
                success = run_mafft_alignment(seed_cat_fasta, seed_cat_msa, config['parameters']['threads'])
                
                if not success:
                    logging.error("Failed to create MSA from custom seed FASTA")
                    sys.exit(1)
                
                cat_hmm = os.path.join(config['paths']['output_dir'], "cat_module.hmm")
                logging.info(f"Creating HMM from MSA: {cat_hmm}")
                success = build_hmm(seed_cat_msa, cat_hmm)
                
                if not success:
                    logging.error("Failed to create HMM from MSA")
                    sys.exit(1)
            else:
                if not all([seed_cat_fasta, seed_cat_msa, cat_hmm]):
                    logging.error("When starting from phase > 1, you must provide --seed_cat_fasta, --seed_cat_msa, and --cat_hmm")
                    sys.exit(1)
                if not all(os.path.exists(f) for f in [seed_cat_fasta, seed_cat_msa, cat_hmm]):
                    logging.error("One or more input files from phase 1 do not exist")
                    sys.exit(1)
        
        if args.start_phase > 2:
            if not domtblout_file:
                logging.error("When starting from phase > 2, you must provide --domtblout_file")
                sys.exit(1)
            if not os.path.exists(domtblout_file):
                logging.error("HMMER domtblout file does not exist")
                sys.exit(1)
                
        if args.start_phase > 3:
            if not all([filtered_validated_fasta, final_fasta]):
                logging.error("When starting from phase > 3, you must provide --filtered_validated_fasta and --final_fasta")
                sys.exit(1)
            if not all(os.path.exists(f) for f in [filtered_validated_fasta, final_fasta]):
                logging.error("One or more input files from phase 3 do not exist")
                sys.exit(1)
                
        if args.start_phase > 4:
            if not annotations_tsv:
                logging.error("When starting from phase > 4, you must provide --annotations_tsv")
                sys.exit(1)
            if not os.path.exists(annotations_tsv):
                logging.error("Annotations TSV file does not exist")
                sys.exit(1)
        
        # Phase 1: Setup, Seed Data Acquisition, and Seed Alignment
        if args.start_phase <= 1 and args.end_phase >= 1:
            seed_cat_fasta, seed_cat_msa, cat_hmm = run_phase1(config)
        else:
            logging.info(f"Skipping phase 1, using provided files: {seed_cat_fasta}, {seed_cat_msa}, {cat_hmm}")
        
        # Phase 2: pHMM Construction and Database Search
        if args.start_phase <= 2 and args.end_phase >= 2:
            domtblout_file = run_phase2(config, cat_hmm)
        else:
            logging.info(f"Skipping phase 2, using provided domtblout file: {domtblout_file}")
        
        # Phase 3: Hit Filtering, Validation, and Sequence Retrieval
        if args.start_phase <= 3 and args.end_phase >= 3:
            filtered_validated_fasta, final_fasta = run_phase3(config, domtblout_file)
        else:
            logging.info(f"Skipping phase 3, using provided files: {filtered_validated_fasta}, {final_fasta}")
        
        # Parse HMMER results for annotation (needed for phases 4 and 5)
        if args.start_phase <= 4 and args.end_phase >= 4:
            e_value_threshold = float(config['parameters']['hmmer_evalue'])
            hits_by_protein = parse_hmmscan_domtblout(domtblout_file, e_value_threshold)
            
            # Phase 4: Annotation
            annotations_tsv = run_phase4(config, final_fasta, hits_by_protein)
        else:
            logging.info(f"Skipping phase 4, using provided annotations file: {annotations_tsv}")
        
        # Phase 5: Clustering and Network Visualization
        if args.start_phase <= 5 and args.end_phase >= 5:
            cluster_file, network_file, network_viz = run_phase5(config, final_fasta, annotations_tsv)
            logging.info(f"Generated clustering and network files: {cluster_file}, {network_file}, {network_viz}")
        else:
            logging.info("Skipping phase 5")
        
        # Log completion time
        end_time = time.time()
        elapsed_time = end_time - start_time
        logging.info(f"Pipeline completed successfully (phases {args.start_phase}-{args.end_phase}) in {elapsed_time:.2f} seconds")
        
    except Exception as e:
        logging.exception(f"Pipeline failed during phase execution: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 