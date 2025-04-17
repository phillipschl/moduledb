#!/usr/bin/env python3
"""
Clustering module for sequence clustering and network generation.
"""

import os
import logging
import subprocess
from typing import Optional, Dict, List, Any, Tuple
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def run_mmseqs_cluster(
    input_fasta: str,
    output_prefix: str,
    temp_dir: str,
    min_seq_id: float = 0.7,
    coverage: float = 0.8,
    cov_mode: int = 0,
    threads: int = 1
) -> Optional[str]:
    """
    Run MMseqs2 easy-cluster on a set of sequences.
    
    Args:
        input_fasta: Path to the input FASTA file.
        output_prefix: Prefix for output files.
        temp_dir: Directory for temporary files.
        min_seq_id: Minimum sequence identity threshold.
        coverage: Coverage threshold.
        cov_mode: Coverage mode (0: bidirectional, 1: query, 2: target).
        threads: Number of CPU threads to use.
        
    Returns:
        Path to the cluster file if clustering was successful, None otherwise.
    """
    logging.info(f"Running MMseqs2 clustering on {input_fasta}")
    
    # Make sure the temporary directory exists
    os.makedirs(temp_dir, exist_ok=True)
    
    # MMseqs2 easy-cluster command
    cmd = [
        "mmseqs", "easy-cluster", 
        input_fasta, 
        output_prefix, 
        temp_dir,
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", str(cov_mode),
        "--threads", str(threads)
    ]
    
    try:
        process = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
        
        cluster_file = f"{output_prefix}_cluster.tsv"
        
        if os.path.exists(cluster_file):
            logging.info(f"MMseqs2 clustering completed successfully, results in {cluster_file}")
            return cluster_file
        else:
            logging.error(f"Expected cluster file {cluster_file} not found")
            return None
        
    except subprocess.CalledProcessError as e:
        logging.error(f"MMseqs2 clustering failed: {e}")
        logging.error(f"STDOUT: {e.stdout}")
        logging.error(f"STDERR: {e.stderr}")
        return None
    except Exception as e:
        logging.error(f"Error during MMseqs2 clustering: {e}")
        return None

def run_mmseqs_search(
    input_fasta: str,
    output_prefix: str,
    temp_dir: str,
    sensitivity: float = 7.0,
    threads: int = 1
) -> Optional[str]:
    """
    Run MMseqs2 all-against-all search on a set of sequences.
    
    Args:
        input_fasta: Path to the input FASTA file.
        output_prefix: Prefix for output files.
        temp_dir: Directory for temporary files.
        sensitivity: Sensitivity parameter (1-9, higher is more sensitive).
        threads: Number of CPU threads to use.
        
    Returns:
        Path to the network edge file if search was successful, None otherwise.
    """
    logging.info(f"Running MMseqs2 all-against-all search on {input_fasta}")
    
    # Make sure the temporary directory exists
    os.makedirs(temp_dir, exist_ok=True)
    
    # Create MMseqs2 DB
    db_prefix = f"{output_prefix}_db"
    
    try:
        # Create the database
        cmd_createdb = ["mmseqs", "createdb", input_fasta, db_prefix]
        subprocess.run(cmd_createdb, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Run all-against-all search
        cmd_search = [
            "mmseqs", "search", 
            db_prefix, db_prefix, 
            f"{output_prefix}_search", 
            temp_dir,
            "-s", str(sensitivity),
            "--threads", str(threads)
        ]
        subprocess.run(cmd_search, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Convert to alignment format
        cmd_convertalis = [
            "mmseqs", "convertalis", 
            db_prefix, db_prefix, 
            f"{output_prefix}_search", 
            f"{output_prefix}_network.m8"
        ]
        subprocess.run(cmd_convertalis, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        network_file = f"{output_prefix}_network.m8"
        
        if os.path.exists(network_file):
            logging.info(f"MMseqs2 search completed successfully, results in {network_file}")
            return network_file
        else:
            logging.error(f"Expected network file {network_file} not found")
            return None
            
    except subprocess.CalledProcessError as e:
        logging.error(f"MMseqs2 search failed: {e}")
        return None
    except Exception as e:
        logging.error(f"Error during MMseqs2 search: {e}")
        return None

def remove_redundancy(
    input_fasta: str,
    output_fasta: str,
    temp_dir: str,
    min_seq_id: float = 0.99,
    threads: int = 1
) -> bool:
    """
    Remove redundant sequences using MMseqs2 linclust.
    
    Args:
        input_fasta: Path to the input FASTA file.
        output_fasta: Path to the output FASTA file.
        temp_dir: Directory for temporary files.
        min_seq_id: Minimum sequence identity threshold.
        threads: Number of CPU threads to use.
        
    Returns:
        True if redundancy removal was successful, False otherwise.
    """
    logging.info(f"Removing redundant sequences from {input_fasta}")
    
    # Make sure the temporary directory exists
    os.makedirs(temp_dir, exist_ok=True)
    
    # MMseqs2 easy-linclust command
    cmd = [
        "mmseqs", "easy-linclust", 
        input_fasta, 
        os.path.splitext(output_fasta)[0], 
        temp_dir,
        "--min-seq-id", str(min_seq_id),
        "--threads", str(threads)
    ]
    
    try:
        process = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
        
        final_fasta = f"{os.path.splitext(output_fasta)[0]}_rep_seq.fasta"
        
        # Rename the output file to the expected name
        if os.path.exists(final_fasta):
            os.rename(final_fasta, output_fasta)
            logging.info(f"MMseqs2 redundancy removal completed successfully, results in {output_fasta}")
            return True
        else:
            logging.error(f"Expected file {final_fasta} not found")
            return False
        
    except subprocess.CalledProcessError as e:
        logging.error(f"MMseqs2 redundancy removal failed: {e}")
        logging.error(f"STDOUT: {e.stdout}")
        logging.error(f"STDERR: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"Error during MMseqs2 redundancy removal: {e}")
        return False

def parse_cluster_file(cluster_file: str) -> Dict[str, List[str]]:
    """
    Parse MMseqs2 cluster file to get cluster memberships.
    
    Args:
        cluster_file: Path to the MMseqs2 cluster file.
        
    Returns:
        Dictionary mapping representative sequence IDs to lists of cluster member IDs.
    """
    clusters = {}
    
    try:
        df = pd.read_csv(cluster_file, sep='\t', header=None, names=['rep', 'member'])
        
        # Group by representative sequence
        for rep, group in df.groupby('rep'):
            clusters[rep] = group['member'].tolist()
            
        logging.info(f"Parsed {len(clusters)} clusters from {cluster_file}")
        return clusters
        
    except Exception as e:
        logging.error(f"Failed to parse cluster file {cluster_file}: {e}")
        return {}

def parse_network_file(
    network_file: str, 
    min_identity: float = 0.7,
    max_evalue: float = 1e-10
) -> List[Tuple[str, str, float]]:
    """
    Parse MMseqs2 network file to get sequence similarities.
    
    Args:
        network_file: Path to the MMseqs2 network file.
        min_identity: Minimum sequence identity threshold for edges.
        max_evalue: Maximum E-value threshold for edges.
        
    Returns:
        List of tuples (query_id, target_id, similarity).
    """
    edges = []
    
    try:
        # Read the network file (BLAST tabular format)
        df = pd.read_csv(
            network_file, 
            sep='\t', 
            header=None, 
            names=['query', 'target', 'identity', 'aln_length', 'mismatches', 
                   'gaps', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits']
        )
        
        # Filter by identity and E-value
        filtered_df = df[(df['identity'] >= min_identity * 100) & (df['evalue'] <= max_evalue)]
        
        # Remove self-matches
        filtered_df = filtered_df[filtered_df['query'] != filtered_df['target']]
        
        # Create edge list
        edges = [(row['query'], row['target'], row['identity'] / 100) 
                 for _, row in filtered_df.iterrows()]
        
        logging.info(f"Parsed {len(edges)} edges from {network_file}")
        return edges
        
    except Exception as e:
        logging.error(f"Failed to parse network file {network_file}: {e}")
        return []

def create_network_graph(
    edges: List[Tuple[str, str, float]], 
    annotations: Dict[str, Dict[str, Any]] = None
) -> nx.Graph:
    """
    Create a NetworkX graph from an edge list.
    
    Args:
        edges: List of tuples (query_id, target_id, similarity).
        annotations: Dictionary mapping sequence IDs to annotation dictionaries.
        
    Returns:
        NetworkX graph object.
    """
    G = nx.Graph()
    
    # Add edges
    for query, target, similarity in edges:
        G.add_edge(query, target, weight=similarity)
    
    # Add node attributes from annotations
    if annotations:
        for node in G.nodes():
            if node in annotations:
                for key, value in annotations[node].items():
                    G.nodes[node][key] = value
    
    logging.info(f"Created network graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    return G

def plot_network(
    G: nx.Graph, 
    output_file: str, 
    color_attribute: str = None, 
    size_attribute: str = None,
    layout = None
) -> bool:
    """
    Plot a network graph and save to file.
    
    Args:
        G: NetworkX graph object.
        output_file: Path to the output image file.
        color_attribute: Node attribute to use for coloring.
        size_attribute: Node attribute to use for sizing.
        layout: Layout algorithm to use (e.g., nx.spring_layout).
        
    Returns:
        True if the plot was created successfully, False otherwise.
    """
    try:
        plt.figure(figsize=(12, 12))
        
        # Set up layout
        if layout is None:
            layout = nx.spring_layout(G)
        elif not callable(layout):
            layout = nx.spring_layout(G)
        else:
            layout = layout(G)
        
        # Set up node sizes
        if size_attribute and any(size_attribute in G.nodes[n] for n in G.nodes()):
            sizes = [G.nodes[n].get(size_attribute, 100) for n in G.nodes()]
        else:
            sizes = 100
        
        # Set up node colors
        if color_attribute and any(color_attribute in G.nodes[n] for n in G.nodes()):
            colors = [G.nodes[n].get(color_attribute, 0) for n in G.nodes()]
            nodes = nx.draw_networkx_nodes(G, layout, node_size=sizes, node_color=colors, cmap=plt.cm.viridis)
            plt.colorbar(nodes)
        else:
            nx.draw_networkx_nodes(G, layout, node_size=sizes)
        
        # Draw edges with alpha based on weight
        edge_weights = [G[u][v].get('weight', 0.1) for u, v in G.edges()]
        nx.draw_networkx_edges(G, layout, alpha=0.3, width=edge_weights)
        
        # Draw labels for larger nodes
        if size_attribute:
            labels = {n: n for n in G.nodes() if G.nodes[n].get(size_attribute, 0) > np.mean(sizes)}
        else:
            labels = {n: n for n in G.nodes()}
        nx.draw_networkx_labels(G, layout, labels=labels, font_size=8)
        
        plt.title("NRPS C-A-T Module Sequence Similarity Network")
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        logging.info(f"Network visualization saved to {output_file}")
        return True
        
    except Exception as e:
        logging.error(f"Failed to create network visualization: {e}")
        return False 