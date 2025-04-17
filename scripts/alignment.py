#!/usr/bin/env python3
"""
Alignment module for sequence alignment and HMM operations.
"""

import os
import logging
import subprocess
from typing import Optional

def run_mafft_alignment(input_fasta: str, output_sto: str, threads: int = 1) -> bool:
    """
    Run MAFFT for multiple sequence alignment and convert output to Stockholm format.
    
    Args:
        input_fasta: Path to the input FASTA file.
        output_sto: Path to the output Stockholm file.
        threads: Number of CPU threads to use.
        
    Returns:
        True if alignment was successful, False otherwise.
    """
    logging.info(f"Running MAFFT alignment on {input_fasta} with {threads} threads")
    
    # Create a temporary FASTA output file
    temp_output_fasta = f"{output_sto}.temp.fasta"
    
    # MAFFT command
    cmd = ["mafft", "--auto", "--thread", str(threads), input_fasta]
    
    try:
        with open(temp_output_fasta, 'w') as outfile:
            process = subprocess.run(
                cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                check=True,
                text=True
            )
        
        logging.info(f"MAFFT alignment completed successfully")
        
        # Convert FASTA format to Stockholm format
        convert_fasta_to_stockholm(temp_output_fasta, output_sto)
        
        # Clean up temporary file
        os.remove(temp_output_fasta)
        
        return True
        
    except subprocess.CalledProcessError as e:
        logging.error(f"MAFFT alignment failed: {e}")
        logging.error(f"STDERR: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"Error during alignment: {e}")
        return False

def convert_fasta_to_stockholm(input_fasta: str, output_sto: str) -> bool:
    """
    Convert FASTA alignment to Stockholm format.
    
    Args:
        input_fasta: Path to the input FASTA file.
        output_sto: Path to the output Stockholm file.
        
    Returns:
        True if conversion was successful, False otherwise.
    """
    try:
        from Bio import AlignIO
        
        # Read the alignment in FASTA format
        alignment = AlignIO.read(input_fasta, "fasta")
        
        # Write the alignment in Stockholm format
        AlignIO.write(alignment, output_sto, "stockholm")
        
        logging.info(f"Converted FASTA to Stockholm format: {output_sto}")
        return True
        
    except Exception as e:
        logging.error(f"Failed to convert FASTA to Stockholm: {e}")
        return False

def build_hmm(alignment_sto: str, output_hmm: str) -> bool:
    """
    Build a profile HMM from a multiple sequence alignment using HMMER's hmmbuild.
    
    Args:
        alignment_sto: Path to the input alignment in Stockholm format.
        output_hmm: Path to the output HMM file.
        
    Returns:
        True if HMM building was successful, False otherwise.
    """
    logging.info(f"Building HMM from alignment: {alignment_sto}")
    
    # hmmbuild command
    cmd = ["hmmbuild", output_hmm, alignment_sto]
    
    try:
        process = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
        
        logging.info(f"HMM building completed successfully: {output_hmm}")
        return True
        
    except subprocess.CalledProcessError as e:
        logging.error(f"hmmbuild failed: {e}")
        logging.error(f"STDOUT: {e.stdout}")
        logging.error(f"STDERR: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"Error during HMM building: {e}")
        return False

def run_hmmsearch(hmm_file: str, database_fasta: str, output_prefix: str, e_value: float = 1e-10, cpu: int = 1) -> Optional[str]:
    """
    Run hmmsearch to search for matches to the profile HMM in a sequence database.
    
    Args:
        hmm_file: Path to the input HMM file.
        database_fasta: Path to the sequence database in FASTA format.
        output_prefix: Prefix for output files.
        e_value: E-value threshold for reporting matches.
        cpu: Number of CPU threads to use.
        
    Returns:
        Path to the domtblout file if search was successful, None otherwise.
    """
    logging.info(f"Running hmmsearch with HMM {hmm_file} against {database_fasta}")
    
    # Output files
    output_log = f"{output_prefix}.log"
    output_domtblout = f"{output_prefix}.domtblout"
    
    # hmmsearch command
    cmd = [
        "hmmsearch",
        "--domtblout", output_domtblout,
        "-E", str(e_value),
        "--cpu", str(cpu),
        hmm_file, database_fasta
    ]
    
    try:
        with open(output_log, 'w') as logfile:
            process = subprocess.run(
                cmd,
                stdout=logfile,
                stderr=subprocess.PIPE,
                check=True,
                text=True
            )
        
        logging.info(f"hmmsearch completed successfully, results in {output_domtblout}")
        return output_domtblout
        
    except subprocess.CalledProcessError as e:
        logging.error(f"hmmsearch failed: {e}")
        logging.error(f"STDERR: {e.stderr}")
        return None
    except Exception as e:
        logging.error(f"Error during hmmsearch: {e}")
        return None

def run_hmmscan(hmm_db: str, query_fasta: str, output_prefix: str, e_value: float = 1e-10, cpu: int = 1) -> Optional[str]:
    """
    Run hmmscan to search for domains in protein sequences.
    
    Args:
        hmm_db: Path to the HMM database.
        query_fasta: Path to the query sequences in FASTA format.
        output_prefix: Prefix for output files.
        e_value: E-value threshold for reporting matches.
        cpu: Number of CPU threads to use.
        
    Returns:
        Path to the domtblout file if scan was successful, None otherwise.
    """
    logging.info(f"Running hmmscan with HMM database {hmm_db} against {query_fasta}")
    
    # Output files
    output_log = f"{output_prefix}.log"
    output_domtblout = f"{output_prefix}.domtblout"
    
    # hmmscan command
    cmd = [
        "hmmscan",
        "--domtblout", output_domtblout,
        "-E", str(e_value),
        "--cpu", str(cpu),
        hmm_db, query_fasta
    ]
    
    try:
        with open(output_log, 'w') as logfile:
            process = subprocess.run(
                cmd,
                stdout=logfile,
                stderr=subprocess.PIPE,
                check=True,
                text=True
            )
        
        logging.info(f"hmmscan completed successfully, results in {output_domtblout}")
        return output_domtblout
        
    except subprocess.CalledProcessError as e:
        logging.error(f"hmmscan failed: {e}")
        logging.error(f"STDERR: {e.stderr}")
        return None
    except Exception as e:
        logging.error(f"Error during hmmscan: {e}")
        return None 