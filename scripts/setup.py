#!/usr/bin/env python3
"""
Setup module for environment checks and directory creation.
"""

import os
import sys
import logging
import subprocess
from pathlib import Path

def check_dependencies():
    """
    Check if required dependencies are installed.
    
    Returns:
        bool: True if all dependencies are installed, False otherwise.
    """
    required_commands = {
        'python': 'python --version',
        'hmmbuild': 'hmmbuild -h',
        'hmmsearch': 'hmmsearch -h',
        'hmmscan': 'hmmscan -h',
        'mmseqs': 'mmseqs -h',
        'mafft': 'mafft --version'
    }
    
    missing = []
    
    for cmd, check_cmd in required_commands.items():
        try:
            subprocess.run(
                check_cmd.split(), 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                check=True
            )
            logging.info(f"{cmd} is installed.")
        except (subprocess.SubprocessError, FileNotFoundError):
            missing.append(cmd)
            logging.error(f"{cmd} is not installed or not in PATH.")
    
    if missing:
        logging.error(f"Missing dependencies: {', '.join(missing)}")
        return False
    
    return True

def create_directories(config):
    """
    Create necessary directories if they don't exist.
    
    Args:
        config (dict): Configuration dictionary.
        
    Returns:
        bool: True if directories were created successfully.
    """
    directories = [
        config['paths']['output_dir'],
        config['paths']['temp_dir'],
        config['paths']['log_dir'],
        os.path.dirname(config['paths']['seed_proteins']),
        config['paths']['pfam_hmms']
    ]
    
    for directory in directories:
        Path(directory).mkdir(parents=True, exist_ok=True)
        logging.info(f"Directory created or already exists: {directory}")
    
    return True

def check_input_files(config):
    """
    Check if required input files exist.
    
    Args:
        config (dict): Configuration dictionary.
        
    Returns:
        bool: True if all required files exist.
    """
    seed_proteins_path = config['paths']['seed_proteins']
    
    if not os.path.isfile(seed_proteins_path):
        logging.error(f"Seed proteins file not found: {seed_proteins_path}")
        return False
    
    uniref90_path = config['paths']['uniref90']
    if not os.path.isfile(uniref90_path):
        logging.error(f"UniRef90 database file not found: {uniref90_path}")
        return False
    
    return True

def setup_environment(config: dict) -> bool:
    """
    Set up the environment for the pipeline.
    
    Args:
        config: Configuration dictionary.
        
    Returns:
        True if setup was successful, False otherwise.
    """
    try:
        # Create output directory
        output_dir = config['paths']['output_dir']
        os.makedirs(output_dir, exist_ok=True)
        logging.info(f"Created output directory: {output_dir}")
        
        # Create temporary directory
        temp_dir = config['paths']['temp_dir']
        os.makedirs(temp_dir, exist_ok=True)
        logging.info(f"Created temporary directory: {temp_dir}")
        
        # Check for required tools
        required_tools = ["hmmsearch", "hmmscan", "hmmbuild", "mafft", "mmseqs"]
        missing_tools = []
        
        for tool in required_tools:
            try:
                subprocess.run(["which", tool], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                logging.info(f"Found required tool: {tool}")
            except subprocess.CalledProcessError:
                logging.error(f"Required tool not found: {tool}")
                missing_tools.append(tool)
        
        if missing_tools:
            logging.error(f"Missing required tools: {', '.join(missing_tools)}")
            logging.error("Please install the missing tools and try again.")
            return False
        
        # Check UniRef90 database
        uniref90_file = config['paths']['uniref90']
        if not os.path.exists(uniref90_file):
            logging.error(f"UniRef90 database not found: {uniref90_file}")
            return False
        
        # Validate UniRef90 file format
        if not validate_uniref90_file(uniref90_file):
            logging.error(f"UniRef90 database format validation failed: {uniref90_file}")
            return False
            
        logging.info(f"Found UniRef90 database: {uniref90_file}")
        
        # Check Pfam HMM database if specified
        if 'pfam_hmm' in config['paths'] and config['paths']['pfam_hmm']:
            pfam_hmm = config['paths']['pfam_hmm']
            if not check_hmm_database(pfam_hmm):
                logging.error(f"Pfam HMM database validation failed: {pfam_hmm}")
                return False
                
        logging.info("Environment setup completed successfully")
        return True
        
    except Exception as e:
        logging.exception(f"Environment setup failed: {e}")
        return False

def validate_uniref90_file(uniref90_file: str, sample_size: int = 10) -> bool:
    """
    Validate the UniRef90 file format.
    
    Args:
        uniref90_file: Path to the UniRef90 FASTA file.
        sample_size: Number of sequence headers to sample.
        
    Returns:
        True if the file format is valid, False otherwise.
    """
    try:
        # Check if the file exists and has a reasonable size
        if not os.path.exists(uniref90_file):
            logging.error(f"UniRef90 file not found: {uniref90_file}")
            return False
            
        file_size = os.path.getsize(uniref90_file)
        if file_size < 1000000:  # 1MB minimum
            logging.warning(f"UniRef90 file is suspiciously small: {file_size} bytes")
            
        # Check if the file starts with a FASTA header
        with open(uniref90_file, 'r') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('>'):
                logging.error(f"UniRef90 file does not start with a FASTA header (>): {first_line}")
                return False
                
            # Check sample headers
            headers = []
            for _ in range(sample_size):
                line = f.readline().strip()
                while line and not line.startswith('>'):
                    line = f.readline().strip()
                if line:
                    headers.append(line)
                    
            # Log sample headers for debugging
            if headers:
                logging.info(f"Sample UniRef90 headers:")
                for header in headers[:5]:
                    logging.info(f"  {header}")
                    
                # Check if headers match UniRef90 format
                valid_headers = [h for h in headers if 'UniRef90_' in h]
                if not valid_headers:
                    logging.warning("No UniRef90_ prefix found in sample headers. This may cause issues with ID matching.")
                    
        return True
            
    except Exception as e:
        logging.exception(f"Error validating UniRef90 file: {e}")
        return False

def check_hmm_database(hmm_db: str) -> bool:
    """
    Check if the HMM database exists and is properly formatted.
    
    Args:
        hmm_db: Path to the HMM database.
        
    Returns:
        True if the database is valid, False otherwise.
    """
    try:
        # Check if the file exists
        if not os.path.exists(hmm_db):
            logging.error(f"HMM database not found: {hmm_db}")
            return False
            
        # Check if the database is properly formatted
        cmd = ["hmmstat", hmm_db]
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # If we get here, the database is valid
        logging.info(f"HMM database validated: {hmm_db}")
        return True
        
    except subprocess.CalledProcessError as e:
        logging.error(f"HMM database validation failed: {e}")
        logging.error(f"STDERR: {e.stderr.decode('utf-8')}")
        return False
    except Exception as e:
        logging.exception(f"Error checking HMM database: {e}")
        return False 