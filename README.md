# NRPS C-A-T Module Sequence Database

A pipeline for generating a high-quality sequence database of continuous C-A-T (Condensation-Adenylation-Thiolation) modules from Non-Ribosomal Peptide Synthetase (NRPS) enzymes.

## Overview

This pipeline:
1. Identifies and extracts sequences representing complete, contiguous C-A-T modules
2. Uses HMMER and MMseqs2 for sequence analysis and clustering
3. Performs robust validation, logging, and annotation
4. Generates visualization networks of sequence relationships

## Installation

### Dependencies

This project uses conda to manage dependencies. To set up the environment:

```bash
# Clone the repository
git clone https://github.com/phillipschl/moduledb.git
cd moduledb

# Create and activate the conda environment
conda env create -f environment.yml
conda activate moduledb

# Install additional performance tools (optional but recommended)
bash install.sh
```

## Performance Optimization

The pipeline includes several optimizations for handling large datasets:

1. **FASTA Index**: The pipeline automatically creates and uses FASTA index files (`.fai`) for faster sequence retrieval
2. **SeqKit Integration**: Uses the ultra-fast SeqKit tool for sequence extraction with automatic index detection
3. **Parallel Processing**: Processes large sequence sets using multiple CPU cores
4. **Ripgrep Support**: Uses ripgrep instead of grep for faster text searching when available
5. **Chunked Processing**: Processes large datasets in optimally sized chunks
6. **Optimized UniRef90 Handling**: Automatically detects UniRef90 databases and uses specialized retrieval methods that work directly with UniRef90 prefixes for optimal performance
7. **SeqKit Version Compatibility**: Automatically detects seqkit capabilities and adapts to work with all versions

To maximize performance:

- Install optional tools via `install.sh`
- Increase the thread count in `config.yaml` based on your system's capabilities
- For very large jobs, pre-index the UniRef90 database: `seqkit faidx uniref90.fasta`
- For systems with limited RAM, reduce the chunk size parameter in `main.py`

### SeqKit Usage Notes

The pipeline works with all versions of SeqKit and uses the `faidx` command to create and use FASTA index files:

- When processing UniRef90 databases, the pipeline automatically creates a `.fai` index file if one doesn't exist
- Once indexed, subsequent retrieval operations are significantly faster
- For large databases (like your 90GB UniRef90), indexed retrieval can be 10-100x faster than direct parsing

To manually pre-index your database for even faster startup:
```bash
# Create FASTA index
seqkit faidx uniref90.fasta

# Verify index was created
ls -la uniref90.fasta.fai
```

The pipeline uses `seqkit grep -f ids.txt uniref90.fasta -p -n -j 4` for efficient sequence retrieval, which:
- `-f ids.txt`: Reads IDs from a file
- `-p`: Enables pattern matching (handles variations in FASTA headers)
- `-n`: Only searches in sequence names/headers (not in sequences)
- `-j 4`: Uses 4 threads for parallel processing

If you encounter any errors with seqkit commands, the pipeline will automatically fall back to alternative retrieval methods.

## Troubleshooting

### Sequence Retrieval Issues

If you encounter problems with sequence retrieval (e.g., "Retrieved 0 sequences from UniRef90"), try the following:

1. **Check the retrieval method**: In `config.yaml`, set `retrieval_method` to one of:
   - `"direct"`: Uses Bio.SeqIO to parse the entire FASTA file directly (best for files <10GB)
   - `"chunked"`: Uses the chunked approach with multiple processes (default for larger files)
   - `"fallback"`: Uses a simple grep-based approach as a last resort
   - `"auto"`: Automatically selects the best method based on file size (default)

2. **Adjust memory settings**: For large FASTA files, adjust the `available_memory_gb` parameter in `config.yaml` to match your system's available RAM.

3. **Manual batch size**: Set `batch_size_override` in `config.yaml` to a smaller number (e.g., 500) if you're experiencing memory issues.

4. **Check file formats**: Ensure your UniRef90 FASTA file is properly formatted with '>' characters at the start of each header line.

5. **Pre-index your database**: Run `seqkit faidx uniref90.fasta` before running the pipeline to create an index file.

6. **Check for specific errors**: Examine the log file in the `logs/` directory for specific error messages.

Example configuration for difficult retrieval cases:
```yaml
parameters:
  threads: 8
  available_memory_gb: 12
  retrieval_method: "direct"  # Use direct method for problematic files
  batch_size_override: 500    # Use smaller batches if memory is limited
```

### UniRef90 Annotation Feature

The pipeline now automatically detects and annotates UniRef90 sequence identifiers:
- For sequences with `UniRef90_` prefix, annotations are fetched from the UniRef API
- The pipeline also provides a dedicated `fetch_uniref90_cluster_info()` function for retrieving detailed cluster information
- If UniRef90 annotation fails, the pipeline falls back to standard UniProt annotation

## Usage

### Basic Usage

```bash
python main.py --uniref90 /path/to/uniref90.fasta --config config.yaml
```

### Command Line Arguments

#### Required Arguments:
- `--uniref90`: Path to the UniRef90 database FASTA file

#### Configuration Options:
- `--config`: Path to the configuration file (default: config.yaml)
- `--seed_proteins`: Path to a file containing seed UniProt IDs (overrides config)
- `--output_dir`: Directory to store results (overrides config)
- `--log_level`: Logging level, one of: DEBUG, INFO, WARNING, ERROR, CRITICAL (default: INFO)

#### Pipeline Control:
- `--start_phase`: Phase to start execution from (1-5, default: 1)
- `--end_phase`: Phase to end execution at (1-5, default: 5)

#### Custom Input Options:
- `--use_custom_seed`: Use a custom seed C-A-T modules FASTA file directly for phase 2 without requiring other phase 1 outputs
- `--seed_cat_fasta`: Path to seed C-A-T modules FASTA file (required if start_phase > 1 or if using --use_custom_seed)
- `--seed_cat_msa`: Path to seed alignment file (required if start_phase > 1 and not using --use_custom_seed)
- `--cat_hmm`: Path to HMM file (required if start_phase > 1 and not using --use_custom_seed)

#### Intermediate Files (for restarting):
- `--domtblout_file`: Path to HMMER domtblout file (required if start_phase > 2)
- `--filtered_validated_fasta`: Path to filtered validated FASTA file (required if start_phase > 3)
- `--final_fasta`: Path to final FASTA file (required if start_phase > 3)
- `--annotations_tsv`: Path to annotations TSV file (required if start_phase > 4)

## Use Cases

### 1. Full Pipeline Run with Default Settings

Run the complete pipeline using the UniRef90 database and the default configuration:

```bash
python main.py --uniref90 uniref90.fasta
```

### 2. Using Custom Seed Proteins

Use a custom list of UniProt IDs as seed proteins instead of the default ones:

```bash
python main.py --uniref90 uniref90.fasta --seed_proteins my_seed_proteins.txt
```

### 3. Using Pre-extracted C-A-T Modules

If you already have a FASTA file with C-A-T modules and want to skip the UniProt extraction step:

```bash
python main.py --uniref90 uniref90.fasta --seed_cat_fasta my_cat_modules.fasta --use_custom_seed --start_phase 2
```

### 4. Resuming from a Specific Phase

Resume the pipeline from phase 3, using existing output files from phases 1 and 2:

```bash
python main.py --uniref90 uniref90.fasta \
    --start_phase 3 \
    --seed_cat_fasta results/seed_cat_modules.fasta \
    --seed_cat_msa results/seed_cat_msa.sto \
    --cat_hmm results/cat_module.hmm \
    --domtblout_file results/uniref90_hits.domtblout
```

### 5. Running Only Clustering and Visualization (Phase 5)

Run only the clustering and visualization phase using existing output files:

```bash
python main.py --uniref90 uniref90.fasta \
    --start_phase 5 \
    --seed_cat_fasta results/seed_cat_modules.fasta \
    --seed_cat_msa results/seed_cat_msa.sto \
    --cat_hmm results/cat_module.hmm \
    --domtblout_file results/uniref90_hits.domtblout \
    --filtered_validated_fasta results/filtered_validated_cat_modules.fasta \
    --final_fasta results/final_cat_modules.fasta \
    --annotations_tsv results/final_cat_modules_annotations.tsv
```

### 6. Running With Pfam-A.hmm Validation

Validate identified modules against the Pfam-A.hmm database to ensure they contain the expected domains:

```bash
python main.py --uniref90 uniref90.fasta \
    --pfam_hmm /path/to/Pfam-A.hmm \
    --config config.yaml
```

This will:
- Run the standard pipeline
- Additionally validate extracted modules against Pfam-A.hmm
- Filter out modules that don't contain legitimate C-A-T domains according to Pfam
- Provide detailed statistics on how many modules pass Pfam validation

## Pipeline Phases

1. **Setup, Seed Data Acquisition, and Seed Alignment**
   - Parse seed protein UniProt IDs
   - Fetch sequences and domain information
   - Extract C-A-T modules
   - Align modules and build HMM profile

2. **pHMM Construction and Database Search**
   - Search UniRef90 with the C-A-T module HMM profile

3. **Hit Filtering, Validation, and Sequence Retrieval**
   - Filter HMMER hits
   - Validate domain order and spacing
   - Extract valid C-A-T modules
   - Remove redundancy

4. **Annotation**
   - Annotate modules with metadata
   - Generate annotation file
   - **Enhanced UniRef90 Support**: Automatically retrieves annotations for UniRef90 identifiers using the UniRef API
   - Provides detailed cluster information including taxonomy distribution and representative members

5. **Clustering and Network Visualization**
   - Cluster modules by sequence similarity
   - Create sequence similarity network
   - Generate network visualization

## Output Files

The pipeline generates several output files:

- `results/seed_cat_modules.fasta`: Extracted C-A-T modules from seed proteins
- `results/seed_cat_msa.sto`: Multiple sequence alignment of seed modules
- `results/cat_module.hmm`: The profile HMM built from seed alignments
- `results/uniref90_hits.domtblout`: HMMER search results against UniRef90
- `results/filtered_validated_cat_modules.fasta`: Validated C-A-T modules
- `results/final_cat_modules.fasta`: Final non-redundant C-A-T module sequences
- `results/final_cat_modules_annotations.tsv`: Annotations for each module
- `results/final_cat_modules_annotated.fasta`: Annotated FASTA file
- `results/cluster_results_clu.tsv`: Clustering results
- `results/network_visualization.png`: Network visualization image
- `logs/pipeline_[timestamp].log`: Detailed execution log

## Project Structure

```
├── main.py                # Main orchestrator script
├── config.yaml            # Configuration file
├── scripts/               # Core logic modules
│   ├── setup.py           # Environment checks
│   ├── data_acquisition.py # Fetching sequences
│   ├── validation.py      # Domain validation
│   ├── alignment.py       # Sequence alignment and HMM operations
│   ├── clustering.py      # Sequence clustering and network visualization
│   └── annotation.py      # Module annotation
├── data/                  # Input data
│   └── seed_proteins.txt  # Default seed UniProt IDs
├── results/               # Output files
├── logs/                  # Log files
├── tmp/                   # Temporary files
└── tests/                 # Unit tests
```

## Troubleshooting

### Common Issues

1. **"Failed to fetch sequence for [UniProt ID]"**
   - Check your internet connection
   - Verify that the UniProt ID exists

2. **"module 'tempfile' has no attribute 'StringIO'"**
   - This is a known issue in older versions
   - Update to the latest version or modify the code to use `io.StringIO` instead

3. **Errors with UniProt queries when using FASTA sequences as input**
   - Use the `--use_custom_seed` flag to bypass UniProt queries
   - Make sure to specify `--start_phase 2` to skip phase 1

4. **UniRef90 Annotation Feature**
   - The pipeline now automatically detects and annotates UniRef90 sequence identifiers
   - For sequences with `UniRef90_` prefix, annotations are fetched from the UniRef API
   - The pipeline also provides a dedicated `fetch_uniref90_cluster_info()` function for retrieving detailed cluster information
   - If UniRef90 annotation fails, the pipeline falls back to standard UniProt annotation

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 