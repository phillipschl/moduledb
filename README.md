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
```

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

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 