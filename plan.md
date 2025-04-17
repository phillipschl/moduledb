# Plan: NRPS C-A-T Module Sequence Database Generation

**Goal:** To create a high-quality sequence database of continuous C-A-T (Condensation-Adenylation-Thiolation) modules from Non-Ribosomal Peptide Synthetase (NRPS) enzymes, leveraging HMMER, UniRef90, InterPro/Pfam, and MMseqs2, with robust validation, logging, and annotation.

**Key Principles:**

* **Focus on C-A-T Continuity:** The primary goal is to identify and extract sequences representing complete, contiguous C-A-T modules.
* **Best Practices:** Employ modular code, clear naming conventions, configuration files, comprehensive logging, unit/integration testing, and assertions.
* **No Over-Engineering:** Start with functional components and essential features; avoid premature optimization or unnecessary complexity.
* **Detailed Logging & Assertions:** Implement thorough logging for traceability and assertions for data validation at critical steps.
* **Reproducibility:** Ensure the pipeline can be rerun with clear inputs and outputs.

---

## Phase 1: Setup, Seed Data Acquisition, and Seed Alignment

**Objective:** Prepare the environment and create a high-quality seed alignment of known, continuous C-A-T modules.

**Tasks:**

1.  **Environment Setup:**
    * Install necessary software:
        * HMMER suite (v3.x: `hmmbuild`, `hmmsearch`, `hmmscan` potentially for refinement)
        * MMseqs2
        * Python (3.8+) with libraries:
            * `biopython` (Sequence handling, parsing)
            * `pandas` (Data manipulation, tabular results)
            * `requests` (API calls, e.g., UniProt)
            * `logging` (Built-in)
            * `pytest` (Testing framework)
            * `pyyaml` / `configparser` (Configuration management)
            * `networkx`, `matplotlib`/`seaborn` (Optional, for direct visualization)
    * Set up project directory structure (see below).

2.  **Define Validation Function (`validate_cat_order`):**
    * Create a Python function `validate_cat_order(domain_hits, max_gap_allowed)` early on.
    * **Input:** `domain_hits`: A list of tuples/objects for a *single protein*, each representing a domain hit, containing `(domain_type, start_pos, end_pos)`. Example: `[('C', 10, 450), ('A', 460, 980), ('T', 995, 1100)]`. List *must be sorted by start position*.
    * **Input:** `max_gap_allowed`: Integer defining the maximum allowable amino acid gap between consecutive domains (e.g., 50 aa).
    * **Logic:**
        * Iterate through the sorted `domain_hits` list looking for consecutive C, A, T triplets.
        * For a potential C-A-T triplet (`C_hit`, `A_hit`, `T_hit`):
            * `assert C_hit['domain_type'] == 'C'`
            * `assert A_hit['domain_type'] == 'A'`
            * `assert T_hit['domain_type'] == 'T'`
            * Check order: `C_hit['start'] < A_hit['start'] < T_hit['start']`.
            * Check continuity (gap):
                * `A_hit['start'] - C_hit['end'] <= max_gap_allowed`
                * `T_hit['start'] - A_hit['end'] <= max_gap_allowed`
                * Ensure gaps are non-negative: `A_hit['start'] > C_hit['end']` and `T_hit['start'] > A_hit['end']`.
        * **Output:** Boolean `True` if *at least one* valid, continuous C-A-T module pattern is found in the input list for that protein, `False` otherwise. Optionally, return the coordinates `(C_start, T_end)` of the *first* valid module found.
    * **Testing:** Write unit tests (`pytest`) for this function covering various scenarios (correct order, incorrect order, gaps too large, overlapping domains, multiple modules).
    * **Logging:** Log within the function which checks are failing for debugging.

3.  **Seed Protein Selection:**
    * Identify a set of well-characterized NRPS proteins known to contain C-A-T modules (e.g., from MIBiG database, literature reviews). Use UniProtKB IDs.
    * *Manual Curation Recommended:* Prioritize proteins where the C-A-T structure is experimentally validated or strongly supported. Aim for diversity if possible.
    * *Alternative (Automated but needs verification):* Query UniProtKB/InterPro for proteins annotated with Pfam IDs PF00668 (C), PF00501 (A), and PF00550 (T). This might yield many hits needing further filtering.

4.  **Seed Domain Data Acquisition:**
    * For each selected seed protein UniProt ID:
        * Fetch the full amino acid sequence (e.g., using Biopython and UniProt API).
        * Fetch domain annotations (Pfam hits with start/end coordinates) via InterProScan results (available via UniProt API or InterPro website/API). Store these as `(protein_id, domain_type, pfam_id, start, end)`.
        * **Crucially:** Filter domain hits to *only* include PF00668 (C), PF00501 (A), PF00550 (T).

5.  **Seed C-A-T Module Identification & Extraction:**
    * For each seed protein:
        * Retrieve its C, A, T domain hits (from step 4). Sort them by start position.
        * Apply the `validate_cat_order` function to identify valid, continuous C-A-T modules.
        * Log which seed proteins contain valid modules and which don't (and why).
        * For each *validated* continuous C-A-T module identified by `validate_cat_order` (potentially returning the start/end coordinates):
            * Extract the amino acid sequence segment from the *start of the C domain* to the *end of the T domain* from the full protein sequence.
            * Create a unique, informative sequence ID (e.g., `UniProtID_CAT_Module_1_Cstart_Tend`).
            * Store these extracted C-A-T sequences in a multi-FASTA file (`seed_cat_modules.fasta`).

6.  **Seed Alignment Generation:**
    * Perform Multiple Sequence Alignment (MSA) on the `seed_cat_modules.fasta` file.
    * Use a reliable MSA tool (e.g., MAFFT, Clustal Omega).
        ```bash
        mafft --auto seed_cat_modules.fasta > seed_cat_msa.sto
        ```
    * Ensure the output format is Stockholm (`.sto`), suitable for `hmmbuild`.
    * **Assertion:** Check if the number of sequences in the output MSA matches the number of extracted seed sequences.

**Logging:** Log number of seed proteins considered, number passing validation, number of C-A-T modules extracted, paths to intermediate and final files.

---

## Phase 2: pHMM Construction and Database Search

**Objective:** Build a profile Hidden Markov Model (pHMM) from the seed alignment and search it against the UniRef90 database.

**Tasks:**

1.  **pHMM Construction:**
    * Use `hmmbuild` from HMMER.
    * Input: `seed_cat_msa.sto`
    * Output: `cat_module.hmm`
    * Command:
        ```bash
        hmmbuild cat_module.hmm seed_cat_msa.sto
        ```
    * **Assertion:** Check if the `cat_module.hmm` file was created and is not empty.

2.  **Prepare UniRef90 Database:**
    * Ensure the UniRef90 FASTA file (provided via command-line argument in the main script) is accessible.
    * If needed (depending on HMMER version/workflow), prepare the database for `hmmsearch` (e.g., `hmmpress` - though often not strictly required for FASTA input).

3.  **Perform pHMM Search:**
    * Use `hmmsearch` to search the `cat_module.hmm` against UniRef90.
    * Optimize for speed if necessary (e.g., `--cpu` option).
    * Output full sequence scores (`-o`), domain scores (`--domtblout`), and potentially per-hit HMMER output (`-A` if needed for alignment details, though often large). Tabular domain output is usually most useful.
    * Command:
        ```bash
        hmmsearch --domtblout results/uniref90_hits.domtblout \
                  -o results/uniref90_hits.log \
                  --cpu <num_cores> \
                  cat_module.hmm <path_to_uniref90.fasta>
        ```
    * **Logging:** Log the HMMER command executed, start/end times, path to UniRef90 file used.

4.  **Parse HMMER Results:**
    * Read the `--domtblout` file (`uniref90_hits.domtblout`). This contains information about domains found in target sequences matching the pHMM.
    * Use Pandas to load the tabular data, skipping comment lines. Pay attention to column headers (target name, accession, query name, E-value, score, bias, domain start/end, etc.).
    * **Assertion:** Check if the output file exists and can be parsed. Check if expected columns are present.

**Logging:** Log the number of raw domain hits parsed from the `domtblout` file.

---

## Phase 3: Hit Filtering, Validation, and Sequence Retrieval

**Objective:** Filter the raw HMMER hits based on quality, validate the C-A-T structure and continuity, and retrieve the corresponding continuous amino acid sequences.

**Tasks:**

1.  **E-value Filtering:**
    * Filter the parsed domain hits based on a stringent domain-specific E-value threshold (e.g., `domain_i-evalue <= 1e-10`). Make this configurable.
    * Keep track of hits passing this filter.

2.  **Group Hits by Target Protein:**
    * Group the filtered domain hits by the target protein ID (UniRef90 ID).

3.  **Identify Potential C-A-T Modules per Protein:**
    * For *each* target protein ID with hits:
        * **Re-evaluate Domain Boundaries (Refinement - Optional but Recommended):** The initial `hmmsearch` used a pHMM built from *full* C-A-T modules. The domain boundaries reported might span the whole module. To get *precise* C, A, T boundaries within the hit protein, it's often better to run `hmmscan` with individual Pfam HMMs (PF00668, PF00501, PF00550) against the *full sequence* of the hit protein (requires fetching the full sequence first based on the UniRef90 ID).
        * **If skipping refinement:** Use the domain boundaries reported by `hmmsearch` (less precise, might treat the whole hit as one domain). This simplifies the process but might compromise accuracy of the C-A-T check. Assume for now we proceed *without* refinement, using the initial `hmmsearch` results directly *if* it provides sufficient per-domain info (which `domtblout` *does*). The `domtblout` gives coordinates *within the target sequence* for the match to the *entire* pHMM. This isn't ideal for validating *internal* C-A-T structure directly from `domtblout`.
        * ***Correction/Better Approach:*** The `hmmsearch` against the *combined* C-A-T HMM will report hits spanning the *entire* modeled region. We need to know the *internal* C, A, T locations within that hit sequence.
            * **Method A (Requires Full Sequences):**
                1. Filter `hmmsearch` hits by E-value.
                2. For each promising UniRef90 ID hit, fetch its full sequence.
                3. Run `hmmscan` using the *individual* Pfam HMMs (PF00668, PF00501, PF00550) against this *single* full sequence.
                4. Parse the `hmmscan` output to get precise C, A, T domain locations *within that specific protein*.
            * **Method B (Approximation - Less Robust):** Assume the relative positions within the seed alignment roughly translate to the hit. This is less reliable.
        * ***Adopt Method A (Refined Scan):***
            1.  Fetch full sequences for UniRef90 IDs passing E-value filter (use UniProt API or local copy, mapping UniRef90 -> UniProtKB representative). Handle potential fetching errors.
            2.  For each fetched sequence, run `hmmscan` with C, A, T Pfam HMMs:
                ```bash
                hmmscan --domtblout temp_protein_scan.domtblout \
                        --cpu 1 \
                        <pfam_c_a_t_hmm_db> \ # Need HMMs for PF00668, PF00501, PF00550
                        temp_protein.fasta
                ```
            3.  Parse `temp_protein_scan.domtblout` to get domain hits `(domain_type, start, end)` for *this specific protein*. Filter these internal domain hits by E-value again if desired.

4.  **Validate C-A-T Order and Continuity:**
    * For each protein processed in the previous step:
        * Gather its filtered C, A, T domain hits from the `hmmscan` result.
        * Sort the domain hits by start position.
        * Apply the `validate_cat_order(domain_hits, max_gap_allowed)` function.
        * Keep track of protein IDs that contain at least one validated, continuous C-A-T module. Record the start/end coordinates of the validated module(s) (e.g., `C_domain_start` to `T_domain_end`).

5.  **Retrieve Continuous C-A-T Module Sequences:**
    * For each validated protein and the corresponding module coordinates (`C_start`, `T_end`):
        * Use the previously fetched full protein sequence.
        * Extract the subsequence from `C_start` to `T_end`.
        * Assign a unique ID (e.g., `UniProtID_CAT_Module_coords_Cstart_Tend`).
        * Store these sequences in a new FASTA file (`filtered_validated_cat_modules.fasta`).

6.  **Final Redundancy Removal:**
    * Use MMseqs2 `easy-linclust` or `cd-hit` on `filtered_validated_cat_modules.fasta` to remove highly similar/identical *extracted C-A-T module sequences*. Choose an appropriate identity threshold (e.g., 0.99 or 1.0 for identical).
    * Output: `final_cat_modules.fasta`.
    * **Assertion:** Ensure the number of sequences in the final file is less than or equal to the number in the input.

**Logging:** Log number of hits after E-value filtering, number of proteins processed for domain scanning, number of proteins failing fetch, number of proteins passing C-A-T validation, number of sequences extracted, number of sequences after redundancy removal. Log any errors during sequence fetching or `hmmscan`.

---

## Phase 4: Annotation

**Objective:** Annotate the final set of C-A-T module sequences with relevant biological information.

**Tasks:**

1.  **Identify Information Sources:**
    * UniProtKB (Source protein): Organism, Taxonomy (full lineage), Protein Name, Gene Name (if available), Evidence Levels.
    * Pfam/InterPro: Specific domain annotations within the module (can re-run `hmmscan` or use initial results carefully).
    * (Optional) A-domain Specificity Prediction: Use external tools/models (e.g., NRPSpredictor2/3) if feasible to predict the activated amino acid.
    * (Optional) MIBiG: Link back to biosynthetic gene clusters if the source UniProt ID is present in MIBiG.

2.  **Gather Annotations:**
    * For each sequence ID in `final_cat_modules.fasta`:
        * Parse the ID to get the source UniProt ID and coordinates.
        * Use the UniProt ID to query UniProt API (or local database) for desired metadata (organism, taxonomy, etc.). Use batch queries if possible.
        * Store the coordinates (`C_start`, `A_start`, `A_end`, `T_end` derived from validation step).
        * (Optional) Run specificity prediction on the A-domain segment.
        * (Optional) Check MIBiG for the UniProt ID.

3.  **Store Annotated Data:**
    * Choose a suitable format:
        * **Option A (Simple):** Enhanced FASTA headers (can become very long and hard to parse).
            ```fasta
            >UniqueID|UniProtID=P12345|Organism=Bacillus_subtilis|Taxonomy=Bacteria;Firmicutes;...|Coords=10-1100|C=10-450|A=460-980|T=995-1100|PredictedAA=Val
            MSD...
            ```
        * **Option B (Recommended):** A separate annotation table (TSV/CSV file) linking the `UniqueID` from the FASTA file to the annotations. This is much cleaner for analysis.
            ```tsv
            UniqueID	SequenceFile	SourceUniProtID	Organism	Taxonomy	ModuleCoords	C_Coords	A_Coords	T_Coords	PredictedAA	MIBiG_Link
            P12345_CAT_1_10_1100	final_cat_modules.fasta	P12345	Bacillus subtilis	Bacteria;...	10-1100	10-450	460-980	995-1100	Val	BGC0000XXX
            ...
            ```
        * **Option C (Advanced):** SQLite database. Good for complex queries but adds dependency.
    * Create the annotation file (`final_cat_modules_annotations.tsv`).

**Logging:** Log progress of annotation gathering, number of sequences successfully annotated, any API errors or missing data.

---

## Phase 5: Clustering and Network Visualization

**Objective:** Cluster the final C-A-T sequences based on similarity and visualize the relationships as a network.

**Tasks:**

1.  **Sequence Clustering:**
    * Use MMseqs2 `easy-cluster` for flexible clustering.
    * Input: `final_cat_modules.fasta`
    * Parameters: Choose appropriate sequence identity (`--min-seq-id`) and coverage (`-c`, `--cov-mode`) thresholds based on desired granularity (e.g., 0.5 for broader families, 0.8 for tighter groups).
    * Output: Cluster membership file (`cluster_results_clu.tsv`).
    * Command:
        ```bash
        mmseqs easy-cluster final_cat_modules.fasta cluster_results tmp_mmseqs \
               --min-seq-id 0.7 -c 0.8 --cov-mode 0
        ```
    * **Assertion:** Check if the cluster output file is created.

2.  **Prepare Network Data:**
    * MMseqs2 can create sequence similarity networks directly, or its intermediate results (like the all-vs-all alignment) can be used.
    * Alternatively, parse the `cluster_results_clu.tsv` file to understand cluster membership.
    * For visualization, often an edge list is needed: `SeqID1 SeqID2 SimilarityScore`. This might require running an all-vs-all search (`mmseqs easy-search`) if not directly generated by the clustering workflow chosen.
    * If using `easy-cluster`, the TSV output links representative sequences to cluster members. A network could show representatives linked to members, or an all-vs-all similarity network could be built.
    * Let's assume we want an all-vs-all network:
        ```bash
        # Create MMseqs DB
        mmseqs createdb final_cat_modules.fasta final_db
        # All-vs-all search
        mmseqs search final_db final_db search_results tmp_search -s 7.0
        # Convert results to format for network (e.g., BLAST tab)
        mmseqs convertalis final_db final_db search_results network_edges.m8
        ```

3.  **Network Visualization:**
    * Parse the network edge list (`network_edges.m8` or similar). Columns typically include query, target, identity, alignment length, E-value, bit score. Filter edges based on identity or score.
    * Parse the annotation file (`final_cat_modules_annotations.tsv`).
    * Use a visualization library/tool:
        * **Python:** `NetworkX` to build the graph object. Add node attributes from the annotation file (organism, predicted AA, etc.). Use `matplotlib`, `seaborn`, or specialized libs like `pyvis` for interactive plots.
        * **External Tools:** Export the graph (edge list, node attributes) to formats compatible with Gephi or Cytoscape for advanced interactive visualization. (`. GML`, `.graphml`, simple TSV edge/node lists).
    * Customize visualization: Color nodes by taxonomy or predicted substrate, adjust layout algorithms (e.g., Fruchterman-Reingold), scale node size by sequence length or cluster size.

**Logging:** Log MMseqs2 commands, number of clusters generated, number of nodes and edges in the network, path to network files/images.

---

## Project Structure

├── main.py                # Main orchestrator script
├── config.yaml            # Configuration file (paths, thresholds, params)
├── scripts/               # Core logic modules
│   ├── init.py
│   ├── setup.py           # Environment checks, directory creation
│   ├── data_acquisition.py # Fetching sequences, domain info
│   ├── alignment.py       # Running MSA tools
│   ├── hmm.py             # Running HMMER commands
│   ├── validation.py      # Contains validate_cat_order, parsing, filtering logic
│   ├── annotation.py      # Annotation gathering logic
│   ├── clustering.py      # Running MMseqs2 commands
│   └── visualization.py   # Network generation/plotting (optional)
├── data/                  # Input data (put large files like UniRef90 outside repo or use paths from config)
│   ├── seed_proteins.txt  # List of UniProt IDs for seed
│   └── pfam_hmms/         # Directory for individual Pfam HMMs (C, A, T)
├── results/               # Output files generated by the pipeline
│   ├── seed_cat_modules.fasta
│   ├── seed_cat_msa.sto
│   ├── cat_module.hmm
│   ├── uniref90_hits.domtblout
│   ├── filtered_validated_cat_modules.fasta
│   ├── final_cat_modules.fasta
│   ├── final_cat_modules_annotations.tsv
│   ├── cluster_results_clu.tsv
│   └── network_visualization.png # Or other network formats
├── logs/                  # Log files
│   └── pipeline.log
├── tests/                 # Unit and integration tests
│   ├── init.py
│   ├── test_validation.py
│   ├── test_parsing.py
│   └── ...                # Other test files
├── tmp/                   # Temporary files (can be deleted after run)
└── README.md              # Project description, usage instructions
└── plan.md                # This file

---

## Main Script (`main.py`) Outline

* Parse command-line arguments (`argparse`): path to UniRef90, path to config file.
* Load configuration (`config.yaml`).
* Set up logging (`logging` module).
* Execute pipeline phases sequentially, calling functions from `scripts/` modules.
* Include `try...except` blocks for error handling and logging.
* Use assertions (`assert`) at critical checkpoints to validate intermediate data.
* Log progress, timings, and key statistics for each phase.

---

## Key Considerations & Potential Challenges

* **UniRef90 Size:** Searching UniRef90 is computationally intensive. Ensure sufficient resources (CPU, RAM, time).
* **Domain Boundary Precision:** Relying solely on the initial `hmmsearch` of the combined C-A-T HMM for internal boundaries is inaccurate. Re-scanning hits with individual C, A, T Pfam HMMs (`hmmscan`) is strongly recommended for accurate validation.
* **Defining "Continuity":** The `max_gap_allowed` parameter in `validate_cat_order` is crucial and might need empirical tuning. What's a biologically reasonable linker length?
* **API Rate Limits/Errors:** When fetching data from UniProt/InterPro, handle potential rate limits, network errors, and missing entries gracefully. Implement retries or batching.
* **Annotation Completeness:** Not all proteins will have rich annotations available.
* **Scalability:** For very large datasets or more complex analyses, consider database backends (SQLite, PostgreSQL) instead of flat files.
* **A-Domain Specificity:** Integrating prediction tools adds complexity and dependencies. Start without it if necessary.

This plan provides a detailed roadmap. The AI agent should focus on implementing each phase modularly, ensuring robust validation and logging throughout the process.
