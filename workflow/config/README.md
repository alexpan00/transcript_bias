The workflow is configured via the `config/config.yml` file. The configuration is organized into the following sections:

### General Settings
*   `experiment`: Name of the experiment (e.g., "isoseq").
*   `data_type`: Type of long-read data (e.g., "pacbio", "nanopore").
*   `output_dir`: Directory where results will be saved.

### Input Data
*   `metadata`: Path to the CSV file containing sample information (paths to FASTQ/BAM files, sample names, conditions). Should always include a sample and condition columns.
*   `factors`: Path to the CSV file defining experimental factors for analysis.
*   `sr_analysis`: Whether to carry or not the short-reads analysis.
*   `sr_fastq`: (Optional) Path to a TSV manifest file for short-read data comparison.

### Reference Files
*   `reference_genome`: Path to the genome FASTA file.
*   `reference_annotation`: Path to the reference GTF annotation file.
*   `cage`: (Optional) Path to CAGE peaks BED file for 5' end validation.
*   `tts`: (Optional) Path to PolyA sites BED file for 3' end validation.

### Tools & Analysis
*   `tools`: List of tools to execute. Options include: `kallisto`, `bambu`, `flair`, `isoquant`, `oarfish`.
*   `sqanti_dir`: Path to the SQANTI3 installation directory.
*   `normalization_methods`: A list with the normalization methods to run. Options
include: CPM, TPM, EDA, cqn
*   `min_count_condition`: default 0. Filter out transcripts that do not have at least min_count_condition counts in all the samples of at least one condition 

### SIRV Controls
*   `sirv_analysis`: Whether to carry or not the SIRVs analysis.
*   `sirvs_info`: Path to CSV file with SIRV spike-in information.
*   `sirv_mixes`: Path to JSON file defining SIRV mix ratios.

### Other analysis (true/false)
*   `iso_per_gene_analysis`: perform isoforms per gene analysis
*   `replicability_analysis`: perform replicability analysis
*   `coverage_analysis`: perform coverage analysis
 
### Reporting
*   `report_factors`: List of metadata columns to use for generating reports and plots (e.g., "Tissue", "Pool").

### Resource Overrides (`config/resources.yaml`)
You can optionally specify per-rule thread, memory (`mem_mb`), and SLURM (`slurm_extra`) resource overrides in `config/resources.yaml`. If a rule or resource key is omitted, the workflow automatically falls back to its default dynamic calculation or preset resource values.

