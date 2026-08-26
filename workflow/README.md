# Transcript Reconstruction and Quantification Workflow

This workflow performs transcript reconstruction and quantification from long-read RNA-seq data. It can also incorporate short-read data for a comparative analysis.

## 1. Configuration

Before running the workflow, you need to configure it by editing the `config/config.yml` file.

**Key configuration options:**

*   `metadata`: Path to the metadata file. This file should contain information about your samples, including their conditions.
*   `output_dir`: The directory where the results will be saved.
*   `reference_genome`: Path to the reference genome in FASTA format.
*   `reference_annotation`: Path to the reference annotation in GTF format.
*   `experiment`: A name for your experiment. It will be used as a prefix for output files.
*   `tools`: A list of transcript reconstruction and quantification tools to use. Available tools are: `kallisto`, `isoquant`, `bambu`, `oarfish`, `flair`.
*   `sr_fastq`: (Optional) Path to a manifest file for short-read FASTQ files. If provided, the workflow will also perform a long-read vs. short-read comparison.

## 2. Running the Workflow

This workflow is managed by Snakemake. Make sure you have Snakemake and the required software (see `envs` directory) installed.

**To run the workflow:**

1.  **Navigate to the workflow directory:**
    ```bash
    cd transcripts_workflow/workflow
    ```

2.  **Perform a dry run (recommended):**
    This will show you the jobs that will be executed without actually running them.
    ```bash
    snakemake -n
    ```

3.  **Execute the workflow:**
    *   **On a local machine:**
        ```bash
        snakemake --cores <number_of_cores>
        ```
    *   **On a cluster with SLURM:**
        The `Snakefile` is pre-configured to use SLURM. You can submit the jobs to the cluster with the following command:
        ```bash
        snakemake --cluster "sbatch --qos=<queue> --cpus-per-task={threads} --mem={resources.mem_mb}" --jobs <number_of_jobs>
        ```

## 3. Output

The results of the workflow will be saved in the directory specified by the `output_dir` option in the configuration file. The output directory will contain subdirectories for each tool used in the analysis.
