# lrRNA-seq Transcript Reconstruction & Quantification Workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A **Snakemake** workflow for transcript reconstruction and quantification from **long-read RNA-seq** data (PacBio or Oxford Nanopore). It benchmarks multiple state-of-the-art tools (kallisto, bambu, IsoQuant, oarfish, FLAIR) and provides a comprehensive downstream analysis suite including normalization, SIRV/ERCC spike-in validation, SQANTI3 quality control, and cross-tool comparison via Unique Junction Chains (UJC).

---

## Table of Contents

1. [Features](#features)
2. [Workflow Overview](#workflow-overview)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Quick Start](#quick-start)
6. [Configuration](#configuration)
7. [Input Data](#input-data)
8. [Running the Workflow](#running-the-workflow)
9. [Output](#output)
10. [Tools & Environments](#tools--environments)
11. [Reproducibility Notes](#reproducibility-notes)
12. [Citation](#citation)
13. [Contributing](#contributing)
14. [License](#license)

---

## Features

- **Multi-tool benchmarking**: Run one or more long-read transcript tools and directly compare their results.
- **Optional short-read comparison**: Integrate a short-read (kallisto) arm to assess concordance with standard RNA-seq.
- **Spike-in validation**: Built-in SIRV and ERCC analysis modules for assessing quantification accuracy.
- **Flexible normalization**: CPM, TPM, TMM, EDA, CQN, and ratio-correction normalization methods.
- **SQANTI3 QC**: Automatic transcript structural classification for all pipeline outputs.
- **UJC cross-tool comparison**: Unique Junction Chain analysis to evaluate isoform-level concordance.
- **Cluster-ready**: Resource allocation functions and SLURM parameters are built in.
- **Reproducible environments**: Each rule uses a dedicated Conda environment specified under `envs/`.

---

## Workflow Overview

```
Long-read input (FASTQ / BAM / FLNC)
        │
        ▼
  ┌─────────────────────────────────────────────────────────────┐
  │                     Preprocessing                           │
  │  • BAM indexing    • TAMA install    • Transcript-to-gene   │
  └──────────────────────────┬──────────────────────────────────┘
                             │
          ┌──────────────────┼──────────────────────┐
          ▼                  ▼                       ▼
   ┌─────────────┐   ┌──────────────┐      ┌───────────────┐
   │  Alignment  │   │ Quantification│      │ Reconstruction│
   │  (oarfish)  │   │  (kallisto)   │      │(isoseq/bambu/ │
   └──────┬──────┘   └──────┬───────┘      │ isoquant/flair)│
          │                 │              └───────┬────────┘
          └─────────────────┼────────────────── ──┘
                            │
                            ▼
              ┌─────────────────────────┐
              │  TAMA merge & SQANTI3   │
              │  structural annotation  │
              └────────────┬────────────┘
                           │
          ┌────────────────┼─────────────────┐
          ▼                ▼                  ▼
  ┌──────────────┐ ┌─────────────┐  ┌────────────────┐
  │ Normalization│ │  SIRV/ERCC  │  │ Short-read arm │
  │ (6 methods)  │ │  validation │  │  (optional)    │
  └──────┬───────┘ └──────┬──────┘  └───────┬────────┘
         └────────────────┼─────────────────┘
                          ▼
              ┌────────────────────────┐
              │   Downstream Analysis  │
              │ • Correlation/PCA      │
              │ • Length/GC bias       │
              │ • Replicability        │
              │ • UJC cross-tool       │
              │ • UpSet plots          │
              └────────────────────────┘
```

### Rule modules

| Module | File | Description |
|---|---|---|
| Preprocessing | `rules/preprocessing.smk` | BAM indexing, SQANTI3 and TAMA install, t2g mapping |
| Bambu | `rules/bambu.smk` | Reconstruction + quantification with Bambu |
| IsoQuant | `rules/isoquant.smk` | Reconstruction + quantification with IsoQuant |
| Kallisto (LR) | `rules/kallisto_lr.smk` | Bus → sort → count → quant-tcc |
| Kallisto (SR) | `rules/kallisto_sr.smk` | Short-read quantification arm |
| Oarfish | `rules/oarfish.smk` | Transcriptome alignment + EM quantification |
| FLAIR | `rules/flair.smk` | Align → correct → collapse → quantify |
| TAMA | `rules/tama.smk` | Condition-level merge with TAMA |
| SIRVs | `rules/SIRVs.smk` | SIRV and ERCC spike-in analysis |
| Analysis | `rules/analysis.smk` | NOISeq objects, normalization, plots |
| Ratio counts | `rules/ratio_counts.smk` | Ratio count generation |
| UJC | `rules/UJC.smk` | Unique Junction Chain cross-tool comparison |

---

## Requirements

- **Snakemake** ≥ 8.0 (with the `snakemake-executor-plugin-slurm` plugin for cluster runs)
- **Conda** / **Mamba** (all per-rule software is managed via dedicated Conda environments)
- **Python** ≥ 3.10 (for the Snakemake environment)
- **pandas** (for metadata parsing in the Snakemake main process)

All other software (R packages, aligners, quantifiers, etc.) is installed automatically via the `envs/` YAML files when you run Snakemake with `--use-conda`.

---

## Installation

```bash
# 1. Clone this repository
git clone https://github.com/<your-org>/transcripts_workflow.git
cd transcripts_workflow

# 2. Create the base Snakemake environment (recommended: mamba for speed)
mamba create -n snakemake -c conda-forge -c bioconda snakemake pandas
conda activate snakemake
```

> **SQANTI3** is cloned automatically by the `prepare_sqanti` rule. You do **not** need to install it manually.
> **TAMA** is cloned automatically by the `prepare_tama` rule.

---

## Quick Start

```bash
cd workflow

# 1. Edit the configuration file
nano config/config.yml

# 2. Edit the metadata file to point to your data
nano metadata.csv

# 3. Dry-run to check the workflow
snakemake -n --use-conda

# 4. Run locally with 8 cores
snakemake --cores 8 --use-conda
```

---

## Configuration

All settings are controlled via `workflow/config/config.yml`. See [`workflow/config/README.md`](workflow/config/README.md) for a full description of every option.

### Minimal required settings

```yaml
experiment: "my_experiment"    # A short label for your experiment
data_type: "pacbio"            # "pacbio" or "nanopore"
output_dir: "results"          # Output directory (relative to workflow/)

metadata: "metadata.csv"       # See "Input Data" section below
factors:  "factors.csv"        # Sample-level metadata for plotting

reference_genome:     "/path/to/genome.fa"
reference_annotation: "/path/to/annotation.gtf"

tools: ["kallisto", "bambu"]   # Subset of: kallisto bambu flair isoquant oarfish

normalization_methods: ["TMM", "CPM", "TPM"]
```

### Optional analyses (set to `true` to enable)

| Option | Description |
|---|---|
| `sirv_analysis` | Quantify and validate SIRV spike-ins |
| `sr_analysis` | Run short-read kallisto arm and compare |
| `iso_per_gene_analysis` | Isoforms-per-gene statistics |
| `replicability_analysis` | Cross-replicate reproducibility plots |
| `coverage_analysis` | Per-sample read coverage plots |
| `mixture_definition` | Validate artificial mixture proportions |

---

## Input Data

### `metadata.csv`

A comma-separated file with **one row per sample**. Required columns:

| Column | Description |
|---|---|
| `fastq` | Absolute path to the FASTQ file |
| `aligned` | Absolute path to the genome-aligned BAM file |
| `flnc` | Absolute path to the FLNC BAM (can be the same as `aligned` for non-IsoSeq tools) |
| `sample` | Unique sample identifier |
| `condition` | Experimental condition label (groups samples for joint reconstruction) |
| `SIRV` | SIRV mix applied to this sample (`E0`, `E1`, or `E2`); set to `E0` if no SIRVs were used |

Example:

```
fastq,aligned,flnc,sample,condition,SIRV
/data/Brain_1.fastq,/data/Brain_1.bam,/data/Brain_1.bam,B31,brain,E0
/data/Brain_2.fastq,/data/Brain_2.bam,/data/Brain_2.bam,B32,brain,E0
/data/Kidney_1.fastq,/data/Kidney_1.bam,/data/Kidney_1.bam,K31,kidney,E1
```

### `factors.csv`

A comma-separated file with **one row per sample** listing any additional metadata variables to stratify plots by (e.g., tissue, pool, batch). Example:

```
sample,Tissue,Pool,SIRV
B31,Brain,pool1,E0
B32,Brain,pool1,E0
K31,Kidney,pool1,E1
```

### Optional input files

| Config key | Format | Description |
|---|---|---|
| `sr_fastq` | TSV manifest | Short-read FASTQ manifest (see below) |
| `cage` | BED | CAGE peak positions for 5′-end validation |
| `tts` | BED | PolyA sites for 3′-end validation |
| `sirvs_info` | CSV | SIRV transcript lengths and expected ratios per mix |
| `sirv_mixes` | JSON | Per-mix SIRV concentration ratios |
| `mixture_definition` | CSV | Artificial mixture composition (e.g., 80% Brain + 20% Kidney) |
| `ERCC_counts` | TSV | Pre-computed ERCC count matrix |

#### Short-read manifest (`sr_fastq`)

A tab-separated file with two columns: `sample` and `path` (no header).

```
B31    /data/sr/Brain_1_R1.fastq.gz,/data/sr/Brain_1_R2.fastq.gz
B32    /data/sr/Brain_2_R1.fastq.gz,/data/sr/Brain_2_R2.fastq.gz
```

---

## Running the Workflow

The workflow must be launched from the `workflow/` directory so that relative paths resolve correctly.

```bash
cd workflow
```

### Local machine

```bash
snakemake --cores <N> --use-conda
```

### SLURM cluster

```bash
snakemake \
  --executor slurm \
  --jobs 200 \
  --use-conda \
  --default-resources slurm_account=<account> slurm_partition=<partition> \
  --latency-wait 60
```

> Alternatively, use the classic `--cluster` interface:
> ```bash
> snakemake \
>   --cluster "sbatch --qos={resources.slurm_extra} --cpus-per-task={threads} --mem={resources.mem_mb}" \
>   --jobs 200 \
>   --use-conda
> ```

### Useful flags

| Flag | Description |
|---|---|
| `-n` / `--dry-run` | Preview jobs without running them |
| `--dag \| dot -Tpng > dag.png` | Generate a directed acyclic graph of jobs |
| `--report report.html` | Generate an interactive HTML report of results |
| `--rerun-incomplete` | Re-run any incomplete jobs from a previous run |
| `--keep-going` | Continue with independent jobs if one fails |

---

## Output

All results are written to `output_dir` (default `results/`). The directory structure is:

```
results/
├── transcript_to_gene.tsv          # Reference t2g mapping
├── transcripts.tsv                 # Transcriptome stats (length, GC)
├── metadata_extended.tsv           # Metadata with file basenames
│
├── kallisto/                       # Kallisto (long-read) results
│   ├── 01_index/
│   ├── 02_bus/
│   ├── 03_bus_count/
│   ├── 04_tcc/
│   └── NOIseq/
│       ├── <exp>_counts.tsv
│       ├── <exp>_NOIseq.rds
│       └── <normalization>/
│           ├── <exp>_heatmap.png
│           ├── <exp>_pca.png
│           ├── <exp>_length_<factor>.png
│           └── …
│
├── bambu/                          # Bambu results (same sub-structure as kallisto/)
├── isoquant/                       # IsoQuant results
├── oarfish/                        # Oarfish results
├── flair/                          # FLAIR results
│
├── SIRVs/                          # SIRV/ERCC spike-in analysis
│
├── <exp>_sqanti_summary.csv        # Per-tool SQANTI3 summary
├── <exp>_sqanti_summary_global.csv
├── <exp>_normalization_length_sum.png
├── <exp>_normalization_sirvs_correlation.png
├── <exp>_UJC_upset.png             # Cross-tool UJC upset plot
└── <exp>_UJC_summary.tsv
```

---

## Tools & Environments

| Environment file | Key software |
|---|---|
| `envs/align.yaml` | samtools, minimap2 |
| `envs/bambu.yaml` | R/bambu, Python |
| `envs/combined_analysis.yaml` | R packages for cross-tool analysis |
| `envs/coverage.yaml` | deeptools |
| `envs/cqn.yaml` | R/cqn |
| `envs/flair3.yaml` | FLAIR |
| `envs/gffread.yaml` | gffread |
| `envs/isoquant.yaml` | IsoQuant |
| `envs/kallisto.yaml` | kallisto, bustools, kb-python |
| `envs/kallisto_counts.yaml` | Python/scipy for sparse matrix parsing |
| `envs/NOIseq.yaml` | R/NOISeq, ggplot2, pheatmap, tidyverse |
| `envs/oarfish.yaml` | oarfish |
| `envs/requant.yaml` | featureCounts / requantification tools |
| `envs/seqkit.yaml` | seqkit |
| `envs/SQANTI3.yml` | SQANTI3 and all its dependencies |
| `envs/tama.yaml` | TAMA |
| `envs/tidy.yaml` | R/tidyverse |

---

## Reproducibility Notes

- All per-rule software is pinned via Conda environments in `envs/`. Use `--use-conda` to ensure exact software versions are used.
- SQANTI3 is pinned to a specific git fork (`alexpan00/SQANTI3`). If you need a different version, change the URL in `rules/preprocessing.smk`.
- **Absolute paths in `metadata.csv` and `config.yml`** must be updated to match your local file system. A template `config/config.yml` with placeholder paths is provided.
- Results are fully reproducible given the same input data, config, and Conda environment lock files.

---

## Citation

If you use this workflow in your research, please cite:

> *[Your paper citation here]*

Please also cite the individual tools you use in your analysis. Key references:

- **kallisto** / **bustools**: Melsted et al., *Nature Biotechnology* (2021)
- **bambu**: Chen et al., *Nature Methods* (2023)
- **IsoQuant**: Prjibelski et al., *Nature Biotechnology* (2023)
- **oarfish**: Zare Jousheghani & Patro, *Genome Biology* (2024)
- **FLAIR**: Tang et al., *Nature Communications* (2020)
- **SQANTI3**: Tardaguila et al., *Genome Research* (2018); Pardo-Palacios et al. (ongoing)
- **TAMA**: Kuo et al., *Nature Methods* (2020)
- **NOISeq**: Tarazona et al., *Nucleic Acids Research* (2015)

---

## Contributing

Contributions are welcome! Please read [CONTRIBUTING.md](CONTRIBUTING.md) before opening a pull request.

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.