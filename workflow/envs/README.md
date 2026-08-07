# Conda Environments & Software Version Matrix

This directory contains the isolated **Conda environment definitions** (`*.yaml` / `*.yml`) used across the Snakemake workflow rules. Each rule executes inside its dedicated Conda environment to ensure 100% computational reproducibility.

---

## 🛠️ Key Software Tools & Version Matrix

| Tool / Package | Version | Conda Environment | Purpose |
|---|---|---|---|
| **kallisto** | `0.51.1` | `envs/kallisto.yaml` | Long-read and short-read pseudoalignment & BUS quantification |
| **Bambu** | `3.12.1` | `envs/bambu.yaml` | Long-read transcript reconstruction & quantification |
| **IsoQuant** | `3.6.0` | `envs/isoquant.yaml` | Reference-based transcript identification & quantification |
| **Oarfish** | `0.7.0` | `envs/oarfish.yaml` | Probabilistic long-read transcript quantification |
| **FLAIR** | `2.0.0` | `envs/flair3.yaml` | Full-length transcript isoform identification |
| **SQANTI3** | `5.3.0` | `envs/SQANTI3.yml` | Quality control & structural categorization of transcript models |
| **NOISeq** | `2.46.0` | `envs/NOIseq.yaml` | Normalization, count matrix construction, and differential analysis |
| **minimap2** | `2.28` | `envs/align.yaml` | Long-read spliced alignment to reference genome |
| **samtools** | `1.20` | `envs/align.yaml` | BAM indexing, sorting, and alignment processing |
| **seqkit** | `2.8.2` | `envs/seqkit.yaml` | FASTQ/FASTA sequence manipulation and statistics |
| **gffread** | `0.12.7` | `envs/gffread.yaml` | GTF/GFF conversion, filtering, and sequence extraction |
| **r-tidyverse** | `2.0.0+` | `envs/tidy.yaml` | Data wrangling, aggregation, and plot generation |

---

## 📦 Automatic Provisioning

When executing Snakemake with the `--use-conda` flag, Snakemake automatically creates and manages these environments inside `.snakemake/conda/`:

```bash
snakemake --cores 8 --use-conda --configfile config/config.yml
```
