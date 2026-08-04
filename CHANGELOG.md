# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added
- Standalone repository structure extracted from the `lr_bias` project.
- Comprehensive `README.md` with installation, configuration, and usage instructions.
- `CONTRIBUTING.md` with guidelines for contributors.
- `CHANGELOG.md` (this file).
- `LICENSE` (MIT).
- `.gitignore` tuned for Snakemake + Conda workflows.
- Template `config/config.yml` with placeholder paths (no hardcoded personal paths).

### Changed
- `config/config.yml` updated to replace all absolute personal paths with descriptive placeholders.

---

## [1.0.0] – Initial release

### Added
- Snakemake workflow for long-read RNA-seq transcript reconstruction and quantification.
- Support for six tools: **kallisto** (LR), **bambu**, **IsoQuant**, **IsoSeq**, **oarfish**, **FLAIR**.
- Optional short-read comparison arm via **kallisto** (SR).
- Normalization module: CPM, TPM, TMM, EDA, CQN, ratio-correction.
- SIRV and ERCC spike-in analysis rules.
- SQANTI3 structural classification integrated for all pipeline outputs.
- TAMA merge for condition-level transcript collapse.
- Unique Junction Chain (UJC) cross-tool comparison and UpSet plots.
- Dynamic memory and time resource allocation functions for SLURM.
- Per-rule Conda environments for full reproducibility.
- Snakemake report integration with captions for all key output figures.
