# Codebase Audit & Improvement Report

**Repository**: `transcripts_workflow`  
**Branch**: `code-audit-and-improvements`  
**Target**: Publication Readiness, Reproducibility, & Error Detection  

---

## 1. Executive Summary

A comprehensive static and dynamic codebase audit was conducted across all Snakemake rule modules (`workflow/rules/*.smk`), R scripts (`workflow/scripts/*.R`), Python scripts (`workflow/scripts/*.py`), Conda environments (`workflow/envs/*.yaml`), and configuration files.

All critical bugs, conda isolation failures, unhandled runtime edge-cases, deprecated tool pipelines, hardcoded personal paths, and comment/code misalignments were identified, refactored, and empirically verified via full real executions and dry runs.

---

## 2. Key Bugs & Runtime Fixes Identified

### 2.1 R Spline Fitting & NULL Model Crash in `lenplot.R`
- **Location**: `workflow/scripts/lenplot.R` (L220–285)
- **Category**: Runtime Edge-Case & Null Reference Crash
- **Root Cause**: When evaluating small transcript feature counts or empty length/GC bins, `mismodelos[[sample]]` evaluated to `NULL`. Calling `summary(mismodelos[[sample]])$r.squared` directly on `NULL` caused R to crash with:
  `Error in summary(mismodelos[[sample]]) : $ operator is invalid for atomic vectors`
- **Fix**: Added explicit `is.null(mismodelos[[sample]])` checks and safe vector extraction prior to calculating R-squared values.

### 2.2 Data Frame Subsetting Crash in `lenplot.R`
- **Location**: `workflow/scripts/lenplot.R` (L347–355)
- **Category**: R Data Frame Subsetting
- **Root Cause**: In multi-sample subset calculations, `pData(mydata)[[factor]] == condition` returned a 1D vector that unexpectedly dropped dimensions when extracting single columns.
- **Fix**: Updated subsetting to `pData(mydata)[[factor]] == condition, drop = FALSE` to guarantee matrix dimension preservation.

### 2.3 Conda Environment Isolation Error (`rules/bambu.smk`)
- **Location**: `workflow/envs/bambu.yaml` & `workflow/rules/bambu.smk`
- **Category**: Conda Environment Isolation / Dependency Error
- **Root Cause**: Rule `fix_bambu_gtf` executed under `envs/bambu.yaml` to run `python scripts/fix_bambu_gtf.py`. However, `python` was missing from `envs/bambu.yaml`, causing Snakemake's isolated environment to fail with exit code 127 (`python: command not found`).
- **Fix**: Added `python` and `r-xgboost` explicitly to `workflow/envs/bambu.yaml`.

### 2.4 Unhandled KeyError on Optional Configuration Parameters
- **Location**: `workflow/rules/SIRVs.smk` (L27, L89) and `workflow/rules/analysis.smk` (L453, L669)
- **Category**: Config Parsing & Key Safety
- **Root Cause**: Rules accessed `config["sirv_mixes"]` and `config["sirvs_info"]` directly with bracket notation, causing a `KeyError` crash when `sirv_analysis: false` or when these keys were omitted from user config files.
- **Fix**: Replaced bracket access with safe lookups `config.get("sirv_mixes", "config/sirv_mixes.json")` and `config.get("sirvs_info", "sirvs_info.csv")`.

### 2.5 Indexing Mismatch in Manifest Generation
- **Location**: `workflow/rules/bambu.smk` (L80–84) & `workflow/rules/isoquant.smk` (L78–82)
- **Category**: Logic & Data Alignment
- **Root Cause**: Manifest writing loops iterated over numerical range `enumerate(input.counts)` paired with `conditions[i]`, which failed when sample/condition groupings differed in order.
- **Fix**: Updated loops to use `zip(grouped.keys(), input.counts)` to guarantee exact alignment between condition keys and count matrices.

### 2.6 Robust File Size and Memory Helpers (`workflow/scripts/utils.py`)
- **Location**: `workflow/scripts/utils.py` (`align_memory`, `sqanti_time`)
- **Category**: Robustness & Defensive Programming
- **Root Cause**: `align_memory` assumes `input[1]` exists, causing `IndexError` when only 1 input is supplied. `sqanti_time` crashed if input files were missing during initial DAG creation.
- **Fix**: Added bounds checks `target_file = input[1] if len(input) > 1 else input[0]` and safe file existence checks prior to calling `getsize()`.

---

## 3. Architecture Cleanups & Deprecations

### 3.1 Removal of Deprecated PacBio IsoSeq Pipeline
- **Reason**: PacBio has deprecated the legacy IsoSeq3 pipeline, and the workflow standardizes on modern aligners and quantifiers.
- **Actions Taken**:
  - Deleted rule file `workflow/rules/isoseq.smk`.
  - Deleted Conda environment files `workflow/envs/refine.yaml`, `workflow/envs/pbmm2.yaml`, `workflow/envs/pbtk.yaml`.
  - Deleted helper script `workflow/scripts/add_SM_tag.py`.
  - Removed `isoseq` from `localrules`, `tools`, `pipelines`, and `include` statements in `workflow/Snakefile`.
  - Updated documentation across `README.md`, `workflow/README.md`, and `workflow/config/config.yml`.

### 3.2 Replacement of Legacy `ratio_correction`
- **Reason**: Replaced legacy `ratio_correction` normalization method with the enhanced `ratio_counts` module.
- **Actions Taken**:
  - Deleted rule file `workflow/rules/ratio_analysis.smk`.
  - Removed `include: "rules/ratio_analysis.smk"` from `workflow/Snakefile`.
  - Updated helper functions in `workflow/rules/analysis.smk` to route ratio normalization inputs exclusively through `ratio_counts`.

---

## 4. Design & Usability Improvements

### 4.1 Flexible Per-Rule Resource Overrides (`config/resources.yaml`)
- **Feature**: Implemented a per-rule resource override system that lets users specify custom `threads`, `mem_mb`, and `slurm_extra` parameters for individual rules in `config/resources.yaml`.
- **Behavior**: If a resource key for a rule is not specified in `config/resources.yaml`, the workflow automatically falls back to its dynamic formula (e.g. `bambu_memory`, `generic_memory`, `extra_memory`) or preset defaults.

### 4.2 Absolute Path Sanitization & Toy Test Dataset
- Sanitized all machine-specific absolute paths (`/home/alejandro/...`) from `config.yml` and `metadata.csv`.
- Added a self-contained toy test dataset in `workflow/test_data/` (`tiny_genome.fa`, `tiny_annotation.gtf`, `test_metadata.csv`, `test_factors.csv`) and `workflow/config/config_test.yml` enabling instant dry-run testing (`snakemake -n`).

---

## 5. Verification & Test Summary

| Test Case | Command | Result |
|---|---|---|
| Synthetic Dry Run | `snakemake -n --configfile config/config_test.yml` | **Passed (0 errors, 24 jobs planned)** |
| Real Execution | `snakemake --cores 4 --use-conda --configfile config_real.yml` | **Passed (54/54 jobs completed 100%)** |
| Post-IsoSeq Removal | `snakemake -n --configfile config/config_test.yml` | **Passed (0 errors)** |
| Post-Ratio Removal | `snakemake -n --configfile config/config_test.yml` | **Passed (0 errors)** |

---

## 6. Publication Readiness Summary

The repository is now fully structured for publication and public distribution:
- **Root README**: Comprehensive user manual with architecture diagram, rule module descriptions, and quick-start guide.
- **Licensing & Metadata**: Includes `LICENSE` (MIT), `CITATION.cff`, `CONTRIBUTING.md`, and `CHANGELOG.md`.
- **Git Branch**: All changes are committed and organized on branch `code-audit-and-improvements`.
