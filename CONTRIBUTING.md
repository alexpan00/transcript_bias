# Contributing to transcripts_workflow

Thank you for your interest in contributing! This document describes how to report issues, propose changes, and submit pull requests.

---

## Table of Contents

1. [Code of Conduct](#code-of-conduct)
2. [How to Report a Bug](#how-to-report-a-bug)
3. [How to Request a Feature](#how-to-request-a-feature)
4. [Development Setup](#development-setup)
5. [Adding a New Tool](#adding-a-new-tool)
6. [Coding Style](#coding-style)
7. [Pull Request Checklist](#pull-request-checklist)

---

## Code of Conduct

Please be respectful and constructive in all interactions. We follow the [Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).

---

## How to Report a Bug

1. Check the [Issues](../../issues) page to see if it has been reported already.
2. Open a **new issue** using the **Bug report** template.
3. Include:
   - Your OS and Snakemake version (`snakemake --version`).
   - The relevant section of `config/config.yml` (with paths anonymised).
   - The full error message and log file (usually in `logs/`).
   - The exact Snakemake command you ran.

---

## How to Request a Feature

Open an issue using the **Feature request** template. Describe the tool or analysis you would like to add and, if possible, link to its publication and repository.

---

## Development Setup

```bash
git clone https://github.com/<your-org>/transcripts_workflow.git
cd transcripts_workflow
mamba create -n snakemake -c conda-forge -c bioconda snakemake pandas
conda activate snakemake
```

We recommend testing with a small dataset before opening a pull request.

---

## Adding a New Tool

The workflow is designed to be modular. To add a new transcript reconstruction or quantification tool:

1. **Create a rule file** in `workflow/rules/<tool>.smk` following the conventions in existing rule files (e.g., `bambu.smk` or `oarfish.smk`).
2. **Add a Conda environment** in `workflow/envs/<tool>.yaml` pinning the tool version.
3. **Register the tool name** in the `tools` list in `workflow/Snakefile`.
4. **Include the rule file** in `workflow/Snakefile` with `include: "rules/<tool>.smk"`.
5. **Add the tool's output** to the `rule all` input list so it is included in the default run.
6. **Update documentation**: add a row to the "Tools & Environments" table in `README.md` and document any new config options in `workflow/config/README.md`.

---

## Coding Style

### Snakemake rules

- Use descriptive rule names in `snake_case`.
- Always include `log:` and `benchmark:` directives using the `LOGS` and `BENCHMARKS` path variables.
- Always specify a `conda:` environment pointing to the correct `envs/` file.
- Use `resources: mem_mb=...` with the helper functions from `scripts/utils.py` for memory-intensive jobs.
- Prefer `shell:` blocks for external tool calls; use `run:` blocks only for simple Python logic that does not require a separate environment.

### Python scripts

- Use Python 3.10+ syntax.
- Add a module-level docstring describing the script's purpose, inputs, and outputs.
- Use `argparse` for command-line argument parsing.
- Follow PEP 8.

### R scripts

- Add a header comment describing the script's purpose, expected arguments, and outputs.
- Use `tidyverse` where appropriate for readability.
- Avoid `setwd()`; use relative paths from the script's perspective or pass paths as command-line arguments.

---

## Pull Request Checklist

Before submitting a pull request, please verify:

- [ ] The workflow runs without errors on a test dataset (`snakemake -n --use-conda`).
- [ ] New rules have `log:` and `benchmark:` directives.
- [ ] New environments are in `workflow/envs/` with pinned versions.
- [ ] `README.md` and `workflow/config/README.md` are updated if configuration options changed.
- [ ] `CHANGELOG.md` has an entry under `[Unreleased]`.
- [ ] No hardcoded absolute paths in any source file.
