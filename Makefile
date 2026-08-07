# ==============================================================================
# Transcripts Workflow - Makefile Shortcuts
# ==============================================================================

SNAKEMAKE := /home/alejandro/miniconda3/envs/snakemake/bin/snakemake
CORES     ?= 4

.PHONY: help test run report clean clean-test

help:
	@echo "Available commands:"
	@echo "  make test      - Run dry-run validation using synthetic test dataset"
	@echo "  make run       - Execute full workflow using workflow/config/config.yml"
	@echo "  make report    - Generate interactive Snakemake HTML report (workflow/report.html)"
	@echo "  make clean     - Clean log and benchmark directories"
	@echo "  make clean-test - Remove synthetic test outputs (workflow/results_test/)"

test:
	cd workflow && $(SNAKEMAKE) -n --configfile config/config_test.yml

run:
	cd workflow && $(SNAKEMAKE) --cores $(CORES) --use-conda --configfile config/config.yml

report:
	cd workflow && $(SNAKEMAKE) --report report.html --configfile config/config.yml
	@echo "Report generated at workflow/report.html"

clean:
	rm -rf workflow/logs/* workflow/benchmarks/* workflow/gurobi.log gurobi.log

clean-test:
	rm -rf workflow/results_test/
