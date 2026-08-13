# ==============================================================================
# Transcripts Workflow - Makefile Shortcuts
# ==============================================================================

# Executable & Configuration variables (can be overridden via environment or CLI)
SNAKEMAKE   ?= snakemake
CORES       ?= 4
CONFIG      ?= config/config.yml
REPORT      ?= report.html
EXTRA_FLAGS ?=

.PHONY: help test run report clean

help:
	@echo "Transcripts Workflow Management"
	@echo ""
	@echo "Usage:"
	@echo "  make test        Run dry-run validation using config (CONFIG=$(CONFIG))"
	@echo "  make run         Execute workflow using config (CONFIG=$(CONFIG), CORES=$(CORES))"
	@echo "  make report      Generate interactive HTML report (REPORT=$(REPORT))"
	@echo "  make clean       Remove log, benchmark, and results directories"
	@echo ""
	@echo "Customization Examples:"
	@echo "  make run CONFIG=config/config_real.yml CORES=16"
	@echo "  make report REPORT=my_report.html CONFIG=config/config.yml"
	@echo "  make run SNAKEMAKE=/path/to/snakemake CORES=8"

test:
	cd workflow && $(SNAKEMAKE) -n --configfile $(CONFIG) $(EXTRA_FLAGS)

run:
	cd workflow && $(SNAKEMAKE) --cores $(CORES) --use-conda --configfile $(CONFIG) $(EXTRA_FLAGS)

report:
	cd workflow && $(SNAKEMAKE) --report $(REPORT) --configfile $(CONFIG) $(EXTRA_FLAGS)
	@echo "Report generated at workflow/$(REPORT)"

clean:
	rm -rf workflow/logs/* workflow/benchmarks/* workflow/results/ workflow/results_test/ workflow/results_real/
