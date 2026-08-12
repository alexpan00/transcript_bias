# ==============================================================================
# Transcripts Workflow - Makefile Shortcuts
# ==============================================================================

# Executable & Configuration variables (can be overridden via environment or CLI)
SNAKEMAKE   ?= snakemake
CORES       ?= 4
CONFIG      ?= config/config.yml
TEST_CONFIG ?= config/config_test.yml
REPORT      ?= report.html
EXTRA_FLAGS ?=

.PHONY: help test run report clean clean-test

help:
	@echo "Transcripts Workflow Management"
	@echo ""
	@echo "Usage:"
	@echo "  make test        Run dry-run validation using test config (TEST_CONFIG=$(TEST_CONFIG))"
	@echo "  make run         Execute workflow using config (CONFIG=$(CONFIG), CORES=$(CORES))"
	@echo "  make report      Generate interactive HTML report (REPORT=$(REPORT))"
	@echo "  make clean       Remove log and benchmark directories"
	@echo "  make clean-test  Remove test output directory (workflow/results_test/)"
	@echo ""
	@echo "Customization Examples:"
	@echo "  make run CONFIG=config/my_config.yml CORES=16"
	@echo "  make report REPORT=my_report.html CONFIG=config/my_config.yml"
	@echo "  make run SNAKEMAKE=/path/to/snakemake CORES=8"

test:
	cd workflow && $(SNAKEMAKE) -n --configfile $(TEST_CONFIG) $(EXTRA_FLAGS)

run:
	cd workflow && $(SNAKEMAKE) --cores $(CORES) --use-conda --configfile $(CONFIG) $(EXTRA_FLAGS)

report:
	cd workflow && $(SNAKEMAKE) --report $(REPORT) --configfile $(CONFIG) $(EXTRA_FLAGS)
	@echo "Report generated at workflow/$(REPORT)"

clean:
	rm -rf workflow/logs/* workflow/benchmarks/* workflow/gurobi.log gurobi.log

clean-test:
	rm -rf workflow/results_test/
