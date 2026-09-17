# TCGA omics framework - common entry points
# Usage: make scan PROJECT=TCGA-BRCA ; make run CONFIG=config/brca_tumor_vs_normal.yml

CONFIG  ?= config/brca_tumor_vs_normal.yml
PROJECT ?= TCGA-BRCA
RSCRIPT ?= Rscript

.PHONY: help conda-env conda-env-update deps deps-full projects scan fetch run rnaseq clean-results clean-cache lint docker

help:
	@echo "make conda-env    create the 'tcga-omics' conda environment"
	@echo "make deps         install R dependencies into the current R (no conda)"
	@echo "make deps-full    also install enrichment extras"
	@echo "make projects     list all TCGA projects"
	@echo "make scan         inventory one project        (PROJECT=TCGA-BRCA)"
	@echo "make fetch        download + cache the layers  (CONFIG=...)"
	@echo "make run          full pipeline                (CONFIG=...)"
	@echo "make rnaseq       DESeq2 step only             (CONFIG=...)"
	@echo "make clean-results  delete results/ (cache is kept)"
	@echo "make lint         static check of the R sources"

conda-env:
	@command -v mamba >/dev/null 2>&1 \
	  && mamba env create -f environment.yml \
	  || conda env create -f environment.yml
	@echo ""
	@echo "Now run:  conda activate tcga-omics"

conda-env-update:
	@command -v mamba >/dev/null 2>&1 \
	  && mamba env update -f environment.yml --prune \
	  || conda env update -f environment.yml --prune

deps:
	$(RSCRIPT) install/install_dependencies.R

deps-full:
	$(RSCRIPT) install/install_dependencies.R --optional

projects:
	$(RSCRIPT) scripts/00_explore_projects.R --list

scan:
	$(RSCRIPT) scripts/00_explore_projects.R --project $(PROJECT)

fetch:
	$(RSCRIPT) scripts/01_fetch_data.R --config $(CONFIG)

run:
	$(RSCRIPT) scripts/02_run_analysis.R --config $(CONFIG)

rnaseq:
	$(RSCRIPT) scripts/02_run_analysis.R --config $(CONFIG) --only rnaseq

lint:
	$(RSCRIPT) -e 'for (f in list.files(c("R","scripts","install"), pattern="[.]R$$", full.names=TRUE)) { p <- try(parse(f), silent=TRUE); if (inherits(p,"try-error")) stop(f, ": ", p) else cat("ok  ", f, "\n") }'

clean-results:
	rm -rf results

clean-cache:
	@echo "This deletes every downloaded file. Ctrl-C to abort."; sleep 5; rm -rf cache

docker:
	docker build -t tcga-omics-framework .
