# Reproducible environment for the TCGA omics framework.
# Build:  docker build -t tcga-omics .
# Run:    docker run --rm -v "$PWD":/work -w /work tcga-omics \
#             Rscript scripts/02_run_analysis.R --config config/brca_tumor_vs_normal.yml
#
# Mount a host directory over /work so cache/ and results/ survive the container.

FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.title="TCGA omics framework" \
      org.opencontainers.image.description="Config-driven multi-omics differential analysis on TCGA" \
      org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install -y --no-install-recommends \
        libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/tcga-omics
COPY install/install_dependencies.R install/
RUN Rscript install/install_dependencies.R --optional

COPY . .

ENV R_LIBS_USER=/usr/local/lib/R/site-library
WORKDIR /work
CMD ["R", "--no-save"]
