FROM continuumio/miniconda3:24.1.2-0

LABEL org.opencontainers.image.source="https://github.com/barbarahelena/osmotool"
LABEL org.opencontainers.image.description="osmotool: screen osmoadaptation genes in metagenomic datasets"
LABEL org.opencontainers.image.licenses="MIT"

# --- system deps ---
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
    && rm -rf /var/lib/apt/lists/*

# --- conda environment ---
COPY environment.yml /opt/osmotool/environment.yml
RUN conda env create -f /opt/osmotool/environment.yml \
    && conda clean -afy

# --- install osmotool package ---
# Version is supplied at build time (docker-publish.yml passes the release
# tag); local `docker build` without --build-arg VERSION falls back to this.
ARG VERSION=0.0.0+unknown
ENV SETUPTOOLS_SCM_PRETEND_VERSION_FOR_OSMOTOOL=${VERSION}
COPY . /opt/osmotool/
RUN conda run -n osmotool pip install --no-deps /opt/osmotool

# Make conda env the default
ENV PATH="/opt/conda/envs/osmotool/bin:$PATH"

ENTRYPOINT ["osmotool"]
CMD ["--help"]
