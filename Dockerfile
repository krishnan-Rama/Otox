# Base image with Ubuntu 20.04
FROM ubuntu:20.04

# Metadata
LABEL maintainer="krishnanR1@cardiff.ac.uk"
LABEL description="Bioinformatics pipeline for comparative ortholog motif and domain analysis"

# Avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install essential tools and HMMER
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    build-essential \
    python3 \
    python3-pip \
    default-jre \
    hmmer \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh
ENV PATH="/opt/conda/bin:$PATH"

# Install bioinformatics tools and Python packages
RUN conda install -c bioconda -c conda-forge -y \
    mafft \
    hmmer \
    weblogo \
    biopython \
    clustalw \
    muscle \
    iqtree \
    fasttree \
    ete3 \
    blast \
    seqkit \
    && conda clean -a -y

RUN pip install --no-cache-dir matplotlib pandas seaborn dendropy

# Add pre-existing Pfam database files into the container
COPY Pfam-A.hmm* /opt/conda/envs/bioinformatics/share/

# Add the pipeline script to the container
COPY bio_pipeline.py /usr/local/bin/bio_pipeline
RUN chmod +x /usr/local/bin/bio_pipeline

# Create a non-root user for reproducibility
RUN useradd -m pipeline_user
USER pipeline_user

# Define the working directory
WORKDIR /workspace

# Entry point for running the pipeline
ENTRYPOINT ["bio_pipeline"]

