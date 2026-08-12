FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    bash \
    wget \
    bzip2 \
    ca-certificates \
    make \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p /opt/conda && \
    rm miniforge.sh

ENV PATH=/opt/conda/bin:$PATH

RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10 \
    numpy \
    scipy \
    matplotlib \
    pysam \
    fastp \
    fastqc \
    bcftools \
    freebayes \
    star \
    snpeff \
    multiqc \
    stringtie \
    rmats2sashimiplot \
    samtools \
    subread \
    rmats \
    gffread \
    libsvm \
    cpc2 \
    bwa \
    ciri2 \
    transdecoder \
    r-base=4.5.2 \
    r-data.table \
    r-ggplot2 \
    r-ggrepel \
    r-scales \
    r-stringr \
    r-circlize \
    r-dplyr \
    r-tidyr \
    r-tibble \
    r-readr \
    r-patchwork \
    r-purrr \
    r-tidyverse \
    r-corrplot \
    bioconductor-complexheatmap \
    bioconductor-rtracklayer \
    r-gridbase \
    bioconductor-deseq2 \
    bioconductor-apeglm \
    && mamba clean -a -y

RUN wget https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz -O CPC2_standalone.tar.gz && \
    tar -xzf CPC2_standalone.tar.gz && \
    cd CPC2_standalone-1.0.1/libs/libsvm && \
    tar -xzf libsvm-3.18.tar.gz && \
    libsvm_dir=$(ls -d libsvm-*) && \
    cd $libsvm_dir && \
    make && \
    cp svm-predict svm-scale ../../../../CPC2_standalone-1.0.1/bin/ && \
    cd / && \
    mv /CPC2_standalone-1.0.1 /opt/cpc2 && \
    ln -s /opt/cpc2/bin/CPC2.py /usr/local/bin/CPC2.py && \
    rm /CPC2_standalone.tar.gz
