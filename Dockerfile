FROM r-base:4.2.1
LABEL Maintainer="Sean Jungbluth, jungbluth.sean@gmail.com" Version=1.0

# build
# docker build --platform=linux/amd64 -t rcrux:latest .

# run
# docker run --platform=linux/amd64 -it rcrux:latest /bin/bash

RUN apt-get -y -m update && apt-get install -y git libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libgit2-dev libpng-dev libtiff5-dev libjpeg-dev libmagick++-dev subversion pandoc

ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH

RUN echo "install.packages(\"devtools\", dependencies=TRUE, repos=\"https://cran.rstudio.com\")" | R --no-save

RUN echo "BiocManager::install(\"phyloseq\")" | R --no-save

RUN apt-get update && apt-get install -y libudunits2-dev libgdal-dev

RUN echo "install.packages(\"plotrix\", dependencies=TRUE, repos=\"https://cran.rstudio.com\")" | R --no-save

RUN echo "remotes::install_github(\"cpauvert/psadd\")" | R --no-save

RUN echo "devtools::install_github(\"CalCOFI/rCRUX\", build_vignettes = TRUE)" | R --no-save

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.1/ncbi-blast-2.14.1+-x64-linux.tar.gz && \
    tar -xvzf ncbi-blast-2.14.1+-x64-linux.tar.gz && \
    rm ncbi-blast-2.14.1+-x64-linux.tar.gz

ENV PATH=/ncbi-blast-2.14.1+/bin:$PATH

WORKDIR "/mnt"
