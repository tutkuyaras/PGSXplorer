# Start with a base image that includes R
FROM rocker/r-ver:4.1.0

# Install dependencies and system utilities
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    unzip \
    build-essential \
    python3 \
    python3-pip \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    zlib1g-dev

# Install PLINK v1.9
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip \
    && unzip plink_linux_x86_64_20210606.zip -d /usr/local/bin/ \
    && rm plink_linux_x86_64_20210606.zip

# Install PLINK v2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20241009.zip \
    && unzip plink2_linux_x86_64_20241009.zip -d /usr/local/bin/ \
    && rm plink2_linux_x86_64_20241009.zip

# Install Python packages
RUN pip3 install --upgrade pip \
    && pip3 install numpy pandas scipy matplotlib

# Install R packages for analysis
RUN Rscript -e 'install.packages(c("data.table", "ggplot2", "dplyr", "optparse"), repos="http://cran.rstudio.com/")'

# Create a directory for your pipeline
WORKDIR /PGSXplorer

# Ana Nextflow dosyasını ve config dosyasını kopyala
COPY main.nf /PGSXplorer/
COPY nextflow.config /PGSXplorer/

# R scriptlerini içeren bin klasörünü kopyala
COPY bin/ /PGSXplorer/bin/

# Nextflow modüllerini kopyala
COPY modules/ /PGSXplorer/modules/

# Set the default command for running the pipeline
CMD ["nextflow", "run", "main.nf"]
