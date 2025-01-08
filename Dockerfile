FROM tutkuyaras/pgsxplorer_image:latest

RUN apt-get update && apt-get install -y \
    wget \
    curl \
    bash \
    bzip2 \
    build-essential \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN wget -O /usr/local/bin/beagle.jar https://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar && \
    echo '#!/bin/bash\njava -jar /usr/local/bin/beagle.jar "$@"' > /usr/local/bin/beagle && \
    chmod +x /usr/local/bin/beagle

RUN wget -O /tmp/eagle.tar.gz https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz && \
    tar -xzvf /tmp/eagle.tar.gz -C /usr/local/bin/ && \
    chmod +x /usr/local/bin/Eagle_v2.4.1/eagle && \
    rm /tmp/eagle.tar.gz

RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
    tar -xjf bcftools-1.17.tar.bz2 && \
    cd bcftools-1.17 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools-1.17 bcftools-1.17.tar.bz2

ENV PATH="/usr/local/bin/Eagle_v2.4.1:$PATH"

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/*
