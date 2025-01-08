#!/bin/bash -ue
zcat phasedvcf_chr1.vcf.gz | grep -v '^#' | awk '{print $1, $2, $3}' > rsid_map_chr1.txt
