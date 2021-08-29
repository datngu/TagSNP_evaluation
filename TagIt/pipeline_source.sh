#!/bin/bash

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 15 Aug 2021
# Dependency: vcftools, Plink 1.9
# this script is used to select tagSNPs - used for TagIt algorithm
# syntax bash this_script.sh your_file.vcf.gz tag_name
# tag_name is typically is name of chromosome
# example: bash this_script.sh EUR_chr10.vcf.gz EUR_chr10

file=`realpath $1`
base_name=$2


out_dir=${base_name}


mkdir $out_dir
cd $out_dir

plink --vcf $file \
    --vcf-half-call 'haploid' \
    --make-bed  --const-fid --out ${base_name}_tem_file \
    --threads 1 \
    --memory 2000


plink --bfile ${base_name}_tem_file \
  --r --ld-window-r2 0.2 \
  --ld-window 10000 \
  --ld-window-kb 1000 \
  --out ${base_name}_ld_normal \
  --threads 8 \
  --memory 2000

# compute fre file
vcftools --gzvcf $file --freq --out ${base_name}_frequency
# remove "chr"
sed -i 's/chr//g' ${base_name}_frequency.frq

# run a custom Rsript to process ld file to tagit input
Rscript /path/to/process_ld.R ${base_name}_frequency.frq ${base_name}_ld_normal.ld $base_name

/path/to/TagIt/src/tagit --af ${base_name}_frequency.frq --ld ${base_name}_processed_LD.txt --r2 0.8 --out-summary ${base_name}_summary.txt.gz --out-tagged ${base_name}_tagged.txt.gz --out-tags ${base_name}_tags.txt.gz

zcat ${base_name}_tags.txt.gz | sort -k 2nr | grep -v ^M | sed 's/:/\t/' | awk '{print $1,$2}' > ${base_name}_tags_cleaned.txt

