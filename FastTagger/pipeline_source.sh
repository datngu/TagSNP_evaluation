#!/bin/bash

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 15 Aug 2021

# this script is used to select tagSNPs - used for FastTagger algorithm
# syntax bash this_script.sh your_file.vcf.gz number_of_tagSNP base_name
# example: bash this_script.sh EUR_chr10.vcf.gz 30000 EUR_chr10


in_vcf=$1
i=$2
base_name=$3

sample_param=/path/to/param_example.txt
cat $sample_param | sed "s/max_tagSNP_num=30000/max_tagSNP_num=${i}/g" > ${i}_param.txt
/path/to/preprocessing_FastTagger.R vcf=${in_vcf} out_prefix=${base_name}
mkdir ./N_${i}
/path/to/FastTaggerV2 ${base_name}.maf ${base_name}.matrix ./N_${i}/${base_name}_${i} ${i}_param.txt
/path/to/export_FastTager.R ${base_name}.maf ./N_${i}/${base_name}_${i}.tagSNP.txt





