#!/usr/bin/env Rscript

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 15 March 2021

# This script is create a naive SNP slection set with pre-set number of target SNP. Procedure: 1. Ranking all SNP positions 2. Pick SNPs by odering with pre-set total number. Reading the implemented codes for more details.
#input parameters are:
    # vcf=path/to/vcf.gz
    # size=number_of_tagSNP : for example: size=20000
    # out=out_file_name.txt : output file name, defaut value is: "naive_array_size.txt"

# ./EQ_MAF.R vcf=/media/datn/data/DatProjects/vn_array/data/vn_504_maf_0.01/MAF_0.01.chr22_vn504_Biallelic.SNP_PHASED_SHAPEIT4_no_chr_version_dbSNP_151.vcf.gz size=20000 out=chr18_EQ_MAF_20000.txt

args = commandArgs(trailingOnly=TRUE)
syntax='\nUsage:\t./EQ_MAF.R vcf=path/to/vcf.gz size=size out=output_name.txt\n\nInput parameters are:\n\tvcf=path/to/vcf.gz\n\tsize=number_of_tagSNP : for example: size=20000\n\tout=out_file_name.txt : output file name, defaut value is: "EQ_MAF_array_size.txt"\n\n'

vcf_path = array_size = out_name = NA

if(length(args) == 0 ){
  cat("\nNo argument, Program stop! \n")
  cat(syntax)
  quit()
}

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1] == "vcf") vcf_path = res[2]
  if (res[1] == "size") array_size = as.numeric(res[2])
  if (res[1] == "out") out_name = res[2]   
}

if(is.na(vcf_path)){
  cat(syntax)
  cat("vcf is not found! Program stop!\n\n")
  quit()
}

if(is.na(array_size)){
  cat(syntax)
  cat("size is not found! Program stop!\n\n")
  quit()
}

# vcf_path = "/media/datn/data/DatProjects/vn_array/reference_imputation/chr_10_ref_Biallelic_SNP_only_maf_0.01/ref.vcf.gz"
# array_size = 30000

# extract snp positions
all_site = paste0(vcf_path, "_all_site.txt")
command = paste0("zcat ", vcf_path, " | grep -v ^# | awk '{print $1, $2, $8}' > ", all_site)
system(command)

#all_site = "/media/datn/data/DatProjects/vn_array/reference_imputation/chr_10_ref_Biallelic_SNP_only_maf_0.01/ref.vcf.gz_all_site.txt"

site = read.table(all_site, header = F)

info_slipt = unlist(strsplit(site$V3[1],";", fixed = T))
pick_pos = grep("AF=", info_slipt, fixed = T)[1]


site$V3 = sapply(site$V3, function(x) unlist(strsplit(x,";"))[pick_pos])
site$V3 = sapply(site$V3, function(x) unlist(strsplit(x,"="))[2])
#site$ID = paste(site$X.CHROM, site$POS, site$REF, site$ALT, sep = ":")
site$V3  = as.numeric(site$V3 )

site$MAF = site$V3
site$MAF[site$V3 > 0.5] = 1 - site$V3[site$V3 > 0.5]

min_pos = min(site$V2)
max_pos = max(site$V2)
window_size = (max_pos - min_pos)/ array_size
end_points = rep(window_size, array_size)
end_points = cumsum(end_points) + min_pos
start_points = end_points - window_size
start_points = round(start_points)
end_points = round(end_points)

array = site[1,]
n_pick = 1
for(i in 1:array_size){
  pick = site$V2 >= start_points[i] & site$V2 < end_points[i]
  if(length(pick) == 0){
    n_pick = n_pick + 1
    next
  }
  tem = site[pick,]
  od = order(tem$MAF, decreasing = T)
  tem = tem[od,]

  if(n_pick <= nrow(tem)){
    res = tem[c(1:n_pick),]
    n_pick = 1
  }else{
    res = tem 
    n_pick = n_pick - nrow(tem) + 1
  }
  array = rbind(array, res)
}

array = array[-1,c(1:2)]

if ( is.na(out_name)) out_name = paste0( "EQ_MAF_array_", (array_size), ".txt")

write.table(array , file = out_name, col.names = F, row.names = F, sep ="\t", quote = F)

cat(paste0("\nOutput: ", out_name, "\n"))
cat("Done!\n")