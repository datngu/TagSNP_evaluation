#!/usr/bin/env Rscript

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 15 Aug 2021

# ./preprocessing_FastTagger.R vcf=/path/file.vcf.gz out_prefix=out

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE)
syntax='\nUsage:\t./preprocessing_FastTagger.R vcf=/path/file.vcf.gz out_prefix=out\n\n'

vcf_path = out_name = NA

if(length(args) == 0 ){
  cat("\nNo argument, Program stop! \n")
  cat(syntax)
  quit()
}

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1] == "vcf") vcf_path = res[2]
  #if (res[1] == "size") array_size = as.numeric(res[2])
  if (res[1] == "out_prefix") out_prefix = res[2]   
}

if(is.na(vcf_path)){
  cat(syntax)
  cat("vcf is not found! Program stop!\n\n")
  quit()
}

if ( is.na(out_prefix)) out_prefix = "preprocessing_FastTagger"
out_matrix = paste0(out_prefix, ".matrix")
out_maf = paste0(out_prefix, ".maf")

# vcf_path = file.vcf.gz

cat("\nProcessing vcf.gz file! \n")

cmd = paste0("zcat ", vcf_path , " | grep -v ^# | sed 's/|/\t/g' | sed 's/chr//g' > encoded_tempo_vcf.txt")

system(cmd)

cat("\nReading tempo file to R! \n")

ref = read.table("encoded_tempo_vcf.txt")
x = head(ref$V8)
AF = grep( "AF=", unlist(strsplit(x[1],";")), fixed = T)[1]
ref$V8 = sapply( ref$V8,function(x) unlist(strsplit(x,";"))[AF])
ref$V8 = gsub("AF=", "", ref$V8)

cat("\nWriting matrix file! \n")
matrix = ref[,-c(1:9)]
write.table(matrix, quote = F, sep = " ", row.names = F, col.names = F, file = out_matrix)

cat(paste0("\nOutput: ", out_matrix, "\n"))

maf = ref[,c(1,2,4,5,8)]  
maf[,1] = paste0("SNP", c(0 : (nrow(maf)-1) ))

cat("\nWriting maf file! \n")
write.table(maf, quote = F, sep = "\t", row.names = F, col.names = F, file = out_maf)


cat(paste0("\nOutput: ", out_maf, "\n"))

system("rm encoded_tempo_vcf.txt")

cat("\nDone!\n")