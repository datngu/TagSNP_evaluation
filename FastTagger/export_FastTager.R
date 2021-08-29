#!/usr/bin/env Rscript

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 15 Aug 2021

# ./export_FastTager.R maf_path FastTager_tagSNP_output_path

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE)

syntax='\nUsage:\t./export_FastTager.R maf_file_path FastTager_tagSNP_file_path\n\n'

if(length(args) == 0 ){
  cat("\nNo argument, Program stop! \n")
  cat(syntax)
  quit()
}

maf_path = args[1]
tag_path = args[2]
maf = read.delim(maf_path, header = F)
tag = read.delim(tag_path, header = F)
chr = 10
pick = tag$V1 +1
final_tag = maf[pick,c(1,2)]
final_tag$V1 = chr
#write.table(final_tag, file = paste0("chr_", chr, "_FastTagger_cleaned.txt"), quote = F, col.names = F, row.names = F, sep = "\t")
#x = dirname(tag_path)
y = basename(tag_path)
out_name = paste0( "cleaned_", y)




write.table(final_tag, file = out_name, quote = F, col.names = F, row.names = F, sep = "\t")

cat(paste0("\nOutput: ", out_name, "\n"))
cat("\nDone!\n")