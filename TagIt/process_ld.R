#!/usr/bin/env R

## functions for parsing frq file from vcftools

#my_agrs = c("chr22.frq","plink.ld" , "chr22")

my_agrs = commandArgs(trailingOnly=TRUE)

for( i in 1: length(my_agrs)){
	if(grepl(".frq", my_agrs[i], fix = TRUE)) frq_file= my_agrs[i]

	if(grepl(".ld", my_agrs[i], fix = TRUE)) ld_file= my_agrs[i]
}

tem = which(my_agrs %in% c(frq_file, ld_file))
out_prefix = my_agrs[-tem]
cat("\n")
cat(ld_file)
cat("\n")
cat(frq_file)


options(stringsAsFactors=FALSE)

read_frq <- function(path){
	all_content = readLines(path)
	all_content = all_content[-1]
	list_content = strsplit(all_content,"\t")
	pick = which(lengths(list_content) == 6)
	list_content=list_content[pick]
	mx = do.call("rbind", list_content)
	df = as.data.frame(mx, stringsAsFactors= F)
	V5 = strsplit(df$V5,":")
	ALLELE1 = unlist(lapply(V5,FUN=function(x) unlist(x)[1]))
	ALLELE1_frq = unlist(lapply(V5,FUN=function(x) unlist(x)[2]))
	V6 = strsplit(df$V6,":")
	ALLELE2 = unlist(lapply(V6,FUN=function(x) unlist(x)[1]))
	ALLELE2_frq = unlist(lapply(V6,FUN=function(x) unlist(x)[2]))

	df = df[,c(1:4)]
	colnames(df) = c("CHROM", "POS", "N_ALLELES", "N_CHR")
	df$ALLELE1 = ALLELE1
	df$ALLELE1_frq = ALLELE1_frq
	df$ALLELE2 = ALLELE2
	df$ALLELE2_frq = ALLELE2_frq

	return(df)
}





process_ld <- function(fq_path, ld_path, out_prefix = "outfile"){

	# fq_path= "/home/datn/Dat/chr22_test/chr22.frq"
	# ld_path= "/home/datn/Dat/chr22_test/plink.ld"
	# out_prefix = "outfile"

	out_name = paste0(out_prefix, "_processed_LD.txt")


	fq = read_frq(fq_path)
	fq$CHROM = gsub("chr", "", fq$CHROM)
	fq$al = paste(fq$ALLELE1, fq$ALLELE2, sep = "/")
	fq$id = paste(fq$CHROM, fq$POS, sep = ":")


	ld = read.table(ld_path, header = TRUE)
	ld$marker1 = paste(ld$CHR_A, ld$BP_A, sep = ":")
	ld$al1 = fq$al[match(ld$marker1, fq$id)]
	ld$AF1 = fq$ALLELE2_frq[match(ld$marker1, fq$id)]



	ld$marker2 = paste(ld$CHR_B, ld$BP_B, sep = ":")
	ld$al2 = fq$al[match(ld$marker2, fq$id)]
	ld$AF2 = fq$ALLELE2_frq[match(ld$marker2, fq$id)]
	ld$rs2 = fq$ALLELE2_frq[match(ld$marker2, fq$id)]

	res = data.frame(
		MARKER1 = paste(ld$marker1, ld$al1, sep = "_"),
		MARKER2 = paste(ld$marker2, ld$al2, sep = "_"),
		AF1 = ld$AF1,
		AF2 = ld$AF2,
		R2 = ld$R^2,
		R = ld$R)

	res = res[!is.na(res$AF2),]

	write.table(res , sep = "\t", row.names = F, col.names = T, quote = F, file = "tem.txt")
	command = paste0("echo -n '#' >> ", out_name)
	system(command)
	command = paste0("cat tem.txt >> ", out_name)
	system(command)
	return(0)
	#system("rm tem.txt")
}

#echo -n '#' >> $out_name
#cat tem.txt >> $out_name

if(length(out_prefix) > 0){
	tes = process_ld(frq_file, ld_file, out_prefix)
}else{
	tes = process_ld(frq_file, ld_file)
}	