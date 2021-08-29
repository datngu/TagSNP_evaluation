#!/bin/bash

# Author: Dat T Nguyen <n.dat@outlook.com>
# Date: 15 March 2021

# This script is to use to run leave one out imputation with pre-built m3vcf reference files, input is path to output directory that created by [create_imputation_ref.sh].
# We assume that bcftools, vcftools, Minimac3, and Minimac4 are properly installed and can be called directly from bash shell
#
#input parameters are:
    # -t|--tag -- input is a tab separated file, with two columns contain chr number and position of SNP and the file should have no header.  
    # -r|--ref -- path to output directory that is output of [create_imputation_ref.sh].
    # -o|--out -- output directory that contains imputation results, pre-created dir is not needed.
    # -p|--thread -- number of thread.

# /media/datn/data/DatProjects/vn_array/wraped_codes/imputation_with_prebuilt_ref.sh -t /media/datn/data/DatProjects/vn_array/build_training_chr18/train_7.txt -r /media/datn/data/DatProjects/vn_array/reference_imputation/chr_18_ref_Biallelic_SNP_only_maf_0.01 -o testing_imputation

###################################
syntax=" ./imputation_with_prebuilt_ref.sh -t [path to your tagSNP file] -r [path to your reference directory (output of create_imputation_ref.sh) ] -o [path to your output dir] -p [number of thread]"

while [[ $# -gt 1 ]]
do
key="$1"

case $key in

     -t|--tag) 
     in_tag=$(readlink -f $2)
     shift
     ;;

     -r|--ref) 
     ref_dir=$(readlink -f $2)
     shift
     ;;  

     -p|--thread) 
     CPUNUM=$2
     shift
     ;; 

     -o|--out)
     res_dir=$2
     shift
     ;;
     *)

esac
shift
done

if [[ -z "$in_tag" ]]; then
   echo ""
   echo "Usage:"
   echo $syntax
   echo ""
   echo "ERROR: no tagSNP file !"
   echo ""
   exit
fi

if [[ -z "$ref_dir" ]]; then
   echo ""
   echo "Usage:"
   echo $syntax
   echo ""
   echo "ERROR: no reference directory !"
   echo ""
   exit
fi

#####################

if [[ -z "$res_dir" ]]; then
  echo "No specific setting for output, so we create [mkdir res_dir] in your current directory.."
  res_dir="res_dir"
fi


if [[ -z "$CPUNUM" ]]; then
   CPUNUM=8
   echo "No specific setting for thread number, so use the default setting (-p 8)..."
fi


########################

mkdir $res_dir
vcftools --gzvcf $ref_dir/ref.vcf.gz --positions $in_tag --recode --stdout | bgzip > $res_dir/array.vcf.gz



function leave_one_out_imputation {
	in_SNP_array=$1
	sample=$2
  ref_dir=$3
  res_dir=$4
  mkdir $res_dir/$sample
	bcftools view $in_SNP_array -s $sample -a -o $res_dir/$sample/${sample}_SNP_array.vcf.gz -Oz
	minimac4 --refHaps $ref_dir/${sample}/${sample}_ref_panel.m3vcf.gz \
         --haps $res_dir/$sample/${sample}_SNP_array.vcf.gz \
         --prefix $res_dir/imputed_${sample} \
         --ignoreDuplicates \
         --cpus 1
}

# parallel through all samples

N=$CPUNUM
(
for i in $(cat $ref_dir/all_sample.txt); 
do 
   ((j=j%N)); ((j++==0)) && wait
   leave_one_out_imputation $res_dir/array.vcf.gz $i $ref_dir $res_dir & 
done
)

sleep 40s
# process ref to read in R
cd $res_dir
zcat $ref_dir/ref.vcf.gz | grep -v "##" | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' > ref.encoded.txt
rm *.info
rm -r `ls -d */`


for file in *dose.vcf.gz
do
	out=`echo $file | cut -d "." -f1`
	zcat $file | grep -v "^##" | awk '//{printf "%s\t%s\n", $3, $10}' | sed 's/0|0://g' | sed 's/0|1://g' | sed 's/1|0://g' | sed 's/1|1://g' > ${out}.txt
done

rm *dose.vcf.gz

echo "Done processing!"



# #######################################################################

# ref = read.delim("ref.encoded.txt")

# #x = head(ref)

# ref$INFO = sapply(ref$INFO, function(x) unlist(strsplit(x,";"))[2])
# ref$INFO = sapply(ref$INFO, function(x) unlist(strsplit(x,"="))[2])
# ref$ID = paste(ref$X.CHROM, ref$POS, ref$REF, ref$ALT, sep = ":")

# sample_names = colnames(ref)              
# sample_names = sample_names[c(10: length(sample_names))]

# ## filter non-imputed SNPs
# file = paste0("imputed_", sample_names[1], ".txt")
# tem = read.delim(file)
# pick = ref$ID %in% tem$ID
# ref = ref[pick,]
# ## load imputed files
# imputed = data.frame(ID = ref$ID)
# for(sample in sample_names){
#   file = paste0("imputed_", sample, ".txt")
#   tem = read.delim(file)
#   imputed[,sample] = tem[match(imputed$ID, tem$ID) ,sample]
# }

# data = ref[, c(1:8)]
# x = as.matrix(ref[, sample_names])
# y =  as.matrix(imputed[, sample_names])
# data$r_2 = 0
# for( i in 1 : nrow(data)){
#   data$r_2[i] = cor(x[i,], y[i,], method = "pearson")
# }
# data$r_2 = data$r_2^2

# data$INFO = as.numeric(data$INFO)
# data$MAF = data$INFO
# data$MAF[data$MAF > 0.5] = 1 - data$INFO[data$MAF > 0.5]


# cutt_off = list(
#   c(0.00, 0.01),
#   c(0.01, 0.02),
#   c(0.02, 0.03),
#   c(0.03 ,0.04),
#   c(0.04, 0.05),
#   c(0.05, 0.075),
#   c(0.075, 0.1),
#   c(0.1, 0.125),
#   c(0.125, 0.15),
#   c(0.15, 0.175),
#   c(0.175, 0.2),
#   c(0.2, 0.225),
#   c(0.225, 0.25),
#   c(0.25, 0.275),
#   c(0.275, 0.3),
#   c(0.3, 0.325),
#   c(0.325, 0.35),
#   c(0.35, 0.375),
#   c(0.375, 0.4),
#   c(0.4, 0.425),
#   c(0.425, 0.45),
#   c(0.45, 0.475),
#   c(0.475, 0.5)
# )

# res = c()
# for (i in 1:length(cutt_off) ){
#   pick = data$MAF > cutt_off[[i]][1] &  data$MAF <= cutt_off[[i]][2]
#   tem = data[pick,]
#   tem_res = mean(tem$r_2, na.rm = T)
#   names(tem_res) = paste( cutt_off[[i]][1],  cutt_off[[i]][2], sep = ":")
#   res = c(res, tem_res)
# }



# res

# save(data, res, file = "impute_res_data_chr1_train0.Rdata")


