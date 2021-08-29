# TagSNP_evaluation
Evaluation tag SNP selection strategies

These scipts tested in Ubuntu v18 environment.


# Clone this project:

```sh

git clone https://github.com/datngu/TagSNP_evaluation

```
# I. RUNNING PIPELINES

Assumming that you have a processed (Bialellic SNPs with MAF >= 1%) vcf.gz file named: chr10_EUR.vcf.gz and want to compare impuatation accuracy with 30000 tag SNP selected.

## 1. Running TagIt

```sh

cd TagSNP_evaluation/TagIt 

bash bash configure.sh

TagSNP_evaluation/TagIt/TagIt_pipeline.sh chr10_EUR.vcf.gz chr10_EUR

# picking top 30000 tag SNP selected from 
head -n 30000 chr10_EUR/chr10_EUR_tags_cleaned.txt > chr10_EUR_tags_30000_cleaned.txt

```
You output should be: `TagSNP_evaluation/TagIt/chr10_EUR_tags_30000_cleaned.txt`

## 2. Running FastTagger

```sh

cd TagSNP_evaluation/FastTagger 

bash bash configure.sh

TagSNP_evaluation/TagIt/FastTagger.sh chr10_EUR.vcf.gz 30000 chr10_EUR

```

You output should be: `TagSNP_evaluation/FastTagger/cleaned_EAS_chr10_EUR_30000.tagSNP.txt`

## 3. Running EQ_Naive

```sh

cd TagSNP_evaluation/EQ_Naive

EQ_Naive.R vcf=chr10_EUR.vcf.gz size=30000 out=EQ_Naive_array_30000.txt

```

You output should be: `TagSNP_evaluation/EQ_Naive/EQ_Naive_array_30000.txt`

## 4. Running EQ_MAF

```sh

cd TagSNP_evaluation/EQ_MAF

EQ_MAF.R vcf=chr10_EUR.vcf.gz size=30000 out=EQ_MAF_array_30000.txt

```

You output should be: `TagSNP_evaluation/EQ_MAF/EQ_MAF_array_30000.txt`

# II. EVALUATION WITH LEAVE ONE OUT CROSS VALIDATION

## 1. Create imputation reference pannel
This step take very long time, depends on your size of chromosome and population.
```sh
cd TagSNP_evaluation
create_imputation_ref.sh -v chr10_EUR.vcf.gz -o chr10_EUR_imputation_ref -p 16

```


## 2. Leave one out imputation with your pre-buit reference panel and computing Pearson's correlation

```sh
cd TagSNP_evaluation
# 1. TagIt
## imputation
imputation_with_prebuilt_ref.sh -t TagIt/chr10_EUR_tags_30000_cleaned.txt -r chr10_EUR_imputation_ref -o TagIt_EUR -p 16
## computing Pearson's correlation
compute_imputation_accuracy.R imputation=TagIt_EUR out=TagIt_EUR.Rdata

# 2. FastTagger
## imputation
imputation_with_prebuilt_ref.sh -t FastTagger/cleaned_EAS_chr10_EUR_30000.tagSNP.txt -r chr10_EUR_imputation_ref -o FastTagger_EUR -p 16
## computing Pearson's correlation
compute_imputation_accuracy.R imputation=FastTagger_EUR out=FastTagger_EUR.Rdata


# 3. EQ_Naive
## imputation
imputation_with_prebuilt_ref.sh -t EQ_Naive/EQ_Naive_array_30000.txt -r chr10_EUR_imputation_ref -o EQ_Naive_EUR -p 16
## computing Pearson's correlation
compute_imputation_accuracy.R imputation=EQ_Naive_EUR out=EQ_Naive_EUR.Rdata

# 4. EQ_MAF
## imputation
imputation_with_prebuilt_ref.sh -t EQ_MAF/EQ_MAF_array_30000.txt -r chr10_EUR_imputation_ref -o EQ_MAF_EUR -p 16
## computing Pearson's correlation
compute_imputation_accuracy.R imputation=EQ_MAF_EUR out=EQ_MAF_EUR.Rdata

```

You output should be: `TagIt_EUR.Rdata, FastTagger_EUR.Rdata, EQ_Naive_EUR.Rdata, EQ_MAF_EUR.Rdata`






















