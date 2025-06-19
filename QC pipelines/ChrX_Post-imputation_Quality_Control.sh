#!/bin/bash

set -e

# Merge imputed genotype data and perform quality control on ChrX,
# including HWE filtering, MAF/missing rate filtering by sex,
# SNP intersection, and merging with external cohort (PPMI)

plink --merge-list merge_genotype_imputation.txt --recode --out merge_genotyped_imputation

cd /public/labdata/liaoyu/after_imputation_QC

plink --bfile /public/labdata/liaoyu/imputation/chr_X/merge_genotyped_imputation \
      --exclude range hg19_X_complicated_region \
      --make-bed \
      --out step1

#######

cd /public/labdata/liaoyu/chrx/last_delete_24/ChrX_Post-imputation_Quality_Control

# Split by sex (3 study names standardized to "MEGA")
plink --bfile step1 --keep male_4020.txt --make-bed --out male
plink --bfile step1 --keep female_4020.txt --make-bed --out female

# HWE test in females
plink --bfile female --hardy --out female_hardy
awk '{if($9 < 0.000001) print $0}' female_hardy.hwe > plinkzoomhwe.hwe     # Exclude 6 SNPs
head -50 female_hardy.hwe
awk '{if($9 < 0.000001) print $2}' female_hardy.hwe > hardy_snp_exclude   # Extract SNP IDs

# Filter female data
plink --bfile female \
      --maf 0.01 \
      --geno 0.05 \
      --mind 0.1 \
      --exclude hardy_snp_exclude \
      --make-bed \
      --out females_filtered   # 224,068 variants, 1,472 samples retained

# Filter male data
plink --bfile male \
      --maf 0.01 \
      --geno 0.05 \
      --mind 0.1 \
      --exclude hardy_snp_exclude \
      --make-bed \
      --out males_filtered     # 224,300 variants, 2,548 samples retained

# Get common SNPs between males and females
awk '{print $2}' males_filtered.bim | sort | uniq >> tmp1
awk '{print $2}' females_filtered.bim | sort | uniq >> tmp2
cat tmp1 tmp2 > tmp
sort tmp | uniq -c | awk '$1 == 2 {print $2}' > common_snps

# Extract common SNPs
plink --bfile females_filtered --extract common_snps --make-bed --out females_common
plink --bfile males_filtered --extract common_snps --make-bed --out males_common

# Merge male and female datasets
plink --bfile females_common --bmerge males_common --make-bed --out merged   # 4,020 individuals, 218,593 SNPs

cd /public/labdata/liaoyu/PPMI_chrx

# Extract PD cases from PPMI cohort
plink --bfile PPMI_23 --keep PPMI_PD_ID.txt --make-bed --out PPMI_PD_23
head PPMI_PD_23.bim
wc -l PPMI_PD_23.bim

# Filter merged MEGA data by PPMI SNP list
plink --bfile /public/labdata/liaoyu/after_imputation_QC/merged \
      --extract ppmi_merged_snp111.txt \
      --make-bed \
      --out merge_filter

cd /public/labdata/liaoyu/PPMI_chrx

# Recode filtered PPMI PD data
plink --bfile ppmi_PD23_snpfilter --recode --out ppmi_PD23_snpfilter

cd /public/labdata/liaoyu/after_imputation_QC

# Recode MEGA filtered data
plink --bfile merge_filter --recode --out merge_filter

# Final merge of PPMI and MEGA samples
plink --merge-list merge_list_sample.txt --recode --out all_sample
