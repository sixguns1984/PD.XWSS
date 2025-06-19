#!/bin/bash
# Re-run QC on merged data, apply sex-stratified filters,
# exclude outliers, and prepare for downstream SNP comparisons and encoding.

# Step 1: Exclude SNPs in PARs or complex X regions
plink --bfile merged_final \
      --exclude range hg19_X_complicated_region \
      --make-bed \
      --out step1

# Step 2: R - Match sample information and assign sex
step1 <- read.table("D:\\lyyy\\chr_X\\cox\\2025_3_2\\step1.fam")
PC <- read.table("MEGACombine-SUBJECT-PC.txt", sep = "\t", header = TRUE)
sample <- read.table("D:\\lyyy\\mega_phenotype_data\\MEGACombine.Samples.txt", header = TRUE)

PC <- merge(PC, sample, by = "ID")
PC$IID <- gsub("PPMISI", "PP-", PC$IID)

# Add 24 additional samples (16 males, 8 females)
male <- subset(PC, SEX == "M")
female <- subset(PC, SEX == "F")
step1_male <- subset(step1, V2 %in% male$IID)
step1_female <- subset(step1, V2 %in% female$IID)

write.table(step1_male[, 1:2], "D:\\lyyy\\chr_X\\cox\\2025_3_2\\male.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(step1_female[, 1:2], "D:\\lyyy\\chr_X\\cox\\2025_3_2\\female.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Step 3: Split by sex and apply standard QC
plink --bfile step1 --keep male.txt --make-bed --out male        # 2,838 males
plink --bfile step1 --keep female.txt --make-bed --out female    # 1,653 females

# Step 4: HWE filter in females
plink --bfile female --hardy --out female_hardy
awk '{if($9 < 0.000001) print $0}' female_hardy.hwe > plinkzoomhwe.hwe
head -50 female_hardy.hwe
awk '{if($9 < 0.000001) print $2}' female_hardy.hwe > hardy_snp_exclude   # 270 SNPs to exclude

# Step 5: MAF, geno, mind filters (same exclusion list for both)
plink --bfile female \
      --maf 0.01 \
      --geno 0.05 \
      --mind 0.1 \
      --exclude hardy_snp_exclude \
      --make-bed \
      --out females_filtered   # 213,946 variants, 1,653 females retained

plink --bfile male \
      --maf 0.01 \
      --geno 0.05 \
      --mind 0.1 \
      --exclude hardy_snp_exclude \
      --make-bed \
      --out males_filtered     # 214,377 variants, 2,838 males retained

# Step 6: Identify common SNPs between sexes
awk '{print $2}' males_filtered.bim | sort | uniq >> tmp1
awk '{print $2}' females_filtered.bim | sort | uniq >> tmp2
cat tmp1 tmp2 > tmp
sort tmp | uniq -c | awk '$1 == 2 {print $2}' > common_snps

plink --bfile females_filtered --extract common_snps --make-bed --out females_common
plink --bfile males_filtered --extract common_snps --make-bed --out males_common

# Step 7: Merge sex-stratified common datasets
plink --bfile females_common --bmerge males_common --make-bed --out merged   # 213,758 SNPs, 4,491 samples

# Step 8: Remove 24 outlier samples
plink --bfile merged --remove sample_exclude.txt --make-bed --out sample_4467

# Step 9: Get SNP lists before/after exclusion
plink --bfile sample_4467 --write-snplist --out sample_snps
plink --bfile /public/labdata/liaoyu/chrx/merge_ppmi_QC_new/merged --write-snplist --out merged_snps

# Step 10: Identify shared and unique SNPs
grep -Fxf sample_snps.snplist merged_snps.snplist > before_after_exclude_common_snps.txt
grep -Fvxf merged_snps.snplist sample_snps.snplist > unique_in_last_exclude.txt

# Step 11: Extract SNPs for downstream encoding
plink --bfile sample_4467 \
      --extract before_after_exclude_common_snps.txt \
      --make-bed \
      --out before_after_exclude_common_snps

plink --bfile sample_4467 \
      --extract unique_in_last_exclude.txt \
      --make-bed \
      --out unique_in_last_exclude

# Step 12: Convert to allele dosage format (012 coding)
plink --bfile before_after_exclude_common_snps \
      --recode A \
      --allow-no-sex \
      --out before_after_exclude_common_snps_012

plink --bfile unique_in_last_exclude \
      --recode A \
      --allow-no-sex \
      --out unique_in_last_exclude_012

plink --bfile sample_4467 \
      --recode A \
      --allow-no-sex \
      --out sample_4467
