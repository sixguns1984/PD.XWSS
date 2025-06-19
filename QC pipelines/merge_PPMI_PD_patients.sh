#!/bin/bash

set -e

# Extract SNPs on chromosome 23, update sample IDs for PPMI,
# match with MEGA samples, filter for shared SNPs, and merge datasets

# Step 1: Extract SNPs on chromosome 23
/public2/data/liaoyu/plink/plink \
  --bfile PPMI.All.SNP \
  --chr 23 \
  --make-bed \
  --out PPMI_ALL_SNP_23chr

# Step 2: Update FAM file sample IDs (keep BED/BIM unchanged)
PPMI <- read.table("D:\\lyyy\\chr_X\\cox\\2025_3_2\\PPMI_ALL_SNP_23chr.fam")
PPMI$IID <- paste("PP-", PPMI$V2, sep = "")
mm <- MEGA_SAMPLE[grepl("^PP", MEGA_SAMPLE$IID), ]
write.table(PPMI[, c(7, 7, 3:6)],
            "D:\\lyyy\\chr_X\\cox\\2025_3_2\\PPMI_ALL_SNP_23chr_newID.fam",
            sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Step 3: Generate sample filter for 471 individuals
mm$IID <- gsub("^PPMISI", "PP-", mm$IID)
TT <- merge(PPMI, mm, by = "IID")

# Step 4: Extract 471 individuals from PPMI
/public2/data/liaoyu/plink/plink \
  --bfile PPMI_ALL_SNP_23chr_newID \
  --keep sample_filter.txt \
  --make-bed \
  --out PPMI_ALL_SNP_23chr_417sample

# Step 5: Generate SNP keys from merged.bim
awk '{
  chr = $1; pos = $4; a1 = $5; a2 = $6;
  if (a1 > a2) { alleles = a1 ":" a2 } else { alleles = a2 ":" a1 };
  key = chr ":" pos ":" alleles;
  print $2, key;
}' merged.bim > merge_snp_keys.txt

# Step 6: Generate SNP keys from PPMI bim file
awk '{
  chr = $1; pos = $4; a1 = $5; a2 = $6;
  if (a1 > a2) { alleles = a1 ":" a2 } else { alleles = a2 ":" a1 };
  key = chr ":" pos ":" alleles;
  print $2, key;
}' PPMI_ALL_SNP_23chr_417sample.bim > ppmi_snp_keys.txt

# Step 7: Identify shared SNP keys
sort -k2 merge_snp_keys.txt > merge_snp_keys_sorted.txt
sort -k2 ppmi_snp_keys.txt > ppmi_snp_keys_sorted.txt
join -1 2 -2 2 merge_snp_keys_sorted.txt ppmi_snp_keys_sorted.txt > common_snp_keys.txt

# Extract matching SNPs from merged.bim
awk 'NR==FNR {keys[$1]; next} $2 in keys {print $1}' common_snp_keys.txt merge_snp_keys.txt > merge_snp_filter.txt

# Extract matching SNPs from ppmi.bim
awk 'NR==FNR {keys[$1]; next} $2 in keys {print $1}' common_snp_keys.txt ppmi_snp_keys.txt > ppmi_snp_filter.txt

# Filter merged and PPMI datasets by shared SNPs
plink --bfile merged \
      --extract merge_snp_filter.txt \
      --make-bed \
      --out merge_common

plink --bfile PPMI_ALL_SNP_23chr_417sample \
      --extract ppmi_snp_filter.txt \
      --make-bed \
      --out ppmi_common

# Step 8: Create SNP name mapping file (from merge to PPMI)
awk '{print $2, $3}' common_snp_keys.txt > snp_name_map.txt

# Step 9: Rename SNPs in merged dataset
plink --bfile merge_common \
      --update-name snp_name_map.txt \
      --make-bed \
      --out merge_common_renamed

# Step 10: Merge renamed merged dataset with PPMI dataset
plink --bfile merge_common_renamed \
      --bmerge ppmi_common \
      --make-bed \
      --out merged_final   # 214,818 SNPs, 4,491 samples retained
