# Colocalization analysis between GWAS SNPs and GTEx eQTL signals (e.g., rs142724191)

library(tidyr)
library(coloc)
library(dplyr)
library(readxl)

# Load filtered eQTL data
combined_data <- read.csv("D:\\lyyy\\chr_X\\cox\\male02\\GTEx\\49tissue\\combined_data_filter.csv", header = TRUE)
combined_data <- separate(combined_data, col = 2, into = c("chr", "pos", "ref", "alt", "human"), sep = "_")

# Load significant SNPs from MMSE and HY3 analyses
mmse_snp <- read.table("D:\\lyyy\\chr_X\\cox\\2025_4_17\\MMSE_fuma\\significant_snp.txt", header = TRUE)
hy_snp <- read.table("D:\\lyyy\\chr_X\\cox\\2025_4_17\\HY3_fuma\\significant_snp.txt", header = TRUE)

combined_snp <- rbind(mmse_snp, hy_snp)
colnames(combined_snp)[3] <- "hg19_location"

# Join with hg38 location mapping
information <- read.table("D:\\lyyy\\chr_X\\cox\\eqtl\\location_information.txt", header = TRUE)
combined_snp <- left_join(combined_snp, information, by = "hg19_location")

# Two SNPs not matched in hg38, manually correct from NCBI
hg38_NA_snp <- combined_snp[is.na(combined_snp$hg38_location), ]
hg38_NA_snp$hg38_location <- c(14055340, 102575658)

hg38_noNA_snp <- combined_snp[!is.na(combined_snp$hg38_location), ]
hg38_snp_all <- rbind(hg38_noNA_snp, hg38_NA_snp)

# Filter GTEx SNPs that match these positions
test <- subset(combined_data, pos %in% hg38_snp_all$hg38_location)
intersect_snp <- subset(hg38_snp_all, hg38_location %in% test$pos)

# rs142724191 located in Xp11.22, associated with syndromic X-linked intellectual disability

### Define ±1MB region centered on rs142724191
pos4 <- information[information$V2 == "rs142724191", 1]
pos4_low <- pos4 - 1000000
pos4_high <- pos4 + 100000

# ========== GWAS region extraction (combined data) ==========
combined_data <- read.table("D:\\lyyy\\chr_X\\cox\\2025_3_2\\cox_summary\\MMSE_cox/merged_cox_result/combined.txt")
bim <- read.table("D:\\lyyy\\chr_X\\cox\\2025_3_2\\cox_summary\\sample_4467.bim", header = FALSE)
colnames(bim)[2] <- "id"

combined_data <- merge(bim, combined_data, by = "id")
information <- read.table("D:\\lyyy\\chr_X\\cox\\2025_3_2\\cox_summary\\significant_SNP_hg19to38\\convert_result.txt", header = TRUE)
colnames(information)[2:3] <- c("hg19_location", "hg38_location")

colnames(combined_data)[4] <- "hg19_location"
combined_data <- merge(combined_data, information, by = "hg19_location")

# Prepare GWAS input
combined_data$beta <- log(combined_data$HR)
combined_data$varbeta <- combined_data$SE^2
colnames(combined_data)[2] <- "rsid"
combined_data <- combined_data[, c("rsid", "hg38_location", "beta", "varbeta", "HR")]

# Subset region around rs142724191
combined_data3 <- subset(combined_data, hg38_location >= pos4_low & hg38_location <= pos4_high)

# ========== Reload GTEx eQTL (full) and extract Nerve_Tibial ==========
combined_data <- read.csv("D:\\lyyy\\chr_X\\cox\\male02\\GTEx\\49tissue\\combined_data_filter.csv", header = TRUE)
combined_data <- separate(combined_data, col = 2, into = c("chr", "pos", "ref", "alt", "human"), sep = "_")

Nerve <- subset(combined_data, type == "Nerve_Tibial.v8.signif_variant_gene_pairs.txt.gz")

# Prepare datasets for colocalization
gwas <- combined_data3
eqtl <- Nerve

input <- merge(
  gwas, eqtl,
  by.x = "hg38_location", by.y = "pos",
  all = FALSE,
  suffixes = c("_gwas", "_eqtl")
)

# Run coloc.abf
result <- coloc.abf(
  dataset1 = list(type = "cc", beta = input$beta, varbeta = input$varbeta),
  dataset2 = list(pvalues = input$pval_nominal, type = "quant", N = 838),
  MAF = input$maf
)

need_result <- result$results

# Extract key SNP result
subset(gwas, rsid == "rs142724191")
subset(need_result, snp == "SNP.23")
# write.xlsx(need_result, "D:\\lyyy\\chr_X\\RDD\\Colocalization.xlsx")

# ========== Repeat for male-specific data ==========

combined_data <- read.table("D:\\lyyy\\chr_X\\cox\\2025_3_2\\cox_summary\\MMSE_cox/merged_cox_result/male.txt")
bim <- read.table("D:\\lyyy\\chr_X\\cox\\2025_3_2\\cox_summary\\sample_4467.bim", header = FALSE)
colnames(bim)[2] <- "id"

combined_data <- merge(bim, combined_data, by = "id")
information <- read.table("D:\\lyyy\\chr_X\\cox\\2025_3_2\\cox_summary\\significant_SNP_hg19to38\\convert_result.txt", header = TRUE)
colnames(information)[2:3] <- c("hg19_location", "hg38_location")

colnames(combined_data)[4] <- "hg19_location"
combined_data <- merge(combined_data, information, by = "hg19_location")

# Prepare GWAS input
combined_data$beta <- log(combined_data$HR)
combined_data$varbeta <- combined_data$SE^2
colnames(combined_data)[2] <- "rsid"
combined_data <- combined_data[, c("rsid", "hg38_location", "beta", "varbeta", "HR")]

# Subset 2MB region
combined_data3 <- subset(combined_data, hg38_location >= pos4_low & hg38_location <= pos4_high)

# Reload eQTL again (GTEx nerve tissue)
combined_data <- read.csv("D:\\lyyy\\chr_X\\cox\\male02\\GTEx\\49tissue\\combined_data_filter.csv", header = TRUE)
combined_data <- separate(combined_data, col = 2, into = c("chr", "pos", "ref", "alt", "human"), sep = "_")
Nerve <- subset(combined_data, type == "Nerve_Tibial.v8.signif_variant_gene_pairs.txt.gz")

# Merge and run coloc again
gwas <- combined_data3
eqtl <- Nerve

input <- merge(
  gwas, eqtl,
  by.x = "hg38_location", by.y = "pos",
  all = FALSE,
  suffixes = c("_gwas", "_eqtl")
)

result <- coloc.abf(
  dataset1 = list(type = "cc", beta = input$beta, varbeta = input$varbeta),
  dataset2 = list(pvalues = input$pval_nominal, type = "quant", N = 838),
  MAF = input$maf
)

need_result <- result$results

# Output
subset(gwas, rsid == "rs142724191")
subset(need_result, snp == "SNP.23")
write.xlsx(need_result, "D:\\lyyy\\chr_X\\RDD\\Colocalization_male.xlsx")
