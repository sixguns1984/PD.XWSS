# Set working directory and load libraries
setwd("/public/labdata/liaoyu/chrx/cox_new/MMSE")

library(data.table)
library(dplyr)
library(tidyr)
library(plyr)
library(survival)
library(stats)

# Load genotype data (PLINK recoded to raw format, 012 encoding)
gws_snp <- fread("merged_012.raw", data.table = FALSE)
gws_snp <- gws_snp[, -c(1, 3, 4, 5, 6)]  # Remove FID, PAT, MAT, SEX, PHENOTYPE

# Load and merge phenotype data
MMSE_pheno <- read.csv("no_haplogroup_result_haplogroup.csv")
MMSE_pheno <- na.omit(MMSE_pheno)  # Remove missing values

sample <- read.table("MEGACombine.Samples.txt", header = TRUE)
MMSE_pheno <- merge(MMSE_pheno, sample, by = "ID")
MMSE_pheno <- na.omit(MMSE_pheno)
MMSE_pheno$IID <- gsub("^PPMISI", "PP-", MMSE_pheno$IID)

# Merge genotype and phenotype by IID
combined <- merge(gws_snp, MMSE_pheno, by = "IID")

# Sex-stratified separation
# Total: 3,763 individuals (male = 2,370; female = 1,393)
male <- subset(combined, SEX == "M")
female <- subset(combined, SEX == "F")

# Display dimensions
dim(male)
dim(female)

# Convert male genotypes: recode 1 → 2
male_matrix <- as.matrix(male[, 2:213759])
male_matrix[male_matrix == 1] <- 2
dataframe_male <- as.data.frame(male_matrix)

# Reattach phenotype columns
male <- cbind(dataframe_male, male[, -c(2:213759)])

# Combine male and female data
combined <- rbind(male, female)
combined$SEX <- ifelse(combined$SEX == "M", 0, 1)  # Recode sex: M=0, F=1

# Initialize result container
result <- data.frame()

##########################
# Cox regression on combined data
for (i in colnames(combined[, 1:213758])) {
  cox <- coxph(
    Surv(DURVISIT, MMSE) ~ get(i) + frailty(STUDY_NAME) + SEX +
      AAO + YEARSEDUC + GBA.APOE4.PHS +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    data = combined,
    method = "breslow"
  )
  
  coxSummary <- summary(cox)
  print(i)
  
  result <- rbind(result, cbind(
    id = i,
    HR = coxSummary$conf.int[1, "exp(coef)"],
    HR.95L = coxSummary$conf.int[1, "lower .95"],
    HR.95H = coxSummary$conf.int[1, "upper .95"],
    SE = coxSummary$coefficients[1, "se(coef)"],
    pvalue = coxSummary$coefficients[1, "p"]
  ))
}

##########################
# Cox regression on males only
for (i in colnames(male[, 1:213758])) {
  cox <- coxph(
    Surv(DURVISIT, MMSE) ~ get(i) + frailty(STUDY_NAME) +
      AAO + YEARSEDUC + GBA.APOE4.PHS +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    data = male,
    method = "breslow"
  )
  
  coxSummary <- summary(cox)
  print(i)
  
  result <- rbind(result, cbind(
    id = i,
    HR = coxSummary$conf.int[1, "exp(coef)"],
    HR.95L = coxSummary$conf.int[1, "lower .95"],
    HR.95H = coxSummary$conf.int[1, "upper .95"],
    SE = coxSummary$coefficients[1, "se(coef)"],
    pvalue = coxSummary$coefficients[1, "p"]
  ))
}

##########################
# Cox regression on females only
female <- female[, -1]  # Remove IID column

for (i in colnames(female[, 1:213758])) {
  cox <- coxph(
    Surv(DURVISIT, MMSE) ~ get(i) + frailty(STUDY_NAME) +
      AAO + YEARSEDUC + GBA.APOE4.PHS +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    data = female,
    method = "breslow"
  )
  
  coxSummary <- summary(cox)
  print(i)
  
  result <- rbind(result, cbind(
    id = i,
    HR = coxSummary$conf.int[1, "exp(coef)"],
    HR.95L = coxSummary$conf.int[1, "lower .95"],
    HR.95H = coxSummary$conf.int[1, "upper .95"],
    SE = coxSummary$coefficients[1, "se(coef)"],
    pvalue = coxSummary$coefficients[1, "p"]
  ))
}
