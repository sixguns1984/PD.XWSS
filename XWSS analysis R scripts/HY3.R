# Set working directory and load required libraries
setwd("/public/labdata/liaoyu/chrx/cox_new/HY3")

library(data.table)
library(dplyr)
library(tidyr)
library(plyr)
library(survival)
library(stats)

# Load genotype data (012 encoded PLINK format)
gws_snp <- fread("merged_012.raw", data.table = FALSE)
gws_snp <- gws_snp[, -c(1, 3, 4, 5, 6)]  # Remove FID, PAT, MAT, SEX, PHENOTYPE

# Load and merge phenotype data
hyresult_pheno <- read.csv("HY3_phenotype.csv")
hyresult_pheno <- na.omit(hyresult_pheno)  # Remove missing values

sample <- read.table("MEGACombine.Samples.txt", header = TRUE)
hyresult_pheno <- merge(hyresult_pheno, sample, by = "ID")
hyresult_pheno <- na.omit(hyresult_pheno)
hyresult_pheno$IID <- gsub("^PPMISI", "PP-", hyresult_pheno$IID)

# Merge genotype and phenotype data
combined <- merge(gws_snp, hyresult_pheno, by = "IID")

# Total: 3,735 individuals (male = 2,404; female = 1,331)
male <- subset(combined, SEX == "M")
female <- subset(combined, SEX == "F")

# Convert male genotypes: recode 1 → 2
male_matrix <- as.matrix(male[, 2:213759])
male_matrix[male_matrix == 1] <- 2
dataframe_male <- as.data.frame(male_matrix)

# Combine converted male genotypes with phenotypes
male <- cbind(dataframe_male, male[, -c(2:3759)])  # Ensure the rest columns match
combined <- rbind(male, female)

# Recode sex: M = 0, F = 1
combined$SEX <- ifelse(combined$SEX == "M", 0, 1)

# Prepare result storage
result <- data.frame()

#############################################
# Cox regression: Combined (male + female)
for (i in colnames(combined[, 1:213758])) {
  cox <- coxph(
    Surv(DURVISIT, HY) ~ get(i) + frailty(STUDY_NAME) + SEX +
      AAO + PC1 + PC2 + PC3 + PC4 + PC5 +
      PC6 + PC7 + PC8 + PC9 + PC10,
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

#############################################
# Cox regression: Males only
# (Assumes `male` is correctly preprocessed above)

for (i in colnames(male[, 1:213758])) {
  cox <- coxph(
    Surv(DURVISIT, HY) ~ get(i) + frailty(STUDY_NAME) +
      AAO + PC1 + PC2 + PC3 + PC4 + PC5 +
      PC6 + PC7 + PC8 + PC9 + PC10,
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

#############################################
# Cox regression: Females only
female <- female[, -1]  # Remove IID

for (i in colnames(female[, 1:213758])) {
  cox <- coxph(
    Surv(DURVISIT, HY) ~ get(i) + frailty(STUDY_NAME) +
      AAO + PC1 + PC2 + PC3 + PC4 + PC5 +
      PC6 + PC7 + PC8 + PC9 + PC10,
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
