# MMSE-related survival analysis (sex-stratified and combined) using selected independent SNPs

library(lme4)
library(data.table)
library(dplyr)
library(tidyr)
library(plyr)
library(survival)
library(stats)
library(qqman)
library(survminer)
library(scales)
library(readxl)

# Set working directory
setwd("D:\\lyyy\\chr_X\\cox\\2025_3_2\\cox_summary\\cov_adjust_cox\\MMSE")

# Load independent genome-wide significant SNPs
independent_snp <- read_xlsx("D:\\lyyy\\chr_X\\cox\\2025_3_2\\cox_summary\\MMSE_cox\\fuma_output\\independent_snp.xlsx")
colnames(independent_snp)[10] <- "type"  # Rename 10th column for grouping

# Load 012-coded genotype + covariates data
gws_snp <- read.table("independent_snp_cor_cox.txt")
gws_snp[, 1:11][gws_snp[, 1:11] == 2] <- 1  # Recode 2 → 1 (simplify additive model)

####==============================
#### Combined analysis
combined <- subset(independent_snp, type == "MMSE_combined")
test1 <- gws_snp[, combined$rsID]
test1 <- cbind(test1, gws_snp[12:31])  # Add phenotype & covariates

result <- data.frame()

for (i in colnames(test1[1:nrow(combined)])) {
  cox <- coxph(
    Surv(DURVISIT, MMSE) ~ get(i) + AAO + SEX + frailty(STUDY_NAME) + 
      YEARSEDUC + GBA.APOE4.PHS + PC1 + PC2 + PC3 + PC4 + PC5 +
      PC6 + PC7 + PC8 + PC9 + PC10,
    data = test1,
    method = "breslow"
  )
  
  coxSummary <- summary(cox)
  print(i)
  
  result <- rbind(result, cbind(
    id = i,
    HR = coxSummary$conf.int[1, "exp(coef)"],
    HR.95L = coxSummary$conf.int[1, "lower .95"],
    HR.95H = coxSummary$conf.int[1, "upper .95"],
    pvalue = coxSummary$coefficients[1, "p"]
  ))
}

# Kaplan-Meier plot for one representative SNP (e.g., rs144112368)
res.cox <- coxph(
  Surv(DURVISIT, MMSE) ~ rs144112368 + AAO + SEX + YEARSEDUC + GBA.APOE4.PHS +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = test1,
  method = "breslow"
)

GG <- ggadjustedcurves(
  res.cox, variable = "rs144112368", data = test1,
  xlab = "Years since onset of PD (years)",
  ylab = "Patients surviving free of global cognitive impairment"
) +
  scale_colour_manual(values = c("#CC79A7", "#0072B2")) +
  geom_hline(aes(yintercept = 0.50), colour = "black", linetype = "dashed") +
  annotate(
    "text", x = 5.0, y = 0.30,
    label = paste(
      "Hazard ratio:", round(as.numeric(result[3, 2]), 3),
      "\n(95% CI:", round(as.numeric(result[3, 3]), 3), "-", round(as.numeric(result[3, 4]), 3), 
      ")\nP =", signif(as.numeric(result[3, 5]), 3)
    ),
    colour = "black", size = 5
  ) +
  scale_x_continuous(breaks = seq(0, 12, 2))

####==============================
#### Male-only analysis
male <- subset(independent_snp, type == "MMSE_male")
test1 <- gws_snp[, male$rsID]
test1 <- cbind(test1, gws_snp[12:31])
test1 <- subset(test1, SEX == "M")  # Filter male samples

result <- data.frame()

for (i in colnames(test1[1:nrow(male)])) {
  cox <- coxph(
    Surv(DURVISIT, MMSE) ~ get(i) + AAO + frailty(STUDY_NAME) + YEARSEDUC +
      GBA.APOE4.PHS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    data = test1,
    method = "breslow"
  )
  
  coxSummary <- summary(cox)
  print(i)
  
  result <- rbind(result, cbind(
    id = i,
    HR = coxSummary$conf.int[1, "exp(coef)"],
    HR.95L = coxSummary$conf.int[1, "lower .95"],
    HR.95H = coxSummary$conf.int[1, "upper .95"],
    pvalue = coxSummary$coefficients[1, "p"]
  ))
}

# Plot male-specific rs144112368
res.cox <- coxph(
  Surv(DURVISIT, MMSE) ~ rs144112368 + AAO + YEARSEDUC + GBA.APOE4.PHS +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = test1,
  method = "breslow"
)

GG <- ggadjustedcurves(
  res.cox, variable = "rs144112368", data = test1,
  xlab = "Years since onset of PD (years)",
  ylab = "Patients surviving free of global cognitive impairment"
) +
  scale_colour_manual(values = c("#CC79A7", "#0072B2")) +
  geom_hline(aes(yintercept = 0.50), colour = "black", linetype = "dashed") +
  annotate(
    "text", x = 5.0, y = 0.30,
    label = paste(
      "Hazard ratio:", round(as.numeric(result[5, 2]), 3),
      "\n(95% CI:", round(as.numeric(result[5, 3]), 3), "-", round(as.numeric(result[5, 4]), 3),
      ")\nP =", signif(as.numeric(result[5, 5]), 3)
    ),
    colour = "black", size = 5
  ) +
  scale_x_continuous(breaks = seq(0, 12, 2))

####==============================
#### Female-only analysis
female <- subset(independent_snp, type == "MMSE_female")
test1 <- gws_snp[, female$rsID]
test1 <- cbind(test1, gws_snp[12:31])
test1 <- subset(test1, SEX == "F")  # Filter female samples

result <- data.frame()

for (i in colnames(test1[1:nrow(female)])) {
  cox <- coxph(
    Surv(DURVISIT, MMSE) ~ get(i) + AAO + frailty(STUDY_NAME) + YEARSEDUC +
      GBA.APOE4.PHS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    data = test1,
    method = "breslow"
  )
  
  coxSummary <- summary(cox)
  print(i)
  
  result <- rbind(result, cbind(
    id = i,
    HR = coxSummary$conf.int[1, "exp(coef)"],
    HR.95L = coxSummary$conf.int[1, "lower .95"],
    HR.95H = coxSummary$conf.int[1, "upper .95"],
    pvalue = coxSummary$coefficients[1, "p"]
  ))
}

# Plot female-specific rs72616437
res.cox <- coxph(
  Surv(DURVISIT, MMSE) ~ rs72616437 + AAO + YEARSEDUC + GBA.APOE4.PHS +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = test1,
  method = "breslow"
)

GG <- ggadjustedcurves(
  res.cox, variable = "rs72616437", data = test1,
  xlab = "Years since onset of PD (years)",
  ylab = "Patients surviving free of global cognitive impairment"
) +
  scale_colour_manual(values = c("#CC79A7", "#0072B2")) +
  geom_hline(aes(yintercept = 0.50), colour = "black", linetype = "dashed") +
  annotate(
    "text", x = 5.0, y = 0.30,
    label = paste(
      "Hazard ratio:", round(as.numeric(result[5, 2]), 3),
      "\n(95% CI:", round(as.numeric(result[5, 3]), 3), "-", round(as.numeric(result[5, 4]), 3),
      ")\nP =", signif(as.numeric(result[5, 5]), 3)
    ),
    colour = "black", size = 5
  ) +
  scale_x_continuous(breaks = seq(0, 12, 2))
