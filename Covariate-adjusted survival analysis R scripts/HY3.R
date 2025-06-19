# HY3-related survival analysis using selected independent SNPs (combined + female stratified)

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
setwd("D:\\lyyy\\chr_X\\cox\\2025_3_2\\cox_summary\\cov_adjust_cox\\HY3")

# Load significant SNPs (from FUMA)
independent_snp <- read_xlsx("D:\\lyyy\\chr_X\\cox\\2025_3_2\\cox_summary\\HY3_cox\\fuma_output\\independent_snp.xlsx")
colnames(independent_snp)[10] <- "type"  # Assign 'type' column

# Load genotype + covariates matrix
gws_snp <- read.table("independent_snp_cor_cox.txt")
gws_snp[, 1:2][gws_snp[, 1:2] == 2] <- 1  # Recode: 2 → 1 for additive model

### ===============================
### Combined analysis (rs3128076)

combined <- subset(independent_snp, type == "HY3_combined")
test1 <- gws_snp[, combined$rsID, drop = FALSE]
test1 <- cbind(test1, gws_snp[3:20])  # Add covariates

result <- data.frame()

for (i in colnames(test1[1:nrow(combined)])) {
  cox <- coxph(
    Surv(DURVISIT, HY) ~ get(i) + AAO + SEX + frailty(STUDY_NAME) +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
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

# KM curve for rs3128076
res.cox <- coxph(
  Surv(DURVISIT, HY) ~ rs3128076 + AAO + SEX +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = test1,
  method = "breslow"
)

GG <- ggadjustedcurves(
  res.cox, variable = "rs3128076", data = test1,
  xlab = "Years since onset of PD (years)",
  ylab = "Patients surviving free of HY stage 3"
) +
  scale_colour_manual(values = c("#CC79A7", "#0072B2")) +
  geom_hline(aes(yintercept = 0.50), colour = "black", linetype = "dashed") +
  annotate(
    "text", x = 5.0, y = 0.30,
    label = paste(
      "Hazard ratio:", round(as.numeric(result[1, 2]), 3),
      "\n(95% CI:", round(as.numeric(result[1, 3]), 3), "-", round(as.numeric(result[1, 4]), 3),
      ")\nP =", signif(as.numeric(result[1, 5]), 3)
    ),
    colour = "black", size = 5
  ) +
  scale_x_continuous(breaks = seq(0, 12, 2))


### ===============================
### Female-only analysis (rs111708875)

female <- subset(independent_snp, type == "HY3_female")
test1 <- gws_snp[, female$rsID]
test1 <- cbind(test1, gws_snp[3:20])
test1 <- subset(test1, SEX == "F")  # Filter for female participants

result <- data.frame()

for (i in colnames(test1[1:nrow(female)])) {
  cox <- coxph(
    Surv(DURVISIT, HY) ~ get(i) + AAO + frailty(STUDY_NAME) +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
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

# KM curve for rs111708875
res.cox <- coxph(
  Surv(DURVISIT, HY) ~ rs111708875 + AAO +
    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = test1,
  method = "breslow"
)

GG <- ggadjustedcurves(
  res.cox, variable = "rs111708875", data = test1,
  xlab = "Years since onset of PD (years)",
  ylab = "Patients surviving free of HY stage 3"
) +
  scale_colour_manual(values = c("#CC79A7", "#0072B2")) +
  geom_hline(aes(yintercept = 0.50), colour = "black", linetype = "dashed") +
  annotate(
    "text", x = 5.0, y = 0.30,
    label = paste(
      "Hazard ratio:", round(as.numeric(result[1, 2]), 3),
      "\n(95% CI:", round(as.numeric(result[1, 3]), 3), "-", round(as.numeric(result[1, 4]), 3),
      ")\nP =", signif(as.numeric(result[1, 5]), 3)
    ),
    colour = "black", size = 5
  ) +
  scale_x_continuous(breaks = seq(0, 12, 2))
