#!/bin/bash

set -e

# Remove SNPs in PARs regions, apply Hardy-Weinberg equilibrium filter, 
# and filter individuals/SNPs based on missing genotype rates

plink --bfile chr23 \                      
      --exclude need_to_exclude \             # Exclude specified SNPs
      --geno 0.05 \                           # Filter SNPs with >5% missing data
      --make-bed \                          
      --mind 0.1 \                            # Filter individuals with >10% missing data
      --out pre_imputation11       