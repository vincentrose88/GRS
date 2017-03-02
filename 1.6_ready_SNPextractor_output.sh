#!/bin/bash
plink19 --vcf geno/all.SNPs.vcf --r2 --out geno/ldprune
plink19 --vcf geno/all.SNPs.vcf --freq --out geno/maf
