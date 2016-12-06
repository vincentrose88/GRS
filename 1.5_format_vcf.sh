#!/bin/bash
zcat geno/imputed/$(head -1 <(ls geno/imputed/)) | awk '{print $0}; $1=="#CHROM" {exit}' > geno/VCFheader
tail -1 geno/VCFheader | sed 's/#//g' | sed 's/\t/"\t"/g' | sed 's/$/"/g' | sed 's/^/"/g' > newHeader

cat geno/VCFheader geno/*snps.on.vcf > geno/all.SNPs.vcf
bcftools annotate -R QUAL,FILTER,INFO,FORMAT/GT,FORMAT/ADS,FORMAT/GP geno/all.SNPs.vcf | grep -v '#' > genoFile.noHead
plink19 --vcf geno/all.SNPs.vcf --r2 --out geno/ldprune
plink19 --vcf geno/all.SNPs.vcf --freq --out geno/maf
