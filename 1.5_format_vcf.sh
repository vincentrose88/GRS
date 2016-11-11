#!/bin/bash
zcat geno/imputed/$(head -1 <(ls geno/imputed/)) | head -n $($(head -1 <(ls geno/imputed/)) | awk '$1=="#CHROM" {print FNR; exit}') > geno/VCFheader
tail -1 VCFheader | sed 's/#//g' | sed 's/\t/"\t"/g' | sed 's/$/"/g' | sed 's/^/"/g' > newHeader

cat geno/VCFheader geno/*snps.on.vcf > geno/all.SNPs.vcf
bcftools annotate -R QUAL,FILTER,INFO,FORMAT/GT,FORMAT/ADS,FORMAT/GP geno/all.SNPs.vcf | grep -v '#' > genoFile.noHead
