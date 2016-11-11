#!/bin/bash
snpListFile=$1; shift

cut -d',' -f3,4 $snpListFile | tr ',' '\t' | sort -nk2 > SNP_CHR.tmp
for i in {1..22}
do
    awk -v "chr=$i" '$2==chr {print $1}' SNP_CHR.tmp > geno/$i.snps
    echo geno/extract.snps.from.vcf.wrapper.sh $i | ./submit_jobarray.py -n extractSNPs$1 -m 4g    
done
 
zcat geno/imputed/$(head -1 <(ls geno/imputed/)) | head -n $($(head -1 <(ls geno/imputed/)) | awk '$1=="#CHROM" {print FNR; exit}') > geno/VCFheader
tail -1 VCFheader | sed 's/#//g' | sed 's/\t/"\t"/g' | sed 's/$/"/g' | sed 's/^/"/g' > newHeader

cat geno/VCFheader geno/*snps.on.vcf > geno/all.SNPs.vcf
bcftools annotate -R QUAL,FILTER,INFO,FORMAT/GT,FORMAT/ADS,FORMAT/GP geno/all.SNPs.vcf | grep -v '#' > genoFile.noHead
