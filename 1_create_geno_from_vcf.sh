#!/bin/bash
snpListFile=$1; shift

cat $snpListFile | tr '\r' '\n' | cut -d',' -f2,3 | tr ',' '\t' | sort -nk2 > SNP_CHR.tmp
for i in {1..22}
do
    awk -v "chr=$i" '$2==chr {print $1}' SNP_CHR.tmp > geno/$i.snps
done

for i in {1..22}
do
    echo geno/extract.snps.from.vcf.wrapper.sh $i
done | ./submit_jobarray.py -n extractSNPs -m 4g    
 
echo ./1.5_format_vcf.sh | ./submit_jobarray.py -w extractSNPs -n formatVCF
