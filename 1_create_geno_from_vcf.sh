#!/bin/bash
snpListFile=$1; shift
imputedFolder=$1; shift

mkdir -p geno/imputed
ln -s $imputedFolder/* geno/imputed/

cut -d',' -f3,4 $snpListFile | tr ',' '\t' | sort -nk2 > SNP_CHR.tmp
for i in {1..22}
do
    awk -v "chr=$i" '$2==chr {print $1}' SNP_CHR.tmp > geno/$i.snps
done


