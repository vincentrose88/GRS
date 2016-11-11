#!/bin/bash
chr=$1
grep -wFf $chr.snps <(zcat imputed/$chr.vcf.gz) > $chr.snps.on.vcf
