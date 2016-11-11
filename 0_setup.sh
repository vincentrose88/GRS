#!/bin/bash
imputedFolder=$1; shift

mkdir -p geno/imputed
ln -s $imputedFolder/* geno/imputed/

chmod +x *.sh
chmod +x *.R
chmod +x *.py
