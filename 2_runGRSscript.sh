#!/bin/bash
snpListFile=$1; shift
specFile=$1; shift
idsFile=$1; shift
otherArguments=$@

./create_GRS.R $snpListFile genoFile.noHead newHeader $specFile $idsFile $otherArguments 
