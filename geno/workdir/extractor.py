#!/usr/bin/python
import os

extractList = open('SNP_CHR_POS.tmp','r')
vcfs = os.listdir('../imputed/')

for line in extractList:
    splitLine = line.split()
    if(splitLine[0]!='SNP'): #BAd coder
        lineSNP = splitLine[0]
        lineChr = int(splitLine[1])
        linePos = int(splitLine[2])
        for j in vcfs:
            i = j.split('.')
            chr = i[2]
            chrNr = int(chr)
            
            chunk = i[3].split('-')
            start = chunk[0]
            end = chunk[1]
            startNr = int(start)
            endNr = int(end)
            
            if(chrNr==lineChr):
                if(linePos <= endNr and linePos >= startNr):
                    wrapper=open(lineSNP+'.wrapper.sh','w')
                    wrapper.write("grep -w "+ lineSNP +" <(zcat ../imputed/"+j+") > ../"+chr+"."+str(linePos)+".snps.on.vcf")
                    wrapper.close()
#            os.system(commandStart + j + ' --snp ' + lineRS + ' --out results/' + lineRS)
#            os.system(commandStart + j + ' --snp ' + lineExRS + ' --out results/' + lineExRS)
#           os.system(commandStart + j + ' --snp ' + lineEx + ' --out results/' + lineEx)
