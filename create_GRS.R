#!/usr/bin/Rscript
args <- commandArgs(trailingOnly=T)

snpListFile <- args[1]
genoFile <- args[2]
headerFile <- args[3]
specFile <- args[4]
idsFile <- args[5]
ldFile <- args[6]
outputFile <- args[7]
pvalCutoff <- args[8]
ldCutoff <- args[9]

if(is.na(pvalCutoff)){
    pvalCutoff <- 5e-8
}
if(is.na(outputFile)){
    outputFile <- 'GRS.output'
}
if(is.na(ldCutoff)){
    ldCutoff <- 0.8
}



#Getting genotypes and adding header myself, as read.table is fuzzy
#geno <- read.table('all.SNPs.forR.noHeader',h=F,as.is=T)
#genoHeader <- read.table('newHeader',h=F,as.is=T)
geno <- read.table(genoFile,h=F,as.is=T)
genoHeader <- read.table(headerFile,h=F,as.is=T)
colnames(geno) <- t(genoHeader)
#SNP list from literature
#lit <- read.csv2('../lists/GRS_leadSNPs_full_list.csv',as.is=T)
lit <- read.csv(snpListFile,as.is=T)
#Justin Case code
lit$N <- as.numeric(lit$N)
lit$Effect <- as.numeric(lit$Effect)
lit$EAF <- as.numeric(lit$EAF)
lit$P.value <- as.numeric(lit$P.value)

#Merging with chromosome and rsID (not position, as different build might be an issue)
t <- merge(geno,lit,by.x=c('CHROM','ID'),by.y=c('Chr','SNP'))

#List of particids
#ids <- read.csv2('../lists/birthdate.csv',as.is=T)
#idsFormatted <- paste0(substr(paste0('57x',ids$subjid),1,5),'-',substr(paste0('57x',ids$subjid),6,nchar(paste0('57x',ids$subjid))))
idsdf <- read.table(idsFile,h=F,as.is=T)
idNames <- idsdf$V1
ids <- which(colnames(final) %in% idNames)

final <- t[,colnames(t) %in% c('ID','Effect','EAF','REF','ALT','Effect.Allele','Other.Allele','CHR','POS','Pos','GRS.type','Trait','Locus','Note','N','P.value',idNames)]


#Strand flip function
strandFlip <- function(x){
    y <- NA
    if(x=='A'){
        y <- 'T'}
    else if(x=='G'){
        y <- 'C'
    }else if(x=='C'){
        y <- 'G'
    }else if(x=='T'){
        y <- 'A'
    }else{
        print('not a normal allele, return NA')
    }
    return(y)
}

#Allele checker function
alleleChecker <- function(x,ref='REF',alt='ALT',eff='Effect.Allele',nonEff='Other.Allele'){
    #5 cases:
    #1: REF==eff &&  ALT==nonEff -> Flip!
    #2: REF==nonEff && ALT==eff -> Do nothing
    #3: REF==strandFlip(eff) && ALT==strandFlip(nonEff) -> Flip!
    #4: REF==strandFlip(nonEff) && ALT==strandFlip(eff) -> No allele flip but a strand flip is nice.
    #5: if REF and ALT doesn't match either of the alleles -> (throw error and/or) return NA
    
    case <- NA
    if(x[ref]==x[eff] & x[alt]==x[nonEff]){
        case <- 1
    }else if(x[ref]==x[nonEff] & x[alt]==x[eff]){
        case <- 2
    }else if(x[ref]==strandFlip(x[eff]) & x[alt]==strandFlip(x[nonEff])){
        case <- 3
    }else if(x[ref]==strandFlip(x[nonEff]) & x[alt]==strandFlip(x[eff])){
        case <- 4
    }else if(x[ref]!=x[eff] & x[ref]!=x[nonEff] & x[ref]!=strandFlip(x[eff]) & x[ref]!=strandFlip(x[nonEff])){
        case <- NA
    }else{
        print('No matching alleles found for the following SNP which is discarded:')
        print(x)
        print('------------Check your list of SNPs - site might be triallelic----------------------------')
        case <- NA
    }
    return(case)    
}

final$case <- apply(final,1,alleleChecker)
ids <- which(colnames(final) %in% idNames)

finalBeforeFlip <- final

#Flip those which have cases 1 or 3
for(snp in which(final$case %in% c(1,3))){   
    final[snp,ids] <- -final[snp,ids]+2
    final[snp,'ALT'] <- final[snp,'Effect.Allele']
    final[snp,'REF'] <- final[snp,'Other.Allele']
}
#Remove those with NA cases
final <- final[!is.na(final$case),]

#Flip strand on case 4:
for(snp in which(final$case==2)){   
    final[snp,'ALT'] <- strandFlip(final[snp,'ALT'])
    final[snp,'REF'] <- strandFlip(final[snp,'REF'])
}

#Now the effect allele is always 2 and equal to ALT and the non effect allele is always 0 and equal to REF

#Flip according to effect direction
for(snp in which(final$Effect < 0)){
    final[snp,'Effect'] <- -final[snp,'Effect']
    final[snp,'EAF'] <- 1-final[snp,'EAF']
    final[snp,'REF'] <- final[snp,'Effect.Allele']
    final[snp,'ALT'] <- final[snp,'Other.Allele']
    final[snp,ids] <- -final[snp,ids]+2
}

#Now the POSITIVE effect allele is always 2 and equal to ALT and the NEGATIVE effect allele is always 0 and equal to REF

# remove duplicated with the same trait association
worst <- NULL
for(trait in unique(final$Trait)){
    traitOnly <- final[final$Trait==trait,]
    cand <- traitOnly[traitOnly$ID %in% traitOnly[duplicated(traitOnly$ID),'ID'],]
    if(dim(cand)[1]>1){ #ie, there is duplicated snps for the same trait
        for(snp in unique(cand$ID)){
            snpCand <- cand[cand$ID == snp,]
            worst <- c(worst,
                       row.names(
                           snpCand[-which(snpCand$N==max(snpCand$N,na.rm=T)),]
                           )
                       )
        }
    }
}

if(!is.null(worst)){
    worst <- as.integer(worst)
    final <- final[-worst,]
}   


# LD prune
ld <- read.table(ldFile,h=T,as.is=T)
ld <- ld[ld$R2 > ldCutoff,] #Only look at SNPs with a R2 above 0.8 (default ldCutoff value)
ldSNPsINFO <- final[final$ID %in% c(ld$SNP_A,ld$SNP_B),c('ID','Trait','P.value')]
tmp <- merge(ld,ldSNPsINFO,by.x='SNP_A',by.y='ID',all.x=T)
ld <- merge(tmp,ldSNPsINFO,by.x='SNP_B',by.y='ID',all.x=T)
#two step prune - within same trait and across traits (for combined GRSs)
## Same trait
SNPLD <- NULL
traitLD <- NULL
if(sum(ld$Trait.x==ld$Trait.y)>0){
    traitLD <- c(traitLD,ld$Trait.y)
    cand <- ld[ld$Trait.x==ld$Trait.y,]
    SNPLD <- c(SNPLD,ifelse(cand$P.value.x<cand$P.value.y,cand$SNP_B,cand$SNP_A))   
}
if(!is.null(SNPLD)){
    outLD <- data.frame(SNP=SNPLD,trait=traitLD)
    for(out in outLD){
        final <- final[-which(final$Trait==outLD$trait && final$ID==outLD$SNP),]
    }
}
## Across traits (dont' remove, but note which SNPs are in LD across which traits)
### Remove already pruned SNPs if any
if(sum(ld$Trait.x==ld$Trait.y)>0){
    ld <- ld[-which(ld$Trait.y==ld$Trait.x),]
}
pruneSNPs <- NULL
pruneReason <- NULL
for(row in 1:nrow(ld)){
    i <- ld[row,]
    pruneSNPs <- c(pruneSNPs,ifelse(i$P.value.x<i$P.value.y,i$SNP_A,i$SNP_B))
    pruneReason <- c(pruneReason,paste(i$Trait.y,i$Trait.x))
}
prune <- data.frame(SNP=pruneSNPs,pruneReason=pruneReason,stringsAsFactors=F)
pruneCSNPs <- NULL
pruneCReason <- NULL
#Only adds one, combined, reason for each SNP
for(snp in unique(prune$SNP)){
    pruneCSNPs <- c(pruneCSNPs,snp)
    pruneCReason <- c(pruneCReason,paste(unique(strsplit(paste(prune[prune$SNP==snp,'pruneReason'],collapse=' '),' ')[[1]]),collapse=' '))
}
pruneC <- data.frame(SNP=pruneCSNPs,pruneReason=pruneCReason,stringsAsFactors=F)
finalWithLDINFO <- merge(final,pruneC,by.x='ID',by.y='SNP',all.x=T)

#Formatting clean up (info (16 columns), then genotypes (17 column to the end) and removing the redundant Pos
final <- finalWithLDINFO[,c(which(!colnames(finalWithLDINFO) %in% ids),which(colnames(finalWithLDINFO) %in% ids))]
final <- final[,-which(colnames(final)=='Pos')]
ids <- which(colnames(final) %in% idNames)

GRS <- NULL
for(trait in unique(final$Trait)){
    tmp <- as.data.frame(apply(final[final$Trait==trait,ids],2,sum))
    colnames(tmp) <- trait
    idRows <- rownames(tmp) #Checked manually that it is the same for each iteration
    GRS <- c(GRS,tmp)
}
GRS <- as.data.frame(GRS)
GRS$IID <- idRows
GRS$FID <- idRows

# Testing case - replace with input file:
#GRStypes <- as.data.frame(c('HDL+LDL','FG+FI','HDL+Chol+TC'),stringsAsFactors=F)
#colnames(GRStypes) <- 'V1'

GRStypes <- read.table(specFile,h=F,as.is=T)

finalGRS <- GRS
#Loops through the combined GRS and findes the SNPs that needs to be substracted to avoid double counting due to LD
for(i in 1:nrow(GRStypes)){
    #This is quite 'hacky' but it works
    traitsInGRS <- strsplit(GRStypes[i,],'+',fixed=T)[[1]]
    tmp <- lapply(strsplit(final$pruneReason,' '),function(x){x %in% traitsInGRS})
    pruneRows <- sapply(tmp,sum)
    newGRS <- NULL
    for(x in 1:nrow(GRS)){
        newGRS <- c(newGRS,
                    sum(GRS[x,colnames(GRS) %in% traitsInGRS])
                    - sum(pruneRows*final[,colnames(final)==GRS[x,'IID']])
                    )
    }
    finalGRS <- cbind(finalGRS,newGRS)
    colnames(finalGRS)[ncol(finalGRS)] <- GRStypes[i,]
}

finalGRS <- finalGRS[,c(which(colnames(finalGRS) %in% c('IID','FID')),which(!colnames(finalGRS) %in% c('IID','FID')))]
#Final formatting
write.table(finalGRS,outputFile,row.names=F,quote=F)
