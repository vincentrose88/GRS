#!/usr/bin/Rscript
args <- commandArgs(trailingOnly=T)

snpListFile <- args[1]
genoFile <- args[2]
headerFile <- args[3]
idsFile <- args[4]
ldFile <- args[5]
specFile <- args[6]
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

print(args)
print('loading in files')


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

#Hacking the god damn exm-snps:
exmGeno <- geno[grep('exm-',geno[,3]),3]
exGenoFix <- substr(exmGeno,5,nchar(exmGeno))
geno[grep('exm-',geno[,3]),3] <- exGenoFix



#Merging with chromosome and rsID (not position, as different build might be an issue)
t <- merge(geno,lit,by.x=c('CHROM','ID'),by.y=c('Chr','SNP'))

#List of particids
#ids <- read.csv2('../lists/birthdate.csv',as.is=T)
#idsFormatted <- paste0(substr(paste0('57x',ids$subjid),1,5),'-',substr(paste0('57x',ids$subjid),6,nchar(paste0('57x',ids$subjid))))
idsdf <- read.table(idsFile,h=F,as.is=T)
idNames <- idsdf$V1

ids <- which(colnames(t) %in% idNames)

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
        print(x[c('ID','POS',ref,alt,eff,nonEff)])
        print('------------Check your list of SNPs - site might be triallelic----------------------------')
        case <- NA
    }
    return(case)    
}

final$case <- apply(final,1,alleleChecker)
ids <- which(colnames(final) %in% idNames)

finalBeforeFlip <- final
print('flipping alleles...')

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

print('removing duplicates...')
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

print('ld pruning step...')
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
    traitLD <- c(traitLD,ld[ld$Trait.x==ld$Trait.y,'Trait.y'])
    cand <- ld[ld$Trait.x==ld$Trait.y,]
    SNPLD <- c(SNPLD,ifelse(cand$P.value.x<cand$P.value.y,cand$SNP_B,cand$SNP_A))   
}
if(!is.null(SNPLD)){
    outLD <- data.frame(SNP=SNPLD,trait=traitLD,stringsAsFactors=F)
    for(out in 1:nrow(outLD)){
        final <- final[-which(final$Trait==outLD[out,'trait'] & final$ID==outLD[out,'SNP']),]
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

#Filter on note - special case
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


#If no combination GRS is requested, ie GRS.spec is NA, then don't try anything funny.
if(!is.na(specFile)){

#Read in GRS.spec
GRStypes <- read.table(specFile,h=F,as.is=T)


#Format GRS spec to traits and sign (+ or -)

finalGRS <- GRS
#Loops through the combined GRS and findes the SNPs that needs to be substracted to avoid double counting due to LD
for(i in 1:nrow(GRStypes)){
    #This is quite 'hacky' but it works
    traitsInGRS <- strsplit(GRStypes[i,],'+',fixed=T)[[1]]
    if(any(grepl('-',traitsInGRS,fixed=T))){ #Finding and splitting off traits with minus for proper formatting
        normalTrait <-  traitsInGRS[-which(grepl('-',traitsInGRS,fixed=T))] #Normal trait is added together
        notNormalTrait <- traitsInGRS[which(grepl('-',traitsInGRS,fixed=T))] #Special traits is substracted from normal Traits (see example below)
        specialTrait <- NULL
        # Three possibilites:
        ## 1. only one minus trait (fx HDL)
        ## 2. two minus traits next to eachother (fx LDL-HDL-TG)
        ## 3. two minus traits far from eachother (fx LDL-HDL+Chol-TG)
        #First case is easy, 2nd og 3rd not so.
        ##Case 1:
        if(length(notNormalTrait)==1){
            minusPos <- regexpr('-',notNormalTrait,fixed=T)[1]
            normalTrait <- c(normalTrait,substr(notNormalTrait,1,minusPos-1))
            leftOver <- substr(notNormalTrait,minusPos+1,nchar(notNormalTrait))
            if(grepl('-',leftOver,fixed=T)){
                ##Case 2:
                while(grepl('-',leftOver,fixed=T)){ #Keeps adding traits to special traits from leftOver as long as there are minuses to find
                    #Just as Case 1, but replacing normalTrait with specialTrait
                    #and specialTrait with Leftover
                    #and loops with while
                    minusPos <- regexpr('-',leftOver,fixed=T)[1]
                    specialTrait <- c(specialTrait,substr(leftOver,1,minusPos-1)) #Now adding to 
                    leftOver <- substr(leftOver,minusPos+1,nchar(leftOver))
                }
                specialTrait <- c(specialTrait,leftOver) #Collect to last special trait from leftOver
            }else{
                ##Case 1: Only one trait in leftOver which is special
                specialTrait <- leftOver
            }
        }else if(length(notNormalTrait)>1){
            ## Case 3 - the longest one - each element can be a case 1 or 2, thus we just do the same as above, and collect smartly:
            caseThree <- notNormalTrait
            for(j in 1:length(caseThree)){
                notNormalTrait <- caseThree[j]
                minusPos <- regexpr('-',notNormalTrait,fixed=T)[1]
                normalTrait <- c(normalTrait,substr(notNormalTrait,1,minusPos-1))
                leftOver <- substr(notNormalTrait,minusPos+1,nchar(notNormalTrait))
                if(grepl('-',leftOver,fixed=T)){
                    ##Case 2:
                    while(grepl('-',leftOver,fixed=T)){ #Keeps adding traits to special traits from leftOver as long as there are minuses to find
                    #Just as Case 1, but replacing normalTrait with specialTrait
                    #and specialTrait with Leftover
                    #and loops with while
                        minusPos <- regexpr('-',leftOver,fixed=T)[1]
                        specialTrait <- c(specialTrait,substr(leftOver,1,minusPos-1)) #Now adding to 
                        leftOver <- substr(leftOver,minusPos+1,nchar(leftOver))
                    }
                    specialTrait <- c(specialTrait,leftOver) #Collect to last special trait from leftOver
                }else{
                    ##Case 1: Only one trait in leftOver which is special
                    specialTrait <- c(specialTrait,leftOver)
                }
            }
        }
    }else{ #No minus/special traits
        normalTrait <- traitsInGRS
        specialTrait <- NULL
    }

    tmp <- lapply(strsplit(final$pruneReason,' '),function(x){x %in% c(normalTrait,specialTrait)})
    pruneRows <- sapply(tmp,sum) #Using sum to substract the same SNP away the right amount of times
    newGRS <- NULL
    if(is.null(specialTrait)){ #No special/minus traits
        for(x in 1:nrow(GRS)){
            newGRS <- c(newGRS,
                        sum(GRS[x,colnames(GRS) %in% normalTrait])
                        - sum(pruneRows*final[,colnames(final)==GRS[x,'IID']])
                        )
        }
    }else{
        for(x in 1:nrow(GRS)){
            newGRS <- c(newGRS,
                        sum(GRS[x,colnames(GRS) %in% normalTrait])
                        - sum(GRS[x,colnames(GRS) %in% specialTrait])
                        - sum(pruneRows*final[,colnames(final)==GRS[x,'IID']])
                        )
        }
    }
    finalGRS <- cbind(finalGRS,newGRS)
    colnames(finalGRS)[ncol(finalGRS)] <- GRStypes[i,]
}

}else{
	finalGRS <- GRS
}


print(paste('final clean up and creating resulting GRS file as',outputFile))

finalGRS <- finalGRS[,c(which(colnames(finalGRS) %in% c('IID','FID')),which(!colnames(finalGRS) %in% c('IID','FID')))]
#Final formatting
write.table(finalGRS,outputFile,row.names=F,quote=F)
