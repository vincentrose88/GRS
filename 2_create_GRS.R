#!/usr/bin/Rscript
args <- commandArgs(trailingOnly=T)

snpListFile <- NA
genoFile <- 'genoFile.noHead'
headerFile <- 'newHeader'
specFile <- NA
idsFile <- NA
ldFile <- 'geno/ldprune.ld'
outputFile <- 'GRS.output'
pvalCutoff <- 5e-8
ldCutoff <- 0.8
SNPout <- NA
SNPextractor <- FALSE

##Args with flags
for(i in 1:length(args)){
    print(args[i])
    if(args[i]=='--snp'){
        snpListFile <- args[i+1]
    }else if(args[i]=='--geno'){
        genoFile <- args[i+1]
    }else if(args[i]=='--head'){
        headerFile <- args[i+1]
    }else if(args[i]=='--spec'){
        specFile <- args[i+1]
    }else if(args[i]=='--ids'){
        idsFile <- args[i+1]
    }else if(args[i]=='--ldfile'){
        ldFile <- args[i+1]
    }else if(args[i]=='--out'){
        outputFile <- args[i+1]
    }else if(args[i]=='--pcut'){
        pvalCutoff <- args[i+1]
    }else if(args[i]=='--ldcut'){
        ldCutoff <- args[i+1]
    }else if(args[i]=='--SNPout'){
        SNPout <- args[i+1]
    }else if(args[i]=='--SNPextractor'){
        SNPextractor <- TRUE
    }else if(grepl('--',args[i])){
        print('flag not recognized or missing. Exiting..')
        q()
    }
}

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


#What did't get extracted?
write.table(lit[-which(lit$SNP %in% geno$ID),],'Litterature_SNPs_not_extracted',row.names=F,quote=F)


maf <- read.table('geno/maf.frq',h=T,as.is=T)
maf <- maf[,-c(3:4,6)]

#Removing tri-allelic sites from the MAF support file to avoid merging articafts. They are still kept in the geno file.
if(length(unique(maf$SNP)) < length(maf$SNP)){
    worstMafCand <- NULL
    for(triSNP in maf$SNP[duplicated(maf$SNP)]){

        mafCand <-  maf[maf$SNP == triSNP,]
        litEAF  <- lit[lit$SNP ==triSNP,'EAF']
        worstMafCand <- c(worstMafCand,
                          as.integer(
                              row.names(
                                  mafCand[which(
                                      abs(mafCand$MAF-litEAF)==max(abs(mafCand$MAF-litEAF)))
                                          #Saving only the SNP which has a MAF closest to that of the EAF
                                          ,]
                                  )
                              )
                          )
    }
    maf <- maf[-worstMafCand,]
}

if(SNPextractor){ #merging on position when extracting with SNPextractor and keeping the names from SNPextractor (chr:pos dummy names for missing names)
    tmp <- merge(geno,lit,by.x=c('CHROM','POS'),by.y=c('Chr','Pos'))
    t <- merge(tmp,maf,by.x=c('CHROM','ID'),by.y=c('CHR','SNP')) #The MAF-file is created from the all.SNPs.vcf, which also have SNPextractor names, and thus this merge should work fine
}else{
#Merging with chromosome and rsID (not position, as different build might be an issue)
    tmp <- merge(geno,lit,by.x=c('CHROM','ID'),by.y=c('Chr','SNP'))
    t <- merge(tmp,maf,by.x=c('CHROM','ID'),by.y=c('CHR','SNP'))
}

#List of particids
#ids <- read.csv2('../lists/birthdate.csv',as.is=T)
#idsFormatted <- paste0(substr(paste0('57x',ids$subjid),1,5),'-',substr(paste0('57x',ids$subjid),6,nchar(paste0('57x',ids$subjid))))
idsdf <- read.table(idsFile,h=F,as.is=T)
idNames <- idsdf$V1

#Merging on pos (SNPextractor=TRUE) or name (SNPextractor=FALSE)
ids <- which(colnames(t) %in% idNames)
if(SNPextractor){
    final <- t[,colnames(t) %in% c('ID','SNP','Effect','EAF','REF','ALT','Effect.Allele','Other.Allele','MAF','CHROM','POS','GRS.type','Trait','Locus','Note','N','P.value',idNames)]
}else{
    final <- t[,colnames(t) %in% c('ID','Effect','EAF','REF','ALT','Effect.Allele','Other.Allele','MAF','CHROM','POS','Pos','GRS.type','Trait','Locus','Note','N','P.value',idNames)]
}
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
    #Special case when: (REF/ALT || eff/nonEff)==(A/T || G/C). Doublecheck with Minor allele==ref/alt when proceeding:
    #1: REF==eff &&  ALT==nonEff -> Flip!
    #2: REF==nonEff && ALT==eff -> Do nothing
    #3: REF==strandFlip(eff) && ALT==strandFlip(nonEff) -> Flip! (allele and strand)
    #4: REF==strandFlip(nonEff) && ALT==strandFlip(eff) -> No allele flip but a strand flip is nice.
    #NA-case: if REF and ALT doesn't match either of the alleles -> (throw error and/or) return NA
    
    case <- NA
    if(((x[ref]=='A' & x[alt]=='T') | (x[ref]=='T' & x[alt]=='A')
        ) | (
        (x[ref]=='G' & x[alt]=='C') | (x[ref]=='C' & x[alt]=='G')
        )
       ){
        if(is.na(x['EAF'])){
            print(x[c('ID','POS',ref,alt,eff,nonEff,'MAF','EAF','Effect','N')])
            print('SNP is a tricky SNP, and does not have any EAF info from the .csv file, thus the script does not know how to flip it. Add EAF information, and run again. Skipping for now')
            case <- NA
        }else{
        
        ## Special case - split into EAF above or below 50%
        ### New variable MAF decide if MAF is in the same 'direction' as EAF
        ##4 cases in each catagory (shown for when EAF<50%)
        #1 maf flip: MafDir = FALSE, eff=ref. ACTION: Allele FLIP (case 1)
        #2 strand flip: MafDir = TRUE, eff=ref. Action: Strand FLIP (case 4)
        #3 Everything checks out: MafDir = TRUE, eff=alt: Do nothing (case 2)
        #4 maf and strand flip: MafDir = FALSE, eff=alt: Allele and strand FLIP (case 3)
        ##When EAF > 50% the alt needs to be the eff, Then case 2<->1 and case 3<->4:
        
        #1 maf flip: MafDir = FALSE, eff=alt. Do Nothing (case 1->2)
        #2 strand flip: MafDir = TRUE, eff=alt. Allele and Strand FLIP (case 4->3)
        #3 Everything checks out: MafDir = TRUE, eff=ref: Allele FLIP (case 2->1)
        #4 maf and strand flip: MafDir = FALSE, eff=ref: Strand FLIP (case 3->4)
        
        if(x['EAF']<0.5){ 
            mafDir <- x['MAF']<0.5
            if(!mafDir & (x[eff]==x[ref] & x[nonEff] == x[alt])){
                case <- 2
            }else if(mafDir & (x[eff]==x[ref] & x[nonEff] == x[alt])){
                case <- 3
            }else if(mafDir & (x[eff]==x[alt] & x[nonEff] == x[ref])){
                case <- 1
            }else if(!mafDir & (x[eff]==x[alt] & x[nonEff] == x[ref])){
                case <- 4
            }else if(x[ref]!=x[eff] & x[ref]!=x[nonEff] & x[ref]!=strandFlip(x[eff]) & x[ref]!=strandFlip(x[nonEff])){
                #Old catcher of mismatches of alleles (regardless of MAF/EAF mismatch)
                case <- NA
            }else{ # Nothing, really NOTHING fits.
                print('Special A/T or G/C allele, but nothing fits')
                print('No matching alleles found for the following SNP which is discarded:')
                print(x[c('ID','POS',ref,alt,eff,nonEff,'MAF','EAF','N','Effect')])
                print('------------Check your list of SNPs - site might be triallelic----------------------------')
                case <- NA
            }
        }else if(x['EAF']>=0.5){
            mafDir <- x['MAF']>=0.5
            if(!mafDir & (x[eff]==x[alt] & x[nonEff] == x[ref])){
                case <- 1
            }else if(mafDir & (x[eff]==x[alt] & x[nonEff] == x[ref])){
                case <- 4
            }else if(mafDir & (x[eff]==x[ref] & x[nonEff] == x[alt])){
                case <- 2
            }else if(!mafDir & (x[eff]==x[ref] & x[nonEff] == x[alt])){
                case <- 3
            }else if(x[ref]!=x[eff] & x[ref]!=x[nonEff] & x[ref]!=strandFlip(x[eff]) & x[ref]!=strandFlip(x[nonEff])){
                #Old catcher of mismatches of alleles (regardless of MAF/EAF mismatch)
                case <- NA
            }else{ # Nothing, really NOTHING fits.
                print('Special A/T or G/C allele, but nothing fits')
                print('No matching alleles found for the following SNP which is discarded:')
                print(x[c('ID','POS',ref,alt,eff,nonEff,'MAF','EAF','N','Effect')])
                print('------------Check your list of SNPs - site might be triallelic----------------------------')
                case <- NA
            }
        }
    }
        #Back to 'normal' SNPs
    }else if(x[ref]==x[eff] & x[alt]==x[nonEff]){
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
        print(x[c('ID','POS',ref,alt,eff,nonEff,'EAF','N','Effect')])
        print('------------Check your list of SNPs - site might be triallelic----------------------------')
        case <- NA
    }
    return(case)    
}




final$case <- apply(final,1,alleleChecker)
ids <- which(colnames(final) %in% idNames)

finalBeforeFlip <- final
print('flipping alleles...')

#Flip those which have cases 1 or 3 (for case 3 this automatically fixes the strandflip)
for(snp in which(final$case %in% c(1,3))){   
    final[snp,ids] <- -final[snp,ids]+2
    final[snp,'ALT'] <- final[snp,'Effect.Allele']
    final[snp,'REF'] <- final[snp,'Other.Allele']
    final[snp,'MAF'] <- 1-final[snp,'MAF']
}
#Remove those with NA cases
final <- final[!is.na(final$case),]

#Flip strand on case 4:
for(snp in which(final$case==4)){   
    final[snp,'ALT'] <- strandFlip(final[snp,'ALT'])
    final[snp,'REF'] <- strandFlip(final[snp,'REF'])
}

#Now the effect allele is always 2 and equal to ALT and the non effect allele is always 0 and equal to REF

#Flip according to effect direction
for(snp in which(final$Effect < 0)){
    final[snp,'Effect'] <- -final[snp,'Effect']
    final[snp,'MAF'] <- 1-final[snp,'MAF']
    final[snp,'EAF'] <- 1-final[snp,'EAF']
    final[snp,'REF'] <- final[snp,'Effect.Allele']
    final[snp,'ALT'] <- final[snp,'Other.Allele']
    final[snp,ids] <- -final[snp,ids]+2
}
colnames(final)[colnames(final)=='MAF'] <- 'EAF.ownData'

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

#Formatting clean up (info (16 columns), then genotypes (17 column to the end)
final <- finalWithLDINFO[,c(which(!colnames(finalWithLDINFO) %in% ids),which(colnames(finalWithLDINFO) %in% ids))]
ids <- which(colnames(final) %in% idNames)

#If SNPout is not NA, save the genofile with flipped SNPs as the SNPout name
if(!is.na(SNPout)){
    print(paste('Saving genotypes for SNPs AFTER flips for other analysis in',SNPout))
    write.table(final[,c(which(colnames(final) %in% c('ID','SNP','CHROM','POS','REF','ALT','Trait','EAF.ownData')),ids)],SNPout,row.names=F,quote=F)
}


#Adding SNPs with weights
wFinal <- final[,ids]*final[,'Effect']
wFinal$ID <- paste(final$ID,'weighted',sep='_')
wFinal$Trait <- paste(final$Trait,'weighted',sep='_')

#copying the rest of the final data-frame
for(col in colnames(final)[!colnames(final) %in% colnames(wFinal)]){
    wFinal[,col] <- final[,col]
}
#Just some ordering to have the same format
final <- final[,order(colnames(final),decreasing=T)]
wFinal <- wFinal[,order(colnames(wFinal),decreasing=T)]

#Binding everything together (IDs are different for SNPs, so should be fine
completeFinal <- rbind(final,wFinal)
ids <- which(colnames(completeFinal) %in% idNames)


print(paste('calculating GRS as specified by',specFile))

GRS <- NULL
for(trait in unique(completeFinal$Trait)){
    inGRS <- completeFinal[completeFinal$Trait==trait,]
    tmp <- as.data.frame(apply(inGRS[,ids],2,sum))
    colnames(tmp) <- trait
    #Normalzing the weighted GRS
    if(grepl('weighted',trait)){
        normalizing <- nrow(inGRS)/sum(inGRS$Effect,na.rm=T)
        tmp <- tmp*normalizing
    }
    idRows <- rownames(tmp) #Checked manually that it is the same for each iteration
    GRS <- c(GRS,tmp)
}
GRS <- as.data.frame(GRS)
GRS$IID <- idRows
GRS$FID <- idRows


#Read in GRS.spec
if(!is.na(specFile)){
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
