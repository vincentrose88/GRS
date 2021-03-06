#Gene Risk Score automatical turn and twists automatically from vcf file and specs
##You'll need the following:
 1. A csv-file (american with ',' as seperator) (referred as `SNPlist.csv`) with the SNPs from the litterature you want in your GRS. **NB: Has to be in the same format and order as 
example.list**
    * Crucial info is: 
    * Effect allele and non-effect allele
    * SNP name (rs-number)
      * Note that if rs-number is not present in the VCF it will not be extracted. Use SNPextractor, which matches on chr and position, and check names.
    * N in study 
      * crucial when several studies report same SNP or SNPs in high LD
    * P-value 
      * used as N, but also when using a p-value cut-off for what to include in your GRS
    * Chr and position
      * used to match SNPs with extractiong from VCF-file when using SNPextractor
      * **NB:** rs names is used to find genotype data with built-in extraction, and thus should be checked

    * Effect (only for a weighted GRS)
      * Should measured in standard deviations (SD) from a rank-normalized transformation. **Check your paper**
    * Effect Allele Frequency (EAF)

  * Non-crucial, but nice to have for later comparison:
    * Locus or nearest gene
    * Note for notes

2. A spec-file (referred as `GRS.spec`) where you specify which GRS or combinations of GRS you want. 
  * Again format is important - see example (only - and + is allowed, no spaces, one GRS pr. line)
    * If you don't specify a spec-file no combinations of GRS is calculated, instead you will just get one GRS for each trait in your csv-file.

3. A id-file (reffered as `IDs.list`) with the IDs of the individuals in your study. One individual pr. line. 
  * **NB: Has to be exactly the same IDs as in your VCF header. Check beforehand**

4. A already formatted VCF-file to the script (see example) or path to the imputed genotypes (in vcf-format and with names as `1.vcf.gz` for chr1) and use script
`1_create_geno_from_vcf.sh`

5. A good relation with me, Vincent, if the case of trouble-shooting.

6. You migth want to update your path to have the same programs as me, so edit your `.bashrc` file by:
 * `nano ~/.bashrc`
 * add `export PATH=$PATH:/home/fng514/bin:/home/cxt155/bin:/home/cxt155/samtools-bcftools-htslib-1.0_x64-linux`

##Running the program

1. Start by running the setup:

  * `./0_setup.sh path-to-impute-folder-with-vcfs`

**NB** I strongly recommend to use my other script: SNPextractor (https://github.com/vincentrose88/SNPextractor) to extract variants with the -g flag for direct usage with GRS**

2. (Optinal) Create the formatted vcf-file for R by (using the names of the SNPs to search through vcfs):

  * `./1_create_geno_from_vcf.sh SNPlist.csv`

  * When finished (quits after submitting two scripts to the grid-engine (check with qstat -f)) two files are created in the folder: `genoFile.noHead` and `newHeader` along some gridengine outputs

  * **Recommended alternative:** Use the SNPextractor to extract SNPs more effectivly (from positions instead of names) and use the -g flag to get GRS ready output (`genoFile.noHead` and `newHeader`) and use the 1.6_ready_SNPextrator_output.sh script (remember to mv the all.SNPs.vcf into your geno/ folder) to make everything ready for the GRS
   * SNP-names might be different, as not all SNPs have names in the vcfs, and thus a dummy name (chr:pos) is used insead

3. Then run the main script with the arguments as follows. If you have created your own geno file, header and/or ld-file, then these can be specified with --geno, --head and --ldfile, repectively.

  * `./2_create_GRS.R --snp SNPlist.csv --ids IDs.list --spec GRS.spec --out example --pcut 1e-7 --ldcut 0.5 --SNPout SNPs.flipped` 

    * Besides the two first arguments (--snp and --ids) everything else is optinal:

      * `--SNPextractor` Flag that tells the GRS to merge on positions when using a SNPextractor output.
    
      * `--spec` Name of file containting specifications for how to combine GRSs. See example `GRS.spec`

      * `--out` (default='GRS.output'). Name of outputfile.

      * `--pcut` (default=5e-8 (GWAS significant)). Threshold for inclusion of SNPs in GRS.

      * `--ldcut` (default 0.8). Threshold for pruning SNPs with a R2 above cutoff.

      * `--SNPout` (default NA, ie. not used). Name of file where the genotypes for all SNPs will be saved *AFTER* being flipped to match the litterature. REF is the non-effect allele, coded as 0, and ALT is the Effect allele, coded as 2. Used for other analysis.

4. Output files will be the GRS, weighted and unweighted for each trait as well as all combinations specified in the .spec file along with `Litterature_SNPs_not_extracted` listing the SNPs not succesfully extracted


By Vincent Appel
