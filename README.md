#Gene Risk Score automatical turn and twists automatically from vcf file and specs
##You'll need the following:
 1. A csv-file (american with ',' as seperator) (referred as `SNPlist.csv`) with the SNPs from the litterature you want in your GRS. **NB: Has to be in the same format and order as 
example.list**
    * Crucial info is: 
    * Effect allele and non-effect allele
    * SNP name (rs-number)
    * N in study 
      * crucial when several studies report same SNP or SNPs in high LD
    * P-value 
      * used as N, but also when using a p-value cut-off for what to include in your GRS
    * Chr and position
      * used to extract genotypes from VCF-file
      * **NB:** rs names is used to find genotype data, and thus should be checked

    * Effect (only for a weighted GRS)
      * measured in standard deviations (SD) from a rank-normalized transformation. **Check your paper**
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

##Running the program

1. Start by running the setup:

  * `./0_setup.sh path-to-impute-folder-with-vcfs`

2. (Optinal) Create the formatted vcf-file for R by:

  * `./1_create_geno_from_vcf.sh SNPlist.csv`

  * When finished (quits after submitting two scripts to the grid-engine (check with qstat -f)) two files are created in the folder: `genoFile.noHead` and `newHeader` along some gridengine outputs

3. Then run the main script with the arguments as follows. If you have created your own geno file, header and/or ld-file, then these can be specified with --geno, --head and --ldfile, repectively.

  * `./2_create_GRS.R --snp SNPlist.csv --ids IDs.list --spec GRS.spec --out example --pcut 1e-7 --ldcut 0.5`

    * Besides the two first arguments (--snp and --ids) everything else is optinal:
    
      * `--spec` Name of file containting specifications for how to combine GRSs. See example `GRS.spec`

      * `--out` (default='GRS.output'). Name of outputfile.

      * `--pcut` (default=5e-8 (GWAS significant)). Threshold for inclusion of SNPs in GRS.

      * `--ldcut` (default 0.8). Threshold for pruning SNPs with a R2 above cutoff.

By Vincent Appel
