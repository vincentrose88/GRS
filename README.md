#Gene Risk Score automatical turn and twists automatically from vcf file and specs
##You'll need the following:
A csv-file (american with ',' as seperator) (referred as `SNPlist.csv`) with the SNPs from the litterature you want in your GRS. **NB: Has to be in the same format as example.list**

A spec-file (referred as `GRS.spec`) where you specify which GRS or combinations of GRS you want. Again format is important - see example (only - and + is allowed, no spaces, one GRS pr. line)

A id-file (reffered as `IDs.list`) with the IDs of the individuals in your study. **NB: Has to be exactly as in your VCF header. Check beforehand**

A already formatted VCF-file to the script (see example) or path to the imputed genotypes (in vcf-format and with names as `1.vcf.gz` for chr1) and use script
`1_create_geno_from_vcf.sh`

A good standing with me, Vincent, if the case of trouble-shooting.

##Running the program

Start by running the setup:

`./0_setup.sh path-to-impute-folder-with-vcfs`

(Optinal) Create the formatted vcf-file for R by:

`./1_create_geno_from_vcf.sh SNPlist.csv`

When finished (quits after submitting two scripts to the grid-engine (check with qstat -f)) two files are created in the folder: `genoFile.noHead` and `newHeader` along some gridengine outputs

Then just run the wrapper for the main script with the arguments as follows (if you didn't automatically create the formatted VCF-file then you must manually create genoFile.noHead and newHeader with those name precisely):

`./runGRSscript.sh SNPlist.csv GRS.spec IDs.list (pval-cutoff outputName)`

At the end you can add the optinal arguments (in brackets) in the order:

`outputName` (default='GRS.output'). Name of outputfile.

`pval-cutoff` (default=5e-8 (GWAS significant)). Threshold for inclusion of SNPs in GRS.

`LD-cutoff` (default 0.8). Threshold for pruning SNPs with a R2 above cutoff.

By Vincent Appel