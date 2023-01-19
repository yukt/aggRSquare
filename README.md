# aggRSquare - A Tool for Calculating Imputation Accuracy
## Installation
Prerequisites: [cget](http://cget.readthedocs.io/en/latest/src/intro.html#installing-cget), [cmake](https://cmake.org/install/).
```
git clone https://github.com/yukt/aggRSquare.git
cd aggRSquare
bash install.sh
```

## Usage
```
aggRSquare -v [Validation.vcf.gz]         [Required] Input Validation File
           -i [Imputation.vcf.gz]         [Required] Input Imputation File
           -o [OutputPrefix]              [Required] Output Prefix
           --validationFormat [GT/DS/GP]  Genotype info format (Default: GT)
           --imputationFormat [GT/DS/GP]  Genotype info format (Default: DS)
           --AF [AlleleFrequency File]    Tab-delimited text file contains SNP and Allele Frequency.
           --bins [Bins File]             A text file contains MAF cutoffs for bins
           --detail                       If ON, an addtional SNP-wise r2 will be output (Default: OFF)
           --help                         If ON, the usage documentation will be displayed.
```
## Input Files
### Bins File
The bins file should contain one MAF cutoff value per line.
The default cutoffs are
```
0.0000
0.0005
0.001
0.002
0.005
0.010
0.015
0.020
0.035
0.050
0.100
0.200
0.300
0.400
0.500
```
which will generate MAF bins `(0,0.005]`,`(0.005,0.001]`, ... , `(0.4,0.5]`.

### Allele Frequency File
The allele frequency file should contain two **tab-delimited** columns including SNP and the corresponding allele frequency. No header is needed, and any line starting with "#" will be ignored.

The SNP should be in the format of CHROM:POS:REF:ALT, and the order of SNPs should be exactly the same as that in the input validation file. The allele frequency should be float numbers, and any other values will be treated as missingness and igored in the analysis.

Example:
```
#SNP	AF
chr1:925875:A:G	NA
chr1:925881:G:A	0.002625
```


## Output Files
Aggregate r2 will be saved in `[OutputPrefix].aggRSquare`, which contains the following columns (tab-delimited).
* Bin.Aggregated.by.MAF : MAF bin in the format of (min, max]
* Average.MAF           : the average MAF of all the variants in the bin
* #Variants             : the number of variants included in the bin
* Imputation.R2         : the aggregate r2 (squared Pearson correlation) 
* Gold.MAF              : the average MAF calculated from genotype data in the input validation file
* Imputed.MAF           : the average MAF calculated from genotype data in the input imputation file

Note: If the allele frequency file (`--AF`) is not detected, MAF will be calculated from the input validation file, and the Average.MAF will be the same as Gold.MAF.

If `--detail` is ON, SNP-wise r2 will be saved in `[OutputPrefix].RSquare`, which contains the following columns (tab-delimited).
* SNP.ID : variant name in the format of CHROM:POS:REF:ALT
* Allele.Frequency	: the allele frequency from the input allele frequency file, or calculated from the input validation file
* No.Samples	: the number of samples (samples with missing values in the validation file will be excluded)
* Imputation.R2	: squared Pearson correlation of genotype data for the variant between the validation file and imputation file
* Validation.AF	: the allele frequency calculated from genotype data in the input validation file
* Imputation.AF : the allele frequency calculated from genotype data in the input imputation file
