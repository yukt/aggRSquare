#include <getopt.h>
#include "RSquare.h"

void aggRSquareVersion();
void helpFile();

int main(int argc, char ** argv)
{
    RSquare myAnalysis;

    int c;
    static struct option loptions[] =
            {

                    {"validation",       required_argument, NULL, 'v'},
                    {"imputation",       required_argument, NULL, 'i'},
                    {"output",           required_argument, NULL, 'o'},
                    {"validationFormat", required_argument, NULL, 'f'},
                    {"imputationFormat", required_argument, NULL, 'g'},
                    {"AF",               required_argument, NULL, 'a'},
                    {"bins",             required_argument, NULL, 'b'},
                    {"help",             no_argument      , NULL, 'h'},
                    {NULL,0,NULL,0}
            };

    while ((c = getopt_long(argc, argv, "v:i:o:f:g:a:b:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'v': myAnalysis.myUserVariables.ValidationFileName = optarg;   break;
            case 'i': myAnalysis.myUserVariables.ImputationFileName = optarg;   break;
            case 'o': myAnalysis.myUserVariables.OutputPrefix       = optarg;   break;
            case 'f': myAnalysis.myUserVariables.formatValidation   = optarg;   break;
            case 'g': myAnalysis.myUserVariables.formatImputation   = optarg;   break;
            case 'a': myAnalysis.myUserVariables.AlleleFreqFileName = optarg;   break;
            case 'b': myAnalysis.myUserVariables.BinsFileName       = optarg;   break;
            case 'h': helpFile(); return 0;
            case '?': helpFile(); return 0;
            default: printf("[ERROR:] Unknown argument: %s\n", optarg);
        }
    }

    if(!myAnalysis.myUserVariables.CheckValidity()) return -1;

    aggRSquareVersion();
    myAnalysis.myUserVariables.Status();

    int start_time = time(0);
    myAnalysis.myUserVariables.CreateCommandLine(argc,argv);


    int time_tot = time(0) - start_time;
    return 0;
}

void aggRSquareVersion()
{
    printf("\n ------------------------------------------------------------------------------ \n");
    printf("             aggRSquare -- A Tool for Evaluating Imputation Accuracy    \n");
    printf(" ------------------------------------------------------------------------------\n");
    printf(" (c) 2019 - Ketian Yu, Sayantan Das, Goncalo Abecasis \n");
    cout<< " Version : " << VERSION<< ";\n Built   : " << DATE << " by " << USER << endl;
//    printf("\n URL = http://genome.sph.umich.edu/wiki/MetaMinimac2");
//    printf("\n GIT = https://github.com/yukt/MetaMinimac2.git\n");
    cout << endl;
}

void helpFile()
{
    aggRSquareVersion();
    printf( " About   : Calculate r2 between imputed dosages and true genotypes.\n");
    printf( " Usage   : aggRSquare [options] \n");
    printf( "\n");
    printf( " Options :\n");
    printf( "   -v [Validation.vcf.gz]         // [Required] Input Validation File\n");
    printf( "   -i [Imputation.vcf.gz]         // [Required] Input Imputation File\n");
    printf( "   -o [OutputPrefix]              // [Required] Output Prefix\n");
    printf( "   --validationFormat [GT/DS/GP]  // [Optional] Genotype info format (Default: GT)\n");
    printf( "   --imputationFormat [GT/DS/GP]  // [Optional] Genotype info format (Default: DS)\n");
    printf( "   --AF [AlleleFrequency File]    // [Optional] See Note (a)\n");
    printf( "   --bins [Bins File]             // [Optional] Default: See defaultBins.txt\n");
    printf( "   -h, --help                          If ON, detailed help on options and usage. \n");
    printf( " Note: (a) AlleleFrequency file must contain two tab-delimited columns with header 'SNP\tAF',\n");
    printf( "           and SNPs should be same as the validation vcf file. \n");
    cout<<endl<<endl;
}
