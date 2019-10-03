//
// Created by Ketian Yu on 9/23/19.
//

#ifndef AGGRSQUARE_USERVARIABLES_H
#define AGGRSQUARE_USERVARIABLES_H

#include <string>
#include <StringBasics.h>

using namespace std;

class UserVariables {
public:
    String ValidationFileName, ImputationFileName;
    String AlleleFreqFileName;
    String BinsFileName;
    String OutputPrefix;

    string formatValidation, formatImputation;
    string CommandLine;

    bool detail;

    UserVariables() {
        ValidationFileName = "";
        ImputationFileName = "";
        AlleleFreqFileName = "";
        BinsFileName       = "";
        OutputPrefix       = "";
        formatValidation   = "GT";
        formatImputation   = "DS";
        detail = false;
    };

    void CreateCommandLine(int argc, char ** argv)
    {
        int len = 0;

        for (int i=0; i<argc; i++)
            len += strlen(argv[i]) + 1;



        char MyCommandLine[len];
        strcpy(MyCommandLine,argv[0]);

        for (int i=1; i<argc; i++)
        {
            strcat(MyCommandLine, " ");
            strcat(MyCommandLine, argv[i]);
        }
        CommandLine=MyCommandLine;
    }

    void Status()
    {
        cout << " Command Line Options: " << endl;
        printf( "   --validation [%s]\n", ValidationFileName.c_str());
        printf( "   --imputation [%s]\n", ImputationFileName.c_str());
        printf( "   --output     [%s]\n", OutputPrefix.c_str());
        printf( "   --AF         [%s]\n", AlleleFreqFileName.c_str());
        printf( "   --bins       [%s]\n", BinsFileName.c_str());
        printf( "   --validationFormat [%s]\n", formatValidation.c_str());
        printf( "   --imputationFormat [%s]\n", formatImputation.c_str());
        printf( "   --detail     [%s]\n", detail?"ON":"OFF");
        printf("\n\n");
    }

    bool CheckValidity()
    {
        if(ValidationFileName == ""){
            cout<< " Missing -v [--validation], a required parameter.\n\n";
            cout<< " Try -h [--help] for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }
        if(ImputationFileName == ""){
            cout<< " Missing -i [--imputation], a required parameter.\n\n";
            cout<< " Try -h [--help] for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }
        if(OutputPrefix == ""){
            cout<< " Missing -o [--output], a required parameter.\n\n";
            cout<< " Try -h [--help] for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }
        if(formatValidation!="DS" and formatValidation!="GT" and formatValidation!="HDS")
        {
            printf("ERROR !!! Invalid argument --validationFormat %s !!!\n", formatValidation.c_str());
            return false;
        }

        if(formatImputation!="DS" and formatImputation!="GT" and formatImputation!="HDS")
        {
            printf("ERROR !!! Invalid argument --validationFormat %s !!!\n", formatValidation.c_str());
            return false;
        }
        return true;
    }
};
#endif //AGGRSQUARE_USERVARIABLES_H
