//
// Created by Ketian Yu on 9/23/19.
//

#ifndef AGGRSQUARE_RSQUARE_H
#define AGGRSQUARE_RSQUARE_H

#include "aggBin.h"
#include "UserVariables.h"
#include "DosageFile.h"

using namespace std;

class RSquare
{
public:
    UserVariables* myUserVariables;
    DosageFile Validation, Imputation;
    int NoAggBins, NoCommonVariants, NoSamples;
    vector<double> BinsCutoffs;
    vector<Bin> aggBins;
    bool hasAlleleFreq;
    IFILE AlleleFreqFile;
    int NoAlleleFreqRead;
    Record CurrentRecord;

    RSquare(UserVariables &ThisUserVariables)
    :myUserVariables(&ThisUserVariables),
    Validation(ThisUserVariables.ValidationFileName, ThisUserVariables.formatValidation),
    Imputation(ThisUserVariables.ImputationFileName, ThisUserVariables.formatImputation),
    BinsCutoffs({0,0.0005,0.001,0.002,0.005,0.010,0.015,0.020,0.035,0.05,0.1,0.2,0.3,0.4,0.5}),
    hasAlleleFreq(false)
    {
        NoAlleleFreqRead = 0;
    }

    String Analyze();
    void OpenStreamInputFiles();
    bool EvaluateAggRSquare();
    void CloseStreamInputFiles();
    bool OutputAggRSquare();

    bool FindCommonVariant();
    bool ProcessCommonVariant();
    double GetAlleleFreq();
    void UpdateAggBins(double freq);
    bool UpdateInputRecords();

    // Sanity Checks
    bool CheckVcfCompatibility();
    bool CompatibleSamples();
    bool CreateAggBins();
    bool LoadBinsFile();
    bool CheckAlleleFreqFile();
    bool CheckSNPNameFormat(char* snp);
    bool OpenOutputFile();
};

#endif //AGGRSQUARE_RSQUARE_H
