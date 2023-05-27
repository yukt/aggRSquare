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
    unsigned long NoCommonVariants, NoCommonVariantsAnalyzed;
    unsigned int NoAggBins, NoSamples;
    vector<double> BinsCutoffs;
    vector<Bin> aggBins;
    bool hasAlleleFreq;
    IFILE AlleleFreqFile;
    IFILE RSquareFile;
//    IFILE CommonSNPsFile;
    unsigned long NoAlleleFreqRead;
    Record CurrentRecord;
    unsigned long ValidationBufferBp;
    vector<Dosage> ValidationBuffer;

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
    void OutputRSquare(double freq) const;

    bool FindSamePosition();
    bool LoadValidationBuffer();
    bool ProcessCommonVariant();
    double ReadAlleleFreq();
    double GetAlleleFreq();
    void UpdateAggBins(double freq);
    bool UpdateInputRecords();

    // Sanity Checks
    bool CheckVcfCompatibility();
    bool CompatibleSamples();
    bool CreateAggBins();
    bool LoadBinsFile();
    bool CheckAlleleFreqFile();
    bool CheckSNPNameFormat(char* snp) const;
    bool OpenOutputFile();
};

#endif //AGGRSQUARE_RSQUARE_H
