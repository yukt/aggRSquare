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
    vector<Bin> aggBins;
    int NoSamples;

    RSquare(UserVariables &ThisUserVariables)
    :myUserVariables(&ThisUserVariables),
    Validation(ThisUserVariables.ValidationFileName, ThisUserVariables.formatValidation),
    Imputation(ThisUserVariables.ImputationFileName, ThisUserVariables.formatImputation)
    {
    }

    String Analyze();
    bool CheckVcfCompatibility();
    bool CompatibleSamples();
};

#endif //AGGRSQUARE_RSQUARE_H
