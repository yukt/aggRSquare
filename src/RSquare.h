//
// Created by Ketian Yu on 9/23/19.
//

#ifndef AGGRSQUARE_RSQUARE_H
#define AGGRSQUARE_RSQUARE_H

#include "aggBin.h"
#include "UserVariables.h"

using namespace std;

class DosageFile
{
    String FileName;
    string Format;
    int NoMarkers;
};


class RSquare
{
public:
    UserVariables myUserVariables;
    DosageFile Validation, Imputation;
    vector<Bin> aggBins;
};

#endif //AGGRSQUARE_RSQUARE_H
