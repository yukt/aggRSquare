//
// Created by Ketian Yu on 9/23/19.
//

#ifndef AGGRSQUARE_AGGBIN_H
#define AGGRSQUARE_AGGBIN_H

#include <algorithm>
#include <cmath>

class Record
{
public:
    unsigned int n;
    double sumX, sumY, sumX2, sumY2, sumXY;
    double sumX2rev, sumY2rev, sumXYrev;
    // X - Validation
    // Y = Imputation

    Record()
    {
        n     = 0;
        sumX  = 0.0;
        sumY  = 0.0;
        sumX2 = 0.0;
        sumY2 = 0.0;
        sumXY = 0.0;
        sumX2rev = 0.0;
        sumY2rev = 0.0;
        sumXYrev = 0.0;
    }

    void clear()
    {
        n     = 0;
        sumX  = 0.0;
        sumY  = 0.0;
        sumX2 = 0.0;
        sumY2 = 0.0;
        sumXY = 0.0;
        sumX2rev = 0.0;
        sumY2rev = 0.0;
        sumXYrev = 0.0;
    }

    void Update(int numGeno, double validationDose, double imputationDose)
    {
        if(validationDose<0 || imputationDose<0)
            return;
        n++;
        double x = 1.0/numGeno*validationDose;
        double y = 1.0/numGeno*imputationDose;

        sumX += x;
        sumY += y;
        sumX2 += x*x;
        sumY2 += y*y;
        sumXY += x*y;
        sumX2rev += (1-x)*(1-x);
        sumY2rev += (1-y)*(1-y);
        sumXYrev += (1-x)*(1-y);
    }

    double GetAlleleFreq() const
    {
        return sumX*1.0/n;
    }
};

class Bin
{
public:
    int BinId;
    int NoVariants;
    double lowerBound, upperBound;
    double MAF, R2, ValidationMAF, ImputationMAF;
    double sumX, sumY, sumX2, sumY2, sumXY;
    unsigned int n;

    Bin()
    {
        NoVariants = 0;
        MAF = 0;
        R2  = 0;
        ValidationMAF = 0;
        ImputationMAF = 0;
        sumX = 0.0;
        sumY = 0.0;
        sumX2 = 0.0;
        sumY2 = 0.0;
        sumXY = 0.0;
        n = 0;
    }
    void init(int id, double low, double up)
    {
        BinId = id;
        lowerBound = low;
        upperBound = up;
    }

    void AddRecord(Record& record, double freq)
    {
        double maf = std::min(freq, 1.0-freq);
        if(freq > maf)
        {
            sumX += record.n - record.sumX;
            sumY += record.n - record.sumY;
            sumX2 += record.sumX2rev;
            sumY2 += record.sumY2rev;
            sumXY += record.sumXYrev;
        }
        else
        {
            sumX += record.sumX;
            sumY += record.sumY;
            sumX2 += record.sumX2;
            sumY2 += record.sumY2;
            sumXY += record.sumXY;
        }
        MAF += maf;
        n += record.n;
        NoVariants++;
    }
    void Summarize()
    {
        if(NoVariants>0)
        {
            MAF = MAF*1.0/NoVariants;

            double EX   = sumX *1.0/n;
            double EY   = sumY *1.0/n;
            double varX = sumX2*1.0/n - EX*EX;
            double varY = sumY2*1.0/n - EY*EY;
            double cov  = sumXY*1.0/n - EX*EY;

            ValidationMAF = EX;
            ImputationMAF = EY;

            if(varX>0 && varY>0)
                R2 = 1.0*cov/varX*cov/varY;
        }
    }
};

#endif //AGGRSQUARE_AGGBIN_H
