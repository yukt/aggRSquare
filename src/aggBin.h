//
// Created by Ketian Yu on 9/23/19.
//

#ifndef AGGRSQUARE_AGGBIN_H
#define AGGRSQUARE_AGGBIN_H

class Record
{
    int n;
    double sumX, sumY, sumX2, sumY2, sumXY;
public:
    Record()
    {
        n     = 0;
        sumX  = 0.0;
        sumY  = 0.0;
        sumX2 = 0.0;
        sumY2 = 0.0;
        sumXY = 0.0;
    }
};

class Bin : private Record
{
    int BinId;
    int NoVariants;
    double lowerBound, upperBound;
public:
    Bin()
    {
        NoVariants = 0;
    }
    void init(int id, double low, double up)
    {
        BinId = id;
        lowerBound = low;
        upperBound = up;
    }
};

#endif //AGGRSQUARE_AGGBIN_H
