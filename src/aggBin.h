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
    int NoVariants;
    double lowerBound, upperBound;
    Bin()
    {
        NoVariants = 0;
    }
};

#endif //AGGRSQUARE_AGGBIN_H
