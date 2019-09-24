//
// Created by Ketian Yu on 9/23/19.
//

#include "RSquare.h"
#include <iostream>
#include <fstream>
#include <exception>

String RSquare::Analyze()
{
    if(!CheckVcfCompatibility())
        return "Input.VCF.Dose.Error";
    if(!CreateAggBins())
        return "Aggregate.Bins.Error";
    return "Success";
}

bool RSquare::CheckVcfCompatibility()
{
    if(!Validation.CheckValidity())
        return false;
    if(!Imputation.CheckValidity())
        return false;
    if(!CompatibleSamples())
        return false;
    return true;
}

bool RSquare::CompatibleSamples()
{
    if(Validation.ChrId != Imputation.ChrId)
    {
        cout << "\n ERROR !!! Input vcf files are on different chromosomes !!! " << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    if(Validation.NoSamples != Imputation.NoSamples)
    {
        cout << "\n ERROR !!! Input vcf files have different number of samples !!! " << endl;
        cout << " -- Validation: #Samples = " << Validation.NoSamples << endl;
        cout << " -- Imputation: #Samples = " << Imputation.NoSamples << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    NoSamples = Validation.NoSamples;

    for(int i=0; i<NoSamples; i++)
    {
        if(Validation.SampleNames[i]!=Imputation.SampleNames[i])
        {
            cout << "\n ERROR !!! Input vcf files have different sample orders !!! " << endl;
            cout << " -- Validation Sample " << i << " = " << Validation.SampleNames[i] << endl;
            cout << " -- Imputation Sample " << i << " = " << Imputation.SampleNames[i] << endl;
            cout << "\n Program Exiting ... \n\n";
            return false;
        }
    }

    Validation.clearSampleName();
    Imputation.clearSampleName();

    return true;
}

bool RSquare::CreateAggBins()
{
    if(!LoadBinsFile())
        return false;

    NoAggBins = BinsCutoffs.size()-1;
    aggBins.resize(NoAggBins);
    for(int i=0; i<NoAggBins; i++)
        aggBins[i].init(i, BinsCutoffs[i], BinsCutoffs[i+1]);
    return true;
}

bool RSquare::LoadBinsFile()
{
    if(myUserVariables->BinsFileName=="")
        return true;

    BinsCutoffs.clear();

    ifstream inFile(myUserVariables->BinsFileName);
    if(!inFile.is_open())
    {
        cout << "\n ERROR !!! Program could NOT open file : " << myUserVariables->BinsFileName << " !!! " << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    double value, lastValue=-1e-6;
    string line;

    while(getline(inFile, line))
    {
        if(line=="")
            continue;
        try
        {
            value = stof(line);
        }
        catch (exception& e)
        {
            cout << "\n ERROR !!! Bins file contains invalid value : " << line << " !!! " << endl;
            cout << " Please check file: " << myUserVariables->BinsFileName << endl;
            cout << "\n Program Exiting ... \n\n";
            inFile.close();
            return false;
        }
        if(value > 0.5)
        {
            cout << "\n ERROR !!! Bins file contains negative MAF cutoff > 0.5 : " << value << " !!! " << endl;
            cout << " Please check file: " << myUserVariables->BinsFileName << endl;
            cout << "\n Program Exiting ... \n\n";
            inFile.close();
            return false;
        }
        if(value < 0)
        {
            cout << "\n ERROR !!! Bins file contains MAF cutoff : " << value << " !!! " << endl;
            cout << " Please check file: " << myUserVariables->BinsFileName << endl;
            cout << "\n Program Exiting ... \n\n";
            inFile.close();
            return false;
        }
        if(value < lastValue)
        {
            cout << "\n ERROR !!! Decreasing cutoff in bin file : " << value << " < " << lastValue << " !!! " << endl;
            cout << " Please check file: " << myUserVariables->BinsFileName << endl;
            cout << "\n Program Exiting ... \n\n";
            inFile.close();
            return false;
        }
        lastValue = value;
        BinsCutoffs.push_back(value);
    }
    inFile.close();

    if(BinsCutoffs.size()<2)
    {
        cout << "\n ERROR !!! Insufficient #cutoffs in bin file : " << myUserVariables->BinsFileName << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }
    return true;
}