//
// Created by Ketian Yu on 9/23/19.
//

#include "RSquare.h"

String RSquare::Analyze()
{
    if(!CheckVcfCompatibility())
        return "Input.VCF.Dose.Error";


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