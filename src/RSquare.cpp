//
// Created by Ketian Yu on 9/23/19.
//

#include "RSquare.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <exception>
#include <time.h>

String RSquare::Analyze()
{
    // Sanity Check
    if(!CheckVcfCompatibility())
        return "Input.VCF.Dose.Error";
    if(!CreateAggBins())
        return "Aggregate.Bins.Error";
    if(!CheckAlleleFreqFile())
        return "Allele.Freq.Error";
    if(!OpenOutputFile())
        return "File.Write.Error";

    // Analysis Starts
    OpenStreamInputFiles();

    if(!EvaluateAggRSquare())
    {
        CloseStreamInputFiles();
        return "Allele.Freq.Error";
    }

    CloseStreamInputFiles();

    if(!OutputAggRSquare())
        return "File.Write.Error";

    return "Success";
}

bool RSquare::EvaluateAggRSquare()
{
    NoCommonVariants = 0;
    NoCommonVariantsAnalyzed = 0;
    ValidationBuffer.clear();

    while(true)
    {
        if(!FindSamePosition())
            return true;
        if(!LoadValidationBuffer())
            return false;
        if(!ProcessCommonVariant())
            return true;
    }
}

bool RSquare::FindSamePosition()
{
    if(Imputation.CurrentVariantName=="No:More:Variant:Period")
        return false;
    if(Validation.CurrentVariantName=="No:More:Variant:Period")
        return false;

    while(true)
    {
        while(Validation.CurrentChr != Imputation.CurrentChr)
        {
            while(Validation.CurrentChr < Imputation.CurrentChr)
            {
                if(!Validation.ReadRecord())
                    return false;
            }
            while(Validation.CurrentChr > Imputation.CurrentChr)
            {
                if(!Imputation.ReadRecord())
                    return false;
            }
        }

        while(Validation.CurrentChr == Imputation.CurrentChr && Validation.CurrentBp != Imputation.CurrentBp)
        {
            if(Validation.CurrentBp < Imputation.CurrentBp)
            {
                if(!Validation.ReadRecord())
                    return false;
            }
            else
            {
                if(!Imputation.ReadRecord())
                    return false;
            }
        }

        if(Validation.CurrentChr == Imputation.CurrentChr && Validation.CurrentBp == Imputation.CurrentBp)
        {
            return true;
        }
    }
}

bool RSquare::LoadValidationBuffer()
{
    ValidationBuffer.clear();
    ValidationBufferBp = Validation.CurrentBp;

    while (Validation.CurrentBp == ValidationBufferBp)
    {
        Dosage CurrentDosage(Validation.CurrentVariantName, Validation.Format);
        if(hasAlleleFreq)
        {
            double freq = ReadAlleleFreq();
            if(freq < 0)
                return false;
            CurrentDosage.LoadAlleleFreq(freq);
        }
        if(CurrentDosage.LoadDosage(*(Validation.GenotypeInfo)))
            ValidationBuffer.push_back(CurrentDosage);
        if(!Validation.ReadRecord())
            break;
    }
    return true;
}

bool RSquare::ProcessCommonVariant()
{
    if (ValidationBuffer.empty())
        return true;
    while(Imputation.CurrentBp == ValidationBufferBp)
    {
        unsigned long NoBufferVariants = ValidationBuffer.size();
        for (int i=0; i<NoBufferVariants; i++)
        {
            if(Imputation.CurrentVariantName == ValidationBuffer[i].VariantName)
            {
                CurrentRecord.clear();
                vector<double>& ValidationDose = ValidationBuffer[i].dose;
                for(unsigned int j=0; j<NoSamples; j++)
                    CurrentRecord.Update(Imputation.GetNumGeno(j), ValidationDose[j], Imputation.GetDosage(j));
                double freq;
                if(hasAlleleFreq)
                    freq = ValidationBuffer[i].AlleleFreq;
                else
                    freq = CurrentRecord.GetAlleleFreq();
                UpdateAggBins(freq);
                break;
            }
        }
        if(!Imputation.ReadRecord())
            return true;
    }
    return true;
}

double RSquare::ReadAlleleFreq()
{
    string line;
    while(NoAlleleFreqRead < Validation.NoMarkers)
    {
        line = "";
        if(AlleleFreqFile->readLine(line)<0)
        {
            cout << "\n ERROR !!! Allele Freq file has fewer records than Validation VCF file !!! " << endl;
            cout << "\n Program Exiting ... \n\n";
            return -1;
        }
        if (line.length()==0 || line.find("#") == 0)
            continue;
        NoAlleleFreqRead++;
    }
    char* pch = strtok((char *)line.c_str(), "\t");
    string VariantName(pch);
    if(VariantName != Validation.CurrentVariantName)
    {
        cout << "\n ERROR !!! Allele Freq file does NOT match Validation VCF file !!! " << endl;
        cout << "\n Program Exiting ... \n\n";
        return -1;
    }
    pch = strtok(NULL, "\t");
    double AlleleFreq;
    try
    {
        AlleleFreq = stod(pch);
    }
    catch (exception& e)
    {
        cout << " Warning !!! " << "Skip " << VariantName << " ( invalid allele frequency " << pch << " ) !!!" << endl;
        return 99;
    }
    return AlleleFreq;
}

double RSquare::GetAlleleFreq()
{
    if(hasAlleleFreq)
        return ReadAlleleFreq();
    return CurrentRecord.GetAlleleFreq();
}

void RSquare::UpdateAggBins(double freq)
{
    double maf = min(freq, 1.0-freq);
    NoCommonVariants++;

    bool analyzed_flag = false;
    if(maf > BinsCutoffs[0])
    {
        for(int i=0; i<NoAggBins; i++)
        {
            if(maf <= BinsCutoffs[i+1])
            {
                aggBins[i].AddRecord(CurrentRecord, freq);
                if(myUserVariables->detail) OutputRSquare(freq);
                NoCommonVariantsAnalyzed++;
                analyzed_flag = true;
                break;
            }
        }
    }
//    if(!analyzed_flag)
//        cout << "    Skip " <<  Imputation.CurrentVariantName << " (MAF = " << freq << ")." << endl;
}

void RSquare::OutputRSquare(double freq) const
{
    unsigned int numGeno = CurrentRecord.n;
    double R2 = 0.0;
    double EX   = CurrentRecord.sumX *1.0/numGeno;
    double EY   = CurrentRecord.sumY *1.0/numGeno;
    double varX = CurrentRecord.sumX2*1.0/numGeno - EX*EX;
    double varY = CurrentRecord.sumY2*1.0/numGeno - EY*EY;
    double cov  = CurrentRecord.sumXY*1.0/numGeno - EX*EY;
    if(varX>0 && varY>0)
        R2 = 1.0*cov/varX*cov/varY;

    ifprintf(RSquareFile, "%s\t%f\t%d\t%f\t%f\t%f\n",Imputation.CurrentVariantName.c_str(), freq, numGeno, R2, EX, EY);
}



bool RSquare::UpdateInputRecords()
{
    if(!Validation.ReadRecord())
        return false;
    if(!Imputation.ReadRecord())
        return false;
    return true;
}


void RSquare::OpenStreamInputFiles()
{
    Validation.OpenStream();
    Imputation.OpenStream();
    if(hasAlleleFreq)
        AlleleFreqFile = ifopen(myUserVariables->AlleleFreqFileName, "r");
}

void RSquare::CloseStreamInputFiles()
{
    Validation.CloseStream();
    Imputation.CloseStream();
//    ifclose(CommonSNPsFile);
    if(hasAlleleFreq)
        ifclose(AlleleFreqFile);
    if(myUserVariables->detail)
        ifclose(RSquareFile);
}

bool RSquare::OutputAggRSquare()
{
    FILE* OutFile = fopen(myUserVariables->OutputPrefix+".aggRSquare", "w");

    time_t rawtime;
    struct tm * timeinfo;

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    fprintf(OutFile, "##filedate=%d.%d.%d\n",(timeinfo->tm_year + 1900),(timeinfo->tm_mon + 1) ,timeinfo->tm_mday);
    fprintf(OutFile, "##source=aggRSquare.v%s\n",VERSION);
    fprintf(OutFile, "##validation=%s\n", Validation.FileName.c_str());
    fprintf(OutFile, "##imputation=%s\n", Imputation.FileName.c_str());
    fprintf(OutFile, "##validationFormat=%s\n", Validation.Format.c_str());
    fprintf(OutFile, "##imputationFormat=%s\n", Imputation.Format.c_str());
    if(hasAlleleFreq) fprintf(OutFile, "##AF=%s\n", myUserVariables->AlleleFreqFileName.c_str());
    else fprintf(OutFile, "##Allele frequencies were calculated from validation file.\n");
    fprintf(OutFile,"##aggRSquare_Command=%s\n", myUserVariables->CommandLine.c_str());

    fprintf(OutFile, "#Bin.Aggregated.by.MAF\tAverage.MAF\t#Variants\tImputation.R2\tGold.MAF\tImputed.MAF\n");
    for(int i=0; i<NoAggBins; i++) {
        Bin &ThisBin = aggBins[i];
        ThisBin.Summarize();
        if (ThisBin.NoVariants > 0) {
            fprintf(OutFile, "(%f,%f]\t%f\t%lu\t%f\t%f\t%f\n",
                    ThisBin.lowerBound, ThisBin.upperBound,
                    ThisBin.MAF,
                    ThisBin.NoVariants,
                    ThisBin.R2,
                    ThisBin.ValidationMAF,
                    ThisBin.ImputationMAF);
        } else {
            printf(" Notice!!! Bin (%f,%f] is empty.\n", ThisBin.lowerBound, ThisBin.upperBound);
        }
    }
    fclose(OutFile);
    return true;
}


// Following are Sanity Checks
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

    for(unsigned int i=0; i<NoSamples; i++)
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
        cout << "\n ERROR !!! Program could NOT open Bins file : " << myUserVariables->BinsFileName << " !!! " << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    double value, lastValue=-1e-6;
    string line;

    while(getline(inFile, line))
    {
        if(line.length()==0 || line.find("#")==0)
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

bool RSquare::CheckAlleleFreqFile()
{
    if(myUserVariables->AlleleFreqFileName=="")
    {
        cout << " WARNING: Missing allele freq file ( --AF ) !!!" << endl;
        cout << "          Variants will be aggregated by allele freq calculated from validation vcf file." << endl;
        return true;
    }


    ifstream inFile(myUserVariables->AlleleFreqFileName);
    if(!inFile.is_open())
    {
        cout << "\n ERROR !!! Program could NOT open Allele Freq file : " << myUserVariables->AlleleFreqFileName << " !!! " << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    string line;
    while(getline(inFile, line)) {
        if (line.length()==0 || line.find("#") == 0)
            continue;
        break;
    }
    inFile.close();
    char* end_str;
    char* pch = strtok_r((char *)line.c_str(), "\t", &end_str);
    if(!CheckSNPNameFormat(pch))
    {
        cout << "\n ERROR !!! Invalid SNP in Allele Freq file : " << pch << " !!! " << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }
    if(end_str==NULL)
    {
        cout << "\n ERROR !!! Insufficient number of columns in Allele Freq file !!! " << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    hasAlleleFreq = true;

    return true;
}

bool RSquare::CheckSNPNameFormat(char* snp) const
{
    int NoColons = 0;
    for(int i=0; i<strlen(snp); i++)
    {
        if(snp[i]==':')
            NoColons++;
        if(NoColons > 3)
            return false;
    }
    if(NoColons < 3)
        return false;

    char* cno = strtok(snp, ":");
    if(cno != Validation.ChrId)
        return false;

    return true;
}

bool RSquare::OpenOutputFile()
{
    ofstream OutFile(myUserVariables->OutputPrefix+".aggRSquare");
    if(!OutFile.is_open())
    {
        cout << "\n ERROR !!! Program could NOT create the output file " << myUserVariables->OutputPrefix+".aggRSquare" << " !!!" << endl;
        cout <<" Please check your write permissions in the output directory\n OR maybe the output directory does NOT exist ...\n";
        cout << "\n Program Exiting ... \n\n";
        return false;
    }
    OutFile.close();

    if(myUserVariables->detail)
    {
        RSquareFile = ifopen(myUserVariables->OutputPrefix+".RSquare", "w");
        if(RSquareFile==NULL)
        {
            cout << "\n ERROR !!! Program could NOT create the output file " << myUserVariables->OutputPrefix+".RSquare" << " !!!" << endl;
            cout <<" Please check your write permissions in the output directory\n OR maybe the output directory does NOT exist ...\n";
            cout << "\n Program Exiting ... \n\n";
            return false;
        }

        ifprintf(RSquareFile, "#SNP.ID\tAllele.Frequency\tNo.Samples\tImputation.R2\tValidation.AF\tImputation.AF\n");
    }
    return true;

}