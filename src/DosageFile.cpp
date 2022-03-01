//
// Created by Ketian Yu on 9/23/19.
//
#include "DosageFile.h"
#define MAXBP 999999999

inline int chr2int(string &chr)
{
    if(chr == "chrX" or chr == "X")
        return 23;
    if(chr == "chrY" or chr == "Y")
        return 24;
    if(chr == "chrMT" or chr == "MT")
        return 25;
    if(chr[0]=='c')
        return stoi(chr.substr(3, chr.size()-1));
    return stoi(chr);
}

double string2dosage(string &temp)
{
    double value = 0.0;
    char *pch = strtok((char *)temp.c_str(),"|/");
    while (pch != NULL)
    {
        if(strcmp(pch, ".")==0)
            return -1;
        try
        {
            value += stod(pch);
        }
        catch(exception& e)
        {
            return -1;
        }
        pch = strtok (NULL, "|/");
    }
    return value;
}

bool Dosage::LoadDosage(VcfRecordGenotype &GenotypeInfo)
{
    int NoSamples = GenotypeInfo.getNumSamples();
    dose.clear();
    dose.resize(NoSamples);
    for (int i=0; i<NoSamples; i++)
    {
        string temp = *GenotypeInfo.getString(Format, i);
        double value = string2dosage(temp);
        dose[i] = value;
    }
    return true;
}

bool DosageFile::CheckValidity()
{
    if(!ValidFileType())
        return false;
    if(!ValidSampleInfo())
        return false;
    return true;
}

bool DosageFile::ValidFileType()
{
    IFILE fileStream = ifopen(FileName, "r");
    if(fileStream)
    {
        string line;
        fileStream->readLine(line);
        ifclose(fileStream);
        if(line.length()>1)
        {
            string tempString;
            tempString=(line.substr(0,17));
            char temp[tempString.length() + 1];
            std::strcpy(temp,tempString.c_str());
            for (char *iter = temp; *iter != '\0'; ++iter)
            {
                *iter = std::tolower(*iter);
            }
            if(((string)temp).compare("##fileformat=vcfv")==0)
                return true;
        }
        cout << "\n ERROR !!! " << FileName << " is not a VCF file !!! " << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }
    cout << "\n ERROR !!! Program could NOT open file : " << FileName << endl;
    cout << "\n Program Exiting ... \n\n";
    return false;
}

bool ValidChrom(string chr)
{
    std::vector<string> ValidChromList {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16",
                                        "17","18","19","20","21","22","23","X","Y", "MT",
                                        "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                        "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                        "chr20","chr21","chr22","chr23","chrX","chrY", "chrMT"};
    for(string valid_chr : ValidChromList)
        if(chr == valid_chr)
            return true;
    return false;
}

bool DosageFile::ValidSampleInfo()
{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    SampleNames.clear();

    if (!inFile.open(FileName, header))
    {
        cout << "\n Program could NOT open file : " << FileName << endl;
        return false;
    }

    NoSamples = header.getNumSamples();
    if(NoSamples==0)
    {
        cout << "\n ERROR !!! None samples read from VCF File " << FileName << " !!! " << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    SampleNames.resize(NoSamples);
    for(int i=0; i<NoSamples; i++)
        SampleNames[i] = header.getSampleName(i);

    inFile.setSiteOnly(false);
    inFile.readRecord(record);
    inFile.close();

    ChrId = record.getChromStr();
    if(!ValidChrom(ChrId))
    {
        cout << "\n ERROR !!! " << FileName << " contains invalid chromosome " << ChrId << " !!! " << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    VcfRecordGenotype &ThisGenotype=record.getGenotypeInfo();
    const string* temp = ThisGenotype.getString(Format, 0);
    if(temp==NULL)
    {
        cout << "\n ERROR !!! " << FileName << " does not contain format " << Format << " !!! " << endl;
        cout << "\n Program Exiting ... \n\n";
        return false;
    }

    return true;
}

void DosageFile::OpenStream()
{
    InputDosageStream = new VcfFileReader();
    CurrentRecord     = new VcfRecord();

    VcfHeader header;
    InputDosageStream->open( FileName.c_str() , header);
    InputDosageStream->setSiteOnly(false);
    NoMarkers = 0;
    ReadRecord();
}

void DosageFile::CloseStream()
{
    delete InputDosageStream;
    delete CurrentRecord;
}

bool DosageFile::ReadRecord()
{
    if(InputDosageStream->readRecord(*CurrentRecord))
    {
        CurrentBp = CurrentRecord->get1BasedPosition();
        string chr(CurrentRecord->getChromStr());
        CurrentChr = chr2int(chr);
        CurrentVariantName = chr+":"+to_string(CurrentBp)+":"+ CurrentRecord->getRefStr()+":"+CurrentRecord->getAltStr();
        GenotypeInfo = &(CurrentRecord->getGenotypeInfo());
        NoMarkers++;
        return true;
    }
    CurrentChr = 99;
    CurrentBp = MAXBP;
    CurrentVariantName = "No:More:Variant:Period";
    return false;
}

double DosageFile::GetDosage(int SampleId)
{
    string temp=*GenotypeInfo->getString(Format,SampleId);
    double dosage = 0.0;
    char *pch = strtok((char *)temp.c_str(),"|/");
    while (pch != NULL)
    {
        if(strcmp(pch, ".")==0)
            return -1;
        try
        {
            dosage += stod(pch);
        }
        catch(exception& e)
        {
            cout << FileName << " : " << CurrentVariantName << " SampleId " << SampleId << "[" << *pch << "]" << endl;
            return -1;
        }
        pch = strtok (NULL, "|/");
    }
    return dosage;
}

double DosageFile::GetNumGeno(int SampleId){
    return CurrentRecord->getNumGTs(SampleId);
}