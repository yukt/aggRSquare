//
// Created by Ketian Yu on 9/23/19.
//
#include "DosageFile.h"
#include "VcfFileReader.h"
#include "VcfHeader.h"

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
