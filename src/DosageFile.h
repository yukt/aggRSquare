//
// Created by Ketian Yu on 9/23/19.
//

#ifndef AGGRSQUARE_DOSAGEFILE_H
#define AGGRSQUARE_DOSAGEFILE_H

#include <string>
#include <StringBasics.h>
#include "VcfFileReader.h"
#include "VcfHeader.h"

using namespace std;

class Dosage
{
public:
    string VariantName;
    string Format;
    vector<double> dose;
    double AlleleFreq;

    Dosage(string &name, string &format)
    {
        VariantName = name;
        Format = format;
    }


    void LoadAlleleFreq(double freq)
    {
        AlleleFreq = freq;
    }
    bool LoadDosage(VcfRecordGenotype& GenotypeInfo);
};

class DosageFile
{
public:
    String FileName;
    string Format;
    unsigned long NoMarkers;
    unsigned int NoSamples;
    vector<string> SampleNames;
    string ChrId;

    VcfFileReader* InputDosageStream;
    VcfRecord*     CurrentRecord;
    VcfRecordGenotype* GenotypeInfo;
    int CurrentChr;
    unsigned long CurrentBp;
    string CurrentVariantName;

    DosageFile(String &filename, string &format)
    {
        FileName = filename;
        Format = format;
    }

    // Check Validity
    bool CheckValidity();
    bool ValidFileType();
    bool ValidSampleInfo();
    void clearSampleName() { SampleNames.clear();}

    // Read Records
    void OpenStream();
    void CloseStream();
    bool ReadRecord();
    double GetDosage(int SampleId) const;
    int GetNumGeno(int SampleId) const;

};

#endif //AGGRSQUARE_DOSAGEFILE_H
