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

class DosageFile
{
public:
    String FileName;
    string Format;
    int NoMarkers, NoSamples;
    vector<string> SampleNames;
    string ChrId;

    VcfFileReader* InputDosageStream;
    VcfRecord*     CurrentRecord;
    VcfRecordGenotype* GenotypeInfo;
    int CurrentBp;
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
    double GetDosage(int SampleId);
    double GetNumGeno(int SampleId);

};

#endif //AGGRSQUARE_DOSAGEFILE_H
