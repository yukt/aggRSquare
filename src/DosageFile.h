//
// Created by Ketian Yu on 9/23/19.
//

#ifndef AGGRSQUARE_DOSAGEFILE_H
#define AGGRSQUARE_DOSAGEFILE_H

#include <string>
#include <StringBasics.h>

using namespace std;

class DosageFile
{
public:
    String FileName;
    string Format;
    int NoMarkers, NoSamples;
    vector<string> SampleNames;
    string ChrId;

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

};

#endif //AGGRSQUARE_DOSAGEFILE_H
