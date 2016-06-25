#ifndef GLOBAL_ALIGN_H
#define GLOBAL_ALIGN_H


#include <seqan/sequence.h>
#include "scoredalignment.h"


using namespace seqan;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {
    char * fullyGlobalAlignment(char * s1, char * s2,
                                int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore,
                                bool useBanding=false, int bandSize=1000);
}



ScoredAlignment * fullyGlobalAlignment(std::string s1, std::string s2,
                                       int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore,
                                       bool useBanding=false, int bandSize=1000);





#endif // GLOBAL_ALIGN_H
