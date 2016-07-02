#ifndef PATH_ALIGN_H
#define PATH_ALIGN_H


#include <seqan/sequence.h>
#include "scoredalignment.h"


using namespace seqan;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {
    char * pathAlignment(char * s1, char * s2,
                         int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore,
                         bool useBanding=false, int bandSize=1000);
}



ScoredAlignment * pathAlignment(std::string s1, std::string s2,
                                int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore,
                                bool useBanding=false, int bandSize=1000);





#endif // PATH_ALIGN_H
