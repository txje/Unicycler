
#ifndef RANDOM_ALIGNMENTS_H
#define RANDOM_ALIGNMENTS_H

#include <seqan/sequence.h>
#include <string>
#include <vector>
#include "semiglobalalignment.h"

using namespace seqan;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {

    char * getRandomSequenceAlignmentScores(int seqLength, int n,
                                            int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore);
    char * getRandomSequenceAlignmentErrorRates(int seqLength, int n,
                                               int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore);
}

SemiGlobalAlignment * fullyGlobalAlignment(std::string s1, std::string s2,
                                           int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore);

std::string getRandomSequence(int seqLength, std::mt19937 & gen, std::uniform_int_distribution<int> & dist);

char getRandomBase(std::mt19937 & gen, std::uniform_int_distribution<int> & dist);

void getMeanAndStDev(std::vector<double> & v, double & mean, double & stdev);


#endif // RANDOM_ALIGNMENTS_H
