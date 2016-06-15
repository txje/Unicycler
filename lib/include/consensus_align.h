
#ifndef CONSENSUS_ALIGN_H
#define CONSENSUS_ALIGN_H

#include <seqan/basic.h>
#include <seqan/score.h>
#include <seqan/consensus.h>

using namespace seqan;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {

    char * multipleSequenceAlignment(char * fullSpanSequences[], char * fullSpanQualities[], int fullSpanCount, 
                                     char * startOnlySequences[], char * startOnlyQualities[], int startOnlyCount, 
                                     char * endOnlySequences[], char * endOnlyQualities[], int endOnlyCount, 
                                     int bandwidth, int matchScore, int mismatchScore,
                                     int gapOpenScore, int gapExtensionScore);

}


char getMostCommonBase(std::vector<char> & bases, std::vector<char> & qualities);

int scoreAlignment(std::string & seq1, std::string & seq2, int startPos, int endPos,
                   int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore);

void fillOutQualities(std::vector<std::string> & sequences, std::vector<std::string> & qualities);
void padToLength(std::vector<std::string> & sequences, std::vector<std::string> & qualities, int length, bool putAtStart);
void cArrayToCppVector(char * seqArray[], char * qualArray[], int count,
                       std::vector<std::string> & seqVector, std::vector<std::string> & qualVector);


#endif // CONSENSUS_ALIGN_H
