
#ifndef SEMIGLOBALALIGNMENT_H
#define SEMIGLOBALALIGNMENT_H


#include <string>
#include <seqan/basic.h>
#include <seqan/align.h>

using namespace seqan;

enum CigarType {MATCH, INSERTION, DELETION, CLIP, NOTHING};


class SemiGlobalAlignment {
public:
    SemiGlobalAlignment(Align<Dna5String, ArrayGaps> & alignment, int refOffset, long long startTime,
              bool startImmediately, bool goToEnd, Score<int, Simple> & scoringScheme);
    std::string getFullString();
    std::string getShortDisplayString();

    int m_readStartPos;
    int m_readEndPos;
    int m_refStartPos;
    int m_refEndPos;
    std::string m_cigar;
    int m_rawScore;
    double m_scaledScore;
    int m_milliseconds;

private:
    CigarType getCigarType(char b1, char b2, bool alignmentStarted);
    std::string getCigarPart(CigarType type, int length);
    int getCigarScore(CigarType type, int length, Score<int, Simple> & scoringScheme);
};

long long getTime();

#endif // SEMIGLOBALALIGNMENT_H
