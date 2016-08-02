
#ifndef ALIGNMENT_H
#define ALIGNMENT_H


#include <string>
#include <seqan/basic.h>
#include <seqan/align.h>

using namespace seqan;

enum CigarType {MATCH, INSERTION, DELETION, CLIP, NOTHING};


class ScoredAlignment {
public:
    ScoredAlignment(Align<Dna5String, ArrayGaps> & alignment, 
                    std::string & readName, std::string & refName,
                    int readLength, int refLength,
                    int refOffset, long long startTime, int bandSize,
                    bool startImmediately, bool goToEndSeq1, bool goToEndSeq2,
                    Score<int, Simple> & scoringScheme);
    std::string getFullString();
    std::string getShortDisplayString();
    bool isRevComp();
    int getReadAlignmentLength() {return m_readEndPos - m_readStartPos;}
    int getRefAlignmentLength() {return m_refEndPos - m_refStartPos;}

    std::string m_readName;
    std::string m_refName;
    int m_readLength;
    int m_refLength;
    int m_readStartPos;
    int m_readEndPos;
    int m_refStartPos;
    int m_refEndPos;
    std::string m_cigar;
    int m_rawScore;
    double m_scaledScore;
    int m_milliseconds;
    int m_bandSize;

private:
    CigarType getCigarType(char b1, char b2, bool alignmentStarted);
    std::string getCigarPart(CigarType type, int length);
    int getCigarScore(CigarType type, int length, Score<int, Simple> & scoringScheme,
                      std::string & readAlignment, std::string & refAlignment,
                      int alignmentPos);
};

long long getTime();

#endif // ALIGNMENT_H
