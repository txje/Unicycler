
#ifndef COMMONKMERSET_H
#define COMMONKMERSET_H

#include <string>
#include "kmers.h"
#include "alignmentline.h"

// This class holds a set of common k-mers for a particular read-ref combination.
class CommonKmerSet {
public:
    CommonKmerSet(std::string & readName, std::string & refName, int readLength, int refLength,
                  float expectedSlope, KmerPositions * kmerPositions);

    std::string m_readName;
    std::string m_refName;

    int m_readLength;
    int m_refLength;
    float m_expectedSlope;

    std::vector<CommonKmer> m_commonKmers;
    
    float m_maxScore;
    int m_maxScoreIndex;

    AlignmentLine * extractAlignmentLine();

private:
    void scorePoints();
    void findMaxScore();
};

float getY(std::vector<CommonKmer> & commonKmers, int i, int count, CommonKmer & c2, CommonKmer & c4);

float getLineLength(float x, float y, float slope, float xSize, float ySize);


#endif // COMMONKMERSET_H

