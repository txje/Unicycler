

#ifndef ALIGNMENTLINE_H
#define ALIGNMENTLINE_H

#include <chrono>
#include <seqan/basic.h>
#include <seqan/seeds.h>
#include "kmers.h"
#include "settings.h"

using namespace seqan;

typedef Seed<Simple> TSeed;
typedef SeedSet<TSeed> TSeedSet;


// AlignmentLine takes a list of common k-mers and builds a seed chain that Seqan can use in a
// banded alignment.
class AlignmentLine {
public:
    AlignmentLine(std::vector<CommonKmer> & commonKmers, int readLength, int refLength, float expectedSlope);

    bool buildSeedChain(int minPointCount, float minAlignmentLength);
    // bool isBadLine() {return hasBadSlope() || hasNoSeedChain();}
    std::string getDescriptiveString();

    std::vector<CommonKmer> m_linePoints;
    int m_readLength;
    int m_refLength;
    float m_expectedSlope;
    // float m_maxScore;
    // float m_areaUnderCurve;
    float m_slope;
    float m_intercept;
    String<TSeed> m_bridgedSeedChain;
    int m_trimmedRefStart;
    int m_trimmedRefEnd;

private:
    void linearRegression(std::vector<CommonKmer> & pts, float * slope, float * intercept);
    void addBridgingSeed(TSeedSet & seedSet, int hStart, int vStart, int hEnd, int vEnd);
    void addSeedMerge(TSeedSet & seedSet, TSeed & seed);
    // bool hasBadSlope() {return m_slope < MIN_ALLOWED_SLOPE || m_slope > MAX_ALLOWED_SLOPE;}
    // bool hasNoSeedChain() {return length(m_bridgedSeedChain) == 0;}
};



// LineFindingResults holds all of the information that results from the findAlignmentLines
// function.
class LineFindingResults {
public:
    LineFindingResults() {}
    ~LineFindingResults() {for (size_t i = 0; i < m_lines.size(); ++i) delete m_lines[i];}

    int m_milliseconds;
    std::vector<AlignmentLine *> m_lines;
};




std::vector<AlignmentLine *> findAlignmentLines(CommonKmerSet * commonKmerSet,
                                                int readLength, int refLength, float expectedSlope,
                                                int verbosity, std::string & output,
                                                float lowScoreThreshold, float highScoreThreshold);

void getMeanAndStDev(std::vector<double> & v, double & mean, double & stdev);

std::string getKmerTable(std::vector<CommonKmer> & commonKmers);

std::string getSeedChainTable(String<TSeed> & seedChain);


#endif // ALIGNMENTLINE_H
