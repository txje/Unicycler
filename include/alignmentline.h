

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
    AlignmentLine(CommonKmer startPoint, int readLength, int refLength);

    void addPoint(CommonKmer newPoint);
    double getRelativeLineError();
    bool buildSeedChain(int minPointCount, float minAlignmentLength);
    std::string getDescriptiveString();

    double m_slope;
    double m_intercept;
    double m_error;
    double m_alignedReadLength;
    double m_alignedRefLength;

    int m_trimmedRefStart;
    int m_trimmedRefEnd;

    std::vector<CommonKmer> m_linePoints;
    String<TSeed> m_bridgedSeedChain;

private:
    int m_readLength;
    int m_refLength;

    // Values used for online calculation of linear regression.
    double m_meanX;
    double m_meanY;
    double m_ssX; // sum of squares for X
    double m_ssY; // sum of squares for Y
    double m_coMoment;
    double m_varX;
    double m_varY;
    double m_covariance;

    // void linearRegression(std::vector<CommonKmer> & pts, float * slope, float * intercept);
    void addBridgingSeed(TSeedSet & seedSet, int hStart, int vStart, int hEnd, int vEnd);
    void addSeedMerge(TSeedSet & seedSet, TSeed & seed);
};



// // LineFindingResults holds all of the information that results from the findAlignmentLines
// // function.
// class LineFindingResults {
// public:
//     LineFindingResults() {}
//     ~LineFindingResults() {for (size_t i = 0; i < m_lines.size(); ++i) delete m_lines[i];}

//     int m_milliseconds;
//     std::vector<AlignmentLine *> m_lines;
// };




// std::vector<AlignmentLine *> findAlignmentLines(CommonKmerSet * commonKmerSet,
//                                                 int readLength, int refLength, float expectedSlope,
//                                                 int verbosity, std::string & output,
//                                                 float lowScoreThreshold, float highScoreThreshold);

void getMeanAndStDev(std::vector<double> & v, double & mean, double & stdev);

std::string getKmerTable(std::vector<CommonKmer> & commonKmers);

std::string getSeedChainTable(String<TSeed> & seedChain);


#endif // ALIGNMENTLINE_H
