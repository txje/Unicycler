

#ifndef ALIGNMENTLINE_H
#define ALIGNMENTLINE_H

#include <chrono>
#include <string>
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
    AlignmentLine(CommonKmer startPoint, std::string readName, std::string refName,
                  int readLength, int refLength);
    AlignmentLine(AlignmentLine * mergeLine1, AlignmentLine * mergeLine2);

    void addPoint(CommonKmer & newPoint);
    double getRelativeLineError();
    bool buildSeedChain(int minPointCount, float minAlignmentLength, int kSize);
    std::string getDescriptiveString();
    bool isNear(AlignmentLine * other);

    std::string m_readName;
    std::string m_refName;

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

    void addBridgingSeed(TSeedSet & seedSet, int hStart, int vStart, int hEnd, int vEnd);
    void addSeedMerge(TSeedSet & seedSet, TSeed & seed);

    void getStartEndPoints(float * xStart, float * yStart, float * xEnd, float * yEnd);
    float segmentsDistance(float x11, float y11, float x12, float y12, float x21, float y21, float x22, float y22);
    bool segmentsIntersect(float x11, float y11, float x12, float y12, float x21, float y21, float x22, float y22);
    float pointSegmentDistance(float px, float py, float x1, float y1, float x2, float y2);
    float hypot(float dx, float dy);
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


std::string getKmerTable(std::vector<CommonKmer> & commonKmers);

std::string getSeedChainTable(String<TSeed> & seedChain);


#endif // ALIGNMENTLINE_H
