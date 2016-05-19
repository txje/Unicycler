

#include "alignmentline.h"
#include "semiglobalalignment.h"
#include <cmath>


AlignmentLine::AlignmentLine(CommonKmer startPoint, std::string readName, std::string refName,
                             int readLength, int refLength) :
    m_readName(readName), m_refName(refName),
    m_slope(0.0), m_intercept(0.0), m_error(0.0), m_alignedReadLength(0.0), m_alignedRefLength(0.0),
    m_trimmedRefStart(0), m_trimmedRefEnd(refLength),
    m_readLength(readLength), m_refLength(refLength),
    m_meanX(0.0), m_meanY(0.0), m_ssX(0.0), m_ssY(0.0), m_coMoment(0.0), m_varX(0.0), m_varY(0.0), m_covariance(0.0)
{
    addPoint(startPoint);
}


AlignmentLine::AlignmentLine(AlignmentLine * mergeLine1, AlignmentLine * mergeLine2) :
    m_slope(0.0), m_intercept(0.0), m_error(0.0), m_alignedReadLength(0.0), m_alignedRefLength(0.0),
    m_meanX(0.0), m_meanY(0.0), m_ssX(0.0), m_ssY(0.0), m_coMoment(0.0), m_varX(0.0), m_varY(0.0), m_covariance(0.0)
{
    m_readName = mergeLine1->m_readName;
    m_refName = mergeLine1->m_refName;
    m_readLength = mergeLine1->m_readLength;
    m_refLength = mergeLine1->m_refLength;
    m_trimmedRefStart = 0;
    m_trimmedRefEnd = m_refLength;

    for (size_t i = 0; i < mergeLine1->m_linePoints.size(); ++i)
        addPoint(mergeLine1->m_linePoints[i]);
    for (size_t i = 0; i < mergeLine2->m_linePoints.size(); ++i)
        addPoint(mergeLine2->m_linePoints[i]);
}


// This function adds a point to the alignment line and updates its linear regression estimates.
// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
// http://onlinestatbook.com/2/regression/accuracy.html
void AlignmentLine::addPoint(CommonKmer & newPoint) {
    m_linePoints.push_back(newPoint);
    double previousMeanY = m_meanY;

    double n = m_linePoints.size();
    double x = newPoint.m_hPosition;
    double y = newPoint.m_vPosition;

    double dX = x - m_meanX;
    double dY = y - m_meanY;

    m_meanX += dX / n;
    m_meanY += dY / n;

    m_ssX += dX * (x - m_meanX);
    m_ssY += dY * (y - m_meanY);
    m_coMoment += (x - m_meanX) * (y - previousMeanY);

    // Calculating the slope and intercept requires at least 2 points.
    if (n > 1.0) {
        m_varX = m_ssX / (n - 1.0);
        m_varY = m_ssY / (n - 1.0);
        m_covariance = m_coMoment / (n - 1.0);

        m_slope = m_covariance / m_varX;
        m_intercept = m_meanY - (m_slope * m_meanX);
    }

    // Calculating the error requires at least 3 points.
    if (n > 2.0) {
        double rSquared = (m_covariance * m_covariance) / (m_varX * m_varY);
        double var = ((1.0 - rSquared) * m_ssY) / (n - 2.0);
        m_error = sqrt(var);
    }

    // Calculate an estimate for the aligned lengths of read and reference.
    if (n > 1.0) {
        float xStart, yStart, xEnd, yEnd;
        getStartEndPoints(&xStart, &yStart, &xEnd, &yEnd);
        m_alignedReadLength = xEnd - xStart;
        m_alignedRefLength = yEnd - yStart;
    }
}


// Returns the line's error over the estimated reference length.
double AlignmentLine::getRelativeLineError() {
    if (m_alignedRefLength == 0.0)
        return 0.0;
    else
        return m_error / m_alignedRefLength;
}


// This function prepares the alignment line for use in an alignment. If it returns false, then the
// line is no good. If it returns true, that means it successfully built the m_bridgedSeedChain
// member.
bool AlignmentLine::buildSeedChain(int minPointCount, float minAlignmentLength) {

    // Exclude alignments which are too short.
    if (m_alignedReadLength < minAlignmentLength || m_alignedRefLength < minAlignmentLength)
        return false;

    // Exclude alignments with bad slopes.
    if (m_slope < MIN_ALLOWED_SLOPE || m_slope > MAX_ALLOWED_SLOPE)
        return false;

    // Throw out any point which is too divergent from its neighbours. To do
    // this, we determine the slope between each point and its nearest neighbours and how divergent
    // it is from the expected slope. We discard the points with the most divergent slopes.

    // Sort by rotated horizontal position so we proceed along the line from start to end. 
    std::sort(m_linePoints.begin(), m_linePoints.end(), [](const CommonKmer & a, const CommonKmer & b) {
        return a.m_rotatedHPosition < b.m_rotatedHPosition;   
    });

    // Store the max slope for each point in the line group.
    std::vector<double> maxSlopes;
    int pointCount = m_linePoints.size();
    maxSlopes.reserve(pointCount);
    for (int i = 0; i < pointCount; ++i) {
        CommonKmer * thisPoint = &(m_linePoints[i]);
        double maxSlope = 0.0;
        for (int j = -WORST_SLOPE_STEPS; j <= WORST_SLOPE_STEPS; ++j) {
            int neighbourIndex = i + j;
            if (neighbourIndex >= 0 && neighbourIndex < pointCount) {
                CommonKmer * neighbour = &(m_linePoints[neighbourIndex]);
                double slope = (thisPoint->m_rotatedVPosition - neighbour->m_rotatedVPosition) /
                               (thisPoint->m_rotatedHPosition - neighbour->m_rotatedHPosition);
                maxSlope = std::max(maxSlope, fabs(slope));
            }
        }
        maxSlopes.push_back(maxSlope);
    }

    // Determine a slope cutoff that will exclude the correct fraction of points.
    std::vector<double> sortedMaxSlopes = maxSlopes;
    std::sort(sortedMaxSlopes.begin(), sortedMaxSlopes.end());
    double slopeCutoff = sortedMaxSlopes[int(pointCount * (1.0 - WORST_SLOPE_FRACTION_TO_DISCARD))];

    // Create a new vector of line points, excluding those with excessively high slopes.
    std::vector<CommonKmer> fixedLinePoints;
    for (int i = 0; i < pointCount; ++i) {
        double maxSlope = maxSlopes[i];
        if (maxSlope <= slopeCutoff)
            fixedLinePoints.push_back(m_linePoints[i]);
    }
    m_linePoints = fixedLinePoints;

    // Remove any line groups with too few points.
    if (int(m_linePoints.size()) < minPointCount)
        return false;

    // When we do a Seqan banded alignment, we won't use the full reference sequence but just the
    // relevant part. Determine here what that relevant part is.
    int approxRefStart = int(round(m_intercept));
    int approxRefEnd = int(round(m_slope * m_readLength + m_intercept));
    m_trimmedRefStart = std::max(0, approxRefStart - PAD_SIZE);
    m_trimmedRefEnd = std::min(m_refLength, approxRefEnd + PAD_SIZE);
    int trimmedRefLength = m_trimmedRefEnd - m_trimmedRefStart;

    // Build a Seqan seed set using our common k-mers, offsetting by the trimmed reference start position.
    TSeedSet seedSet;
    for (size_t i = 0; i < m_linePoints.size(); ++i) {
        TSeed seed(m_linePoints[i].m_hPosition, m_linePoints[i].m_vPosition - m_trimmedRefStart, KMER_SIZE);
        addSeedMerge(seedSet, seed);
    }

    // We now get a Seqan global chain of the seeds.
    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());
    int seedsInChain = length(seedChain);
    if (seedsInChain == 0) {
        // if (verbosity > 4) 
        //     output += "Global chaining failed";
        return false;
    }

    // if (verbosity > 4) {
    //     output += "  Globally chained seeds before bridging\n";
    //     output += getSeedChainTable(seedChain);
    // }

    // Now we create a new seed chain with all of the gaps bridged. This will help keep alignment
    // in a narrow band, even when the seeds are spaced apart.
    TSeedSet bridgedSeedSet;

    // Create a seed bridge for the start of the chain by following the slope backwards from the
    // first seed.
    TSeed firstSeed = seedChain[0];
    int firstSeedReadStart = beginPositionH(firstSeed);
    int firstSeedRefStart = beginPositionV(firstSeed);
    if (firstSeedReadStart > 0 && firstSeedRefStart > 0) {
        double refPosAtStartOfRead = firstSeedRefStart - (m_slope * firstSeedReadStart);
        if (refPosAtStartOfRead >= 0.0)
            addBridgingSeed(bridgedSeedSet, 0, std::round(refPosAtStartOfRead), firstSeedReadStart, firstSeedRefStart);
        else
            addBridgingSeed(bridgedSeedSet, std::round(-refPosAtStartOfRead / m_slope), 0, firstSeedReadStart, firstSeedRefStart);
    }
    addSeedMerge(bridgedSeedSet, firstSeed);
    
    // Fill in any gaps in the middle of the seed chain.
    for (int i = 1; i < seedsInChain; ++i) {    
        TSeed seed1 = seedChain[i-1];
        TSeed seed2 = seedChain[i];
        int seed1ReadEnd = endPositionH(seed1);
        int seed1RefEnd = endPositionV(seed1);
        int seed2ReadStart = beginPositionH(seed2);
        int seed2RefStart = beginPositionV(seed2);
        addBridgingSeed(bridgedSeedSet, seed1ReadEnd, seed1RefEnd, seed2ReadStart, seed2RefStart);
        addSeedMerge(bridgedSeedSet, seed2);
    }

    // Create a seed bridge for the end of the chain by following the slope forwards from the
    // last seed.
    TSeed lastSeed = seedChain[seedsInChain - 1];
    int lastSeedReadEnd = endPositionH(lastSeed);
    int lastSeedRefEnd = endPositionV(lastSeed);
    if (lastSeedReadEnd < m_readLength && lastSeedRefEnd < trimmedRefLength) {
        double refPosAtStartOfRead = lastSeedRefEnd - (m_slope * lastSeedReadEnd);
        double refPosAtEndOfRead = (m_slope * m_readLength) + refPosAtStartOfRead;
        if (refPosAtEndOfRead <= trimmedRefLength)
            addBridgingSeed(bridgedSeedSet, lastSeedReadEnd, lastSeedRefEnd, m_readLength,
                            std::round(refPosAtEndOfRead));
        else
            addBridgingSeed(bridgedSeedSet, lastSeedReadEnd, lastSeedRefEnd,
                            std::round((trimmedRefLength - refPosAtStartOfRead) / m_slope), trimmedRefLength);
    }

    chainSeedsGlobally(m_bridgedSeedChain, bridgedSeedSet, SparseChaining());
    return true;
}

std::string AlignmentLine::getDescriptiveString() {
    return m_refName + ", slope = " + std::to_string(m_slope) + ", intercept = " + std::to_string(m_intercept);
}



// This function uses the slope, intercept and read/ref lengths to determine (x, y) coordinates for
// the start and end of the alignment line.
void AlignmentLine::getStartEndPoints(float * xStart, float * yStart, float * xEnd, float * yEnd) {
    if (m_intercept >= 0.0) {
        *xStart = 0.0;
        *yStart = m_intercept;
    }
    else {
        *xStart = -m_intercept / m_slope;
        *yStart = 0.0;
    }

    float yAtXEnd = (m_slope * m_readLength) + m_intercept;
    if (yAtXEnd <= m_refLength) {
        *xEnd = m_readLength;
        *yEnd = yAtXEnd;
    }
    else {
        *xEnd = (m_refLength - m_intercept) / m_slope;
        *yEnd = m_refLength;
    }
}


// This function takes the seed chain which should end at the point hStart, vStart. New seeds will be
// added to the chain to reach the point hEnd, vEnd.
void AlignmentLine::addBridgingSeed(TSeedSet & seedSet, int hStart, int vStart, int hEnd, int vEnd) {
    int hDiff = hEnd - hStart;
    int vDiff = vEnd - vStart;

    // If the start and end are against each other, then there's nothing to bridge.
    if (hDiff == 0 || vDiff == 0)
        return;

    TSeed seed(hStart, vStart, hEnd, vEnd);
    addSeedMerge(seedSet, seed);
}



// This function adds a seed to a seed set. First it tries to merge the seed, and if that doesn't
// work it adds it as a separate seed.
void AlignmentLine::addSeedMerge(TSeedSet & seedSet, TSeed & seed) {
    if (!addSeed(seedSet, seed, 1, Merge()))
        addSeed(seedSet, seed, Single());
}



std::string getSeedChainTable(String<TSeed> & seedChain) {
    std::string table = "\tH start\tH end\tV start\tV end\tLower diag\tUpper diag\n";
    int seedsInChain = length(seedChain);
    for (int i = 0; i < seedsInChain; ++i) {
        table += "\t" + std::to_string(beginPositionH(seedChain[i])) +
                 "\t" + std::to_string(endPositionH(seedChain[i])) +
                 "\t" + std::to_string(beginPositionV(seedChain[i])) +
                 "\t" + std::to_string(endPositionV(seedChain[i])) +
                 "\t" + std::to_string(lowerDiagonal(seedChain[i])) +
                 "\t" + std::to_string(upperDiagonal(seedChain[i])) + "\n";
    }
    return table;
}



std::string getKmerTable(std::vector<CommonKmer> & commonKmers) {
    std::string table = "\tSeq 1 pos\tSeq 2 pos\tRotated seq 1 pos\tRotated seq 2 pos\tScore\n";
    for (size_t i = 0; i < commonKmers.size(); ++i) {
        table += "\t" + std::to_string(commonKmers[i].m_hPosition) +
                 "\t" + std::to_string(commonKmers[i].m_vPosition) + 
                 "\t" + std::to_string(commonKmers[i].m_rotatedHPosition) +
                 "\t" + std::to_string(commonKmers[i].m_rotatedVPosition) +
                 "\t" + std::to_string(commonKmers[i].m_score) + "\n";
    }
    return table;
}



// Returns true if this alignment line is close to the other one.
bool AlignmentLine::isNear(AlignmentLine * other) {
    if (m_readName != other->m_readName || m_refName != other->m_refName)
        return false;
    float x11, y11, x12, y12;
    getStartEndPoints(&x11, &y11, &x12, &y12);
    float x21, y21, x22, y22;
    other->getStartEndPoints(&x21, &y21, &x22, &y22);
    return segmentsDistance(x11, y11, x12, y12, x21, y21, x22, y22) <= ALIGNMENT_LINE_MERGE_DISTANCE;
}

// Functions for determining distance between two line segments.
//   one segment is (x11, y11) to (x12, y12)
//   the other is   (x21, y21) to (x22, y22)
// http://stackoverflow.com/questions/2824478/shortest-distance-between-two-line-segments
float AlignmentLine::segmentsDistance(float x11, float y11, float x12, float y12,
                                      float x21, float y21, float x22, float y22) {
    if (segmentsIntersect(x11, y11, x12, y12, x21, y21, y22, y22))
        return 0.0;
    float distance = pointSegmentDistance(x11, y11, x21, y21, x22, y22);
    distance = std::min(distance, pointSegmentDistance(x12, y12, x21, y21, x22, y22));
    distance = std::min(distance, pointSegmentDistance(x21, y21, x11, y11, x12, y12));
    distance = std::min(distance, pointSegmentDistance(x22, y22, x11, y11, x12, y12));
    return distance;
}

bool AlignmentLine::segmentsIntersect(float x11, float y11, float x12, float y12,
                                      float x21, float y21, float x22, float y22) {
    float dx1 = x12 - x11;
    float dy1 = y12 - y11;
    float dx2 = x22 - x21;
    float dy2 = y22 - y21;
    float delta = (dx2 * dy1) - (dy2 * dx1);
    if (delta == 0.0)
        return false; // parallel segments
    float s = (dx1 * (y21 - y11) + dy1 * (x11 - x21)) / delta;
    float t = (dx2 * (y11 - y21) + dy2 * (x21 - x11)) / (-delta);
    return (0.0 <= s <= 1.0) && (0.0 <= t <= 1.0);
}

float AlignmentLine::pointSegmentDistance(float px, float py,
                                          float x1, float y1, float x2, float y2) {
    float dx = x2 - x1;
    float dy = y2 - y1;
    if (dx == 0.0 && dy == 0.0)
        return hypot(px - x1, py - y1);
    float t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy);
    if (t < 0.0) {
        dx = px - x1;
        dy = py - y1;
    }
    else if (t > 1.0) {
        dx = px - x2;
        dy = py - y2;
    }
    else {
        float nearX = x1 + t * dx;
        float nearY = y1 + t * dy;
        dx = px - nearX;
        dy = py - nearY;
    }
    return hypot(dx, dy);
}

float AlignmentLine::hypot(float dx, float dy) {
    return sqrt( (dx * dx) + (dy * dy) );
}


