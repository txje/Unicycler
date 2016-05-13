

#include "alignmentline.h"
#include "semiglobalalignment.h"



AlignmentLine::AlignmentLine(std::vector<CommonKmer> & commonKmers, int readLength, int refLength, float expectedSlope) :
    m_linePoints(commonKmers), m_readLength(readLength), m_refLength(refLength), m_expectedSlope(expectedSlope),
    m_slope(0.0), m_intercept(0.0), m_trimmedRefStart(0), m_trimmedRefEnd(0)
    // m_maxScore(0.0), m_areaUnderCurve(0.0)
{
    // int pointCount = m_linePoints.size();
    // for (int i = 0; i < pointCount; ++i)
    // {
    //     float score = m_linePoints[i].m_score;
    //     m_maxScore = std::max(score, m_maxScore);
    //     float rotatedVSize = 0.0;
    //     if (i > 0)
    //         rotatedVSize += (m_linePoints[i].m_rotatedVPosition - m_linePoints[i-1].m_rotatedVPosition) / 2.0;
    //     if (i < pointCount - 1)
    //         rotatedVSize += (m_linePoints[i+1].m_rotatedVPosition - m_linePoints[i].m_rotatedVPosition) / 2.0;
    //     m_areaUnderCurve += score * rotatedVSize;
    // }
}


// This function prepares the alignment line for use in an alignment. If it returns false, then the
// line is no good. If it returns true, that means it successfully built the m_bridgedSeedChain
// member.
bool AlignmentLine::buildSeedChain(int minPointCount, float minAlignmentLength) {

    // We are only interested in line groups which seem to span their full possible length (i.e.
    // not line groups caused by short, local alignments) and for which the band is reasonably 
    // long (to avoid things like 3 bp alignments).

    // Determine the mean position for the points and use it to calculate the band length.
    int pointCount = m_linePoints.size();
    double vSum = 0.0, hSum = 0.0;
    for (int i = 0; i < pointCount; ++i) {
        hSum += m_linePoints[i].m_hPosition;
        vSum += m_linePoints[i].m_vPosition;
    }
    double meanH = hSum / pointCount;
    double meanV = vSum / pointCount;
    double bandLength = getLineLength(meanH, meanV, m_expectedSlope, m_readLength, m_refLength);

    // Exclude alignments which are too short.
    if (bandLength < minAlignmentLength)
        return false;

    // Now we want to test whether the CommonKmers in the band seem to span the full band. To
    // do so, we get the std dev of the rotated H position and compare it to the expected std
    // dev of a uniform distribution.
    std::vector<double> rotatedHPositions;
    rotatedHPositions.reserve(pointCount);
    for (int i = 0; i < pointCount; ++i)
        rotatedHPositions.push_back(m_linePoints[i].m_rotatedHPosition);
    double meanRotatedH, rotatedHStdDev;
    getMeanAndStDev(rotatedHPositions, meanRotatedH, rotatedHStdDev);
    double uniformStdDev = bandLength / 3.464101615137754; // sqrt(12)

    // At least half of the uniform distribution's std dev is required.
    if (rotatedHStdDev < 0.5 * uniformStdDev)
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

    // Perform a simple least-squares linear regression. If the slope is too far from 1.0, we 
    // don't bother continuing to the seed chaining step as we'll throw the line out.
    linearRegression(m_linePoints, &m_slope, &m_intercept);
    if (m_slope < MIN_ALLOWED_SLOPE || m_slope > MAX_ALLOWED_SLOPE)
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
    return  "slope: " + std::to_string(m_slope) + ", intercept = " + std::to_string(m_intercept);
}




// Adapted from:
// http://stackoverflow.com/questions/11449617/how-to-fit-the-2d-scatter-data-with-a-line-with-c
void AlignmentLine::linearRegression(std::vector<CommonKmer> & pts, float * slope, float * intercept) {
    int n = pts.size();
    double sumH = 0.0, sumV = 0.0, sumHV = 0.0, sumHH = 0.0;
    for (int i = 0; i < n; ++i) {
        double hPos = pts[i].m_hPosition;
        double vPos = pts[i].m_vPosition;
        sumH += hPos;
        sumV += vPos;
        sumHV += hPos * vPos;
        sumHH += hPos * hPos;
    }
    double hMean = sumH / n;
    double vMean = sumV / n;
    double denominator = sumHH - (sumH * hMean);
    *slope = float((sumHV - sumH * vMean) / denominator);
    *intercept = float(vMean - (*slope * hMean));
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






// This function searches for lines in the 2D read-ref space that represent likely semi-global
// alignments.
std::vector<AlignmentLine *> findAlignmentLines(CommonKmerSet * commonKmerSet,
                                                int readLength, int refLength,
                                                int verbosity, std::string & output,
                                                float lowScoreThreshold, float highScoreThreshold,
                                                float mergeDistance) {
    std::vector<AlignmentLine *> returnedLines;

    int commonKmerCount = commonKmerSet->m_commonKmers.size();

    if (commonKmerCount < 2) {
        // if (verbosity > 3)
        //     output += "  No lines found, too few common k-mers (" + std::to_string(getTime() - startTime) + " ms)\n";
        return returnedLines;
    }

    if (verbosity > 4) {
        output += "  Common k-mer positions:\n";
        output += getKmerTable(commonKmerSet->m_commonKmers);
        output += "  Max score: " + std::to_string(commonKmerSet->m_maxScore) + "\n";
    }

    // Now group all of the line points. For a line group to form, a point must score above the
    // high threshold. The group will continue (in both directions) until the score drops below
    // the low threshold.
    std::vector<std::vector<CommonKmer> > lineGroups;
    bool lineInProgress = false;
    for (int i = 0; i < commonKmerCount; ++i) {
        if (lineInProgress) {
            if (commonKmerSet->m_commonKmers[i].m_score >= lowScoreThreshold)
                lineGroups.back().push_back(commonKmerSet->m_commonKmers[i]);
            else // This line group is done.
                lineInProgress = false;
        }
        else if (commonKmerSet->m_commonKmers[i].m_score >= highScoreThreshold) {
            // It's time to start a new line group!
            lineGroups.push_back(std::vector<CommonKmer>());
            lineInProgress = true;

            // Step backwards to find where the group should start (the first point over the low
            // threshold).
            int groupStartPoint = i;
            while (groupStartPoint >= 0 &&
                   commonKmerSet->m_commonKmers[groupStartPoint].m_score >= lowScoreThreshold)
                --groupStartPoint;
            ++groupStartPoint;

            // Add the initial group points.
            for (int j = groupStartPoint; j <= i; ++j)
                lineGroups.back().push_back(commonKmerSet->m_commonKmers[j]);
        }
    }
    if (verbosity > 4)
        output += "  Number of potential lines: " + std::to_string(lineGroups.size()) + "\n";

    // It's possible for one actual line group to be broken into multiple pieces because the scores
    // dipped low in the middle. So now we go back through our line groups and merge together those
    // that are sufficiently close to each other.
    std::vector<std::vector<CommonKmer> > mergedLineGroups;
    if (lineGroups.size() > 0)
        mergedLineGroups.push_back(lineGroups[0]);
    for (size_t i = 1; i < lineGroups.size(); ++i) {
        std::vector<CommonKmer> * previousGroup = &(mergedLineGroups.back());
        std::vector<CommonKmer> * thisGroup = &(lineGroups[i]);

        if (thisGroup->front().m_rotatedVPosition - previousGroup->back().m_rotatedVPosition <= mergeDistance)
            previousGroup->insert(previousGroup->end(), thisGroup->begin(), thisGroup->end());
        else
            mergedLineGroups.push_back(*thisGroup);
    }
    lineGroups = mergedLineGroups;
    if (verbosity > 4)
        output += "  Number of potential lines after merging: " + std::to_string(lineGroups.size()) + "\n";

    // Prepare the returned vector.
    for (size_t i = 0; i < lineGroups.size(); ++i)
        returnedLines.push_back(new AlignmentLine(lineGroups[i], readLength, refLength, commonKmerSet->m_expectedSlope));
    return returnedLines;
}



void getMeanAndStDev(std::vector<double> & v, double & mean, double & stdDev) {
    mean = 0.0;
    stdDev = 0.0;
    int count = v.size();
    if (count < 1)
        return;
    double devSum = 0.0;
    for (int i = 0; i < count; ++i)
        mean += v[i];
    mean /= count;
    for (int i = 0; i < count; ++i) {
        double dev = v[i] - mean;
        devSum += dev * dev;
    }
    stdDev = sqrt(devSum / v.size());
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


