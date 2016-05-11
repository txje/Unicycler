

#include "alignmentline.h"
#include "semiglobalalignment.h"

AlignmentLine::AlignmentLine(std::vector<CommonKmer> & commonKmers, int readLength, int refLength,
                             int verbosity, std::string & output) {

    // Perform a simple least-squares linear regression. If the slope is too far from 1.0, we 
    // don't bother continuing to the seed chaining step as we'll throw the line out.
    linearRegression(commonKmers, &m_slope, &m_intercept);
    if (m_slope < MIN_ALLOWED_SLOPE || m_slope > MAX_ALLOWED_SLOPE)
        return;

    // When we do a Seqan banded alignment, we won't use the full reference sequence but just the
    // relevant part. Determine here what that relevant part is.
    int approxRefStart = int(round(m_intercept));
    int approxRefEnd = int(round(m_slope * readLength + m_intercept));
    m_trimmedRefStart = std::max(0, approxRefStart - PAD_SIZE);
    m_trimmedRefEnd = std::min(refLength, approxRefEnd + PAD_SIZE);
    int trimmedRefLength = m_trimmedRefEnd - m_trimmedRefStart;


    // Build a Seqan seed set using our common k-mers, offsetting by the trimmed reference start position.
    TSeedSet seedSet;
    for (size_t i = 0; i < commonKmers.size(); ++i) {
        TSeed seed(commonKmers[i].m_hPosition, commonKmers[i].m_vPosition - m_trimmedRefStart, KMER_SIZE);
        addSeedMerge(seedSet, seed);
    }

    // We now get a Seqan global chain of the seeds.
    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());
    int seedsInChain = length(seedChain);
    if (seedsInChain == 0) {
        if (verbosity > 4) 
            output += "Global chaining failed";
        return;
    }

    if (verbosity > 4) {
        output += "  Globally chained seeds before bridging\n";
        output += getSeedChainTable(seedChain);
    }

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
    if (lastSeedReadEnd < readLength && lastSeedRefEnd < trimmedRefLength) {
        double refPosAtStartOfRead = lastSeedRefEnd - (m_slope * lastSeedReadEnd);
        double refPosAtEndOfRead = (m_slope * readLength) + refPosAtStartOfRead;
        if (refPosAtEndOfRead <= trimmedRefLength)
            addBridgingSeed(bridgedSeedSet, lastSeedReadEnd, lastSeedRefEnd, readLength,
                            std::round(refPosAtEndOfRead));
        else
            addBridgingSeed(bridgedSeedSet, lastSeedReadEnd, lastSeedRefEnd,
                            std::round((trimmedRefLength - refPosAtStartOfRead) / m_slope), trimmedRefLength);
    }

    chainSeedsGlobally(m_bridgedSeedChain, bridgedSeedSet, SparseChaining());
}

std::string AlignmentLine::getDescriptiveString() {
    return  "slope: " + std::to_string(m_slope) + ", intercept = " + std::to_string(m_intercept);
}




// Adapted from:
// http://stackoverflow.com/questions/11449617/how-to-fit-the-2d-scatter-data-with-a-line-with-c
void AlignmentLine::linearRegression(std::vector<CommonKmer> & pts, double * slope, double * intercept) {
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
    *slope = (sumHV - sumH * vMean) / denominator;
    *intercept = vMean - (*slope * hMean);
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
LineFindingResults * findAlignmentLines(std::string & readName, std::string & refName,
                                        int readLength, int refLength,
                                        double expectedSlope, int verbosity,
                                        KmerPositions * kmerPositions, std::string & output) {
    long long startTime = getTime();

    std::vector<CommonKmer> commonKmers = getCommonKmers(readName, refName, expectedSlope, kmerPositions);

    if (commonKmers.size() < 2) {
        if (verbosity > 3)
            output += "  No lines found, too few common k-mers (" + std::to_string(getTime() - startTime) + " ms):\n";
        return 0;
    }

    // We scale the scores relative to the expected k-mer density.
    double expectedDensity = 1.0 / pow(4.0, KMER_SIZE);
    if (verbosity > 4)
        output += "  Expected k-mer density: " + std::to_string(expectedDensity) + "\n";

    // Sort by rotated vertical position so lines should be roughly horizontal.
    std::sort(commonKmers.begin(), commonKmers.end(), [](const CommonKmer & a, const CommonKmer & b) {
        return a.m_rotatedVPosition < b.m_rotatedVPosition;   
    });

    // Score each point based on the number of other points in its band.
    int halfBandSize = BAND_SIZE / 2;
    double maxScore = 0.0;

    // There are four corners of the alignment rectangle which we also need to rotate.
    double rotationAngle = CommonKmer::getRotationAngle(expectedSlope);
    CommonKmer c1("", 0, 0, rotationAngle);
    CommonKmer c2("", 0, refLength, rotationAngle);
    CommonKmer c3("", readLength, refLength, rotationAngle);
    CommonKmer c4("", readLength, 0, rotationAngle);
    double c1Y = c1.m_rotatedVPosition;
    double c3Y = c3.m_rotatedVPosition;
    double c1BandLength = getLineLength(c1.m_hPosition, c1.m_vPosition,
                                        expectedSlope, readLength, refLength);
    double c3BandLength = getLineLength(c3.m_hPosition, c3.m_vPosition,
                                        expectedSlope, readLength, refLength);

    // Now we loop through the CommonKmer points, calculating their k-mer density (score) along the way.
    int commonKmerCount = commonKmers.size();
    for (int i = 0; i < commonKmerCount; ++i) {
        int bandStartIndex = std::max(i - halfBandSize, 0);
        int bandEndIndex = std::min(i + halfBandSize, commonKmerCount - 1);
        int thisBandSize = bandEndIndex - bandStartIndex;

        // Get the Y coordinates for the start and end of the band.
        double bandStartY, bandEndY;
        if (bandStartIndex < 0)
            bandStartY = c4.m_rotatedVPosition;
        else
            bandStartY = commonKmers[bandStartIndex].m_rotatedVPosition;
        if (bandEndIndex >= commonKmerCount)
            bandEndY = c2.m_rotatedVPosition;
        else
            bandEndY = commonKmers[bandEndIndex].m_rotatedVPosition;

        // Now we need to calculate the area of the band.
        double bandArea;

        // We'll need the starting point band length for all area possibilities, so we calculate
        // that now.
        double bandStartLength;
        if (bandStartIndex < 0)
            bandStartLength = 0.0;
        else
            bandStartLength = getLineLength(commonKmers[bandStartIndex].m_hPosition,
                                            commonKmers[bandStartIndex].m_vPosition,
                                            expectedSlope, readLength, refLength);

        // If both the start and end are in the middle of the rotated rectangle, then the area is a
        // parallelogram and calculating its area is easy.
        if (bandStartY >= c1Y && bandEndY <= c3Y)
            bandArea = (bandEndY - bandStartY) * bandStartLength;

        // Other cases are more complex, and we'll need the ending point band length too.
        else {
            double bandEndLength;
            if (bandEndIndex >= commonKmerCount)
                bandEndLength = 0.0;
            else
                bandEndLength = getLineLength(commonKmers[bandEndIndex].m_hPosition,
                                              commonKmers[bandEndIndex].m_vPosition,
                                              expectedSlope, readLength, refLength);

            // If both the start and end are in the bottom or top of the rotated rectangle, then the
            // area is a triangle/trapezoid.
            if ( (bandStartY <= c1Y && bandEndY <= c1Y) || (bandStartY >= c3Y && bandEndY >= c3Y) )
                bandArea = (bandEndY - bandStartY) * ((bandStartLength + bandEndLength) / 2.0);

            // If the start and end span a rectangle's corner, then the area is more complex. We need
            // to add both the parallelogram and trapezoid components.
            else if (bandStartY <= c1Y && bandEndY >= c1Y && bandEndY <= c3Y) {
                double trapezoidArea = (c1Y - bandStartY) * ((bandStartLength + c1BandLength) / 2.0);
                double parallelogramArea = (bandEndY - c1Y) * bandEndLength;
                bandArea = trapezoidArea + parallelogramArea;
            }
            else if (bandStartY >= c1Y && bandStartY <= c3Y && bandEndY >= c3Y) {
                double trapezoidArea = (bandEndY - c3Y) * ((c3BandLength + bandEndLength) / 2.0);
                double parallelogramArea = (c3Y - bandStartY) * bandStartLength;
                bandArea = trapezoidArea + parallelogramArea;
            }

            // The final, most complex scenario is when the band start and end span both C1 and C3.
            // This would be unusual, as it would require either a very large band or a very sparse
            // set of CommonKmers.
            else {
                double trapezoidArea1 = (c1Y - bandStartY) * ((bandStartLength + c1BandLength) / 2.0);
                double parallelogramArea = (c3Y - c1Y) * c1BandLength;
                double trapezoidArea2 = (bandEndY - c3Y) * ((c3BandLength + bandEndLength) / 2.0);
                bandArea = trapezoidArea1 + parallelogramArea + trapezoidArea2;
            }
        }

        // Now that we have the band area, we can get the density of CommonKmers in the band. Also,
        // we'll scale this to the expected level of CommonKmers (given a random sequence).
        double kmerDensity = thisBandSize / bandArea;
        double score = kmerDensity / expectedDensity;
        commonKmers[i].m_bandArea = bandArea;
        commonKmers[i].m_score = score;
        maxScore = std::max(maxScore, score);
    }

    if (verbosity > 4) {
        output += "  Common k-mer positions:\n";
        output += getKmerTable(commonKmers);
        output += "  Max score: " + std::to_string(maxScore) + "\n";
    }

    // Now group all of the line points. For a line group to form, a point must score above the
    // high threshold. The group will continue (in both directions) until the score drops below
    // the low threshold.
    std::vector<std::vector<CommonKmer> > lineGroups;
    bool lineInProgress = false;
    for (int i = 0; i < commonKmerCount; ++i) {
        if (lineInProgress) {
            if (commonKmers[i].m_score >= LOW_SCORE_THRESHOLD)
                lineGroups.back().push_back(commonKmers[i]);
            else // This line group is done.
                lineInProgress = false;
        }
        else if (commonKmers[i].m_score >= HIGH_SCORE_THRESHOLD) {
            // It's time to start a new line group!
            lineGroups.push_back(std::vector<CommonKmer>());
            lineInProgress = true;

            // Step backwards to find where the group should start (the first point over the low
            // threshold).
            int groupStartPoint = i;
            while (groupStartPoint >= 0 &&
                   commonKmers[groupStartPoint].m_score >= LOW_SCORE_THRESHOLD)
                --groupStartPoint;
            ++groupStartPoint;

            // Add the initial group points.
            for (int j = groupStartPoint; j <= i; ++j)
                lineGroups.back().push_back(commonKmers[j]);
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

        if (thisGroup->front().m_rotatedVPosition - previousGroup->back().m_rotatedVPosition <= MERGE_DISTANCE)
            previousGroup->insert(previousGroup->end(), thisGroup->begin(), thisGroup->end());
        else
            mergedLineGroups.push_back(*thisGroup);
    }
    lineGroups = mergedLineGroups;
    if (verbosity > 4)
        output += "  Number of potential lines after merging: " + std::to_string(lineGroups.size()) + "\n";

    // We are only interested in line groups which seem to span their full possible length (i.e.
    // not line groups caused by short, local alignments) and for which the band is reasonably 
    // long (to avoid things like 3 bp alignments).
    std::vector<std::vector<CommonKmer> > lengthFilteredLineGroups;
    for (size_t i = 0; i < lineGroups.size(); ++i) {
        // Determine the mean position for the line group and use it to calculate the band length.
        int groupCount = lineGroups[i].size();
        double vSum = 0.0, hSum = 0.0;
        for (int j = 0; j < groupCount; ++j) {
            hSum += lineGroups[i][j].m_hPosition;
            vSum += lineGroups[i][j].m_vPosition;
        }
        double meanH = hSum / groupCount;
        double meanV = vSum / groupCount;
        double bandLength = getLineLength(meanH, meanV, expectedSlope, readLength, refLength);

        // Exclude alignments which are too short.
        if (bandLength < MIN_ALIGNMENT_LENGTH) {
            if (verbosity > 4)
                output += "    Band too short: " + std::to_string(bandLength) + "\n";
            continue;
        }

        // Now we want to test whether the CommonKmers in the band seem to span the full band. To
        // do so, we get the std dev of the rotated H position and compare it to the expected std
        // dev of a uniform distribution.
        std::vector<double> rotatedHPositions;
        rotatedHPositions.reserve(groupCount);
        for (int j = 0; j < groupCount; ++j)
            rotatedHPositions.push_back(lineGroups[i][j].m_rotatedHPosition);
        double meanRotatedH, rotatedHStdDev;
        getMeanAndStDev(rotatedHPositions, meanRotatedH, rotatedHStdDev);
        double uniformStdDev = bandLength / 3.464101615137754; // sqrt(12)

        // At least half of the uniform distribution's std dev is required.
        if (rotatedHStdDev < 0.5 * uniformStdDev) {
            if (verbosity > 4)
                output += "    Distribution too narrow: " + std::to_string(rotatedHStdDev) + ", uniform std dev = "  + std::to_string(uniformStdDev) + "\n";
            continue;
        }

        lengthFilteredLineGroups.push_back(lineGroups[i]);
    }
    lineGroups = lengthFilteredLineGroups;
    if (verbosity > 4)
        output += "  Number of potential lines after length/span filtering: " + std::to_string(lineGroups.size()) + "\n";

    // For each line group, throw out any point which is too divergent from its neighbours. To do
    // this, we determine the slope between each point and its nearest neighbours and how divergent
    // it is from the expected slope. We discard the points with the most divergent slopes.

    // Parameters for this filtering step. TO DO: MAKE THESE ADJUSTABLE?
    double fractionToDiscard = 0.05;
    int steps = 3;

    for (size_t i = 0; i < lineGroups.size(); ++i) {
        std::vector<CommonKmer> * lineGroup = &(lineGroups[i]);

        // Sort by rotated horizontal position so we proceed along the line from start to end. 
        std::sort(lineGroup->begin(), lineGroup->end(), [](const CommonKmer & a, const CommonKmer & b) {
            return a.m_rotatedHPosition < b.m_rotatedHPosition;   
        });

        // Store the max slope for each point in the line group.
        std::vector<double> maxSlopes;
        maxSlopes.reserve(lineGroup->size());
        int groupSize = lineGroup->size();
        for (int j = 0; j < groupSize; ++j) {
            CommonKmer * thisPoint = &((*lineGroup)[j]);
            double maxSlope = 0.0;
            for (int k = -steps; k <= steps; ++k) {
                int neighbourIndex = j + k;
                if (neighbourIndex >= 0 && neighbourIndex < groupSize) {
                    CommonKmer * neighbour = &((*lineGroup)[neighbourIndex]);
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
        double slopeCutoff = sortedMaxSlopes[int(groupSize * (1.0 - fractionToDiscard))];
  
        // Create a new line group, excluding with excessively high slopes.
        std::vector<CommonKmer> fixedLineGroup;
        for (int j = 0; j < groupSize; ++j) {
            double maxSlope = maxSlopes[j];
            if (maxSlope <= slopeCutoff)
                fixedLineGroup.push_back((*lineGroup)[j]);
        }
        lineGroups[i] = fixedLineGroup;
    }

    // Remove any line groups with too few points.
    lineGroups.erase(std::remove_if(lineGroups.begin(), lineGroups.end(), 
                                    [](std::vector<CommonKmer> i) {return i.size() < MIN_POINT_COUNT;}),
                     lineGroups.end());
    if (verbosity > 4)
        output += "  Number of lines after point count filtering: " + std::to_string(lineGroups.size()) + "\n";

    if (lineGroups.size() == 0) {
        if (verbosity > 3) 
            output += "  No lines found (" + std::to_string(getTime() - startTime) + " ms):\n";
        return 0;
    }
    
    // Prepare the returned object.
    LineFindingResults * results = new LineFindingResults();
    results->m_milliseconds = getTime() - startTime;
    for (size_t i = 0; i < lineGroups.size(); ++i) {
        if (verbosity > 4) {
            output += "  Line " + std::to_string(i+1) + "\n";
            output += getKmerTable(lineGroups[i]);
        }

        AlignmentLine * line = new AlignmentLine(lineGroups[i], readLength, refLength, verbosity, output);
        if (line->isBadLine())
            delete line;
        else {
            results->m_lines.push_back(line);
        }
    }

    if (verbosity > 4)
        output += "  Number of lines after slope/chain filtering: " + std::to_string(lineGroups.size()) + "\n";
    if (verbosity > 3) {
        if (results->m_lines.size() == 0) {
            output += "  No lines found (" + std::to_string(results->m_milliseconds) + " ms):\n";
            delete results;
            return 0;
        }
        else {
            output += "  Lines found (" + std::to_string(results->m_milliseconds) + " ms):\n";
            for (size_t i = 0; i < results->m_lines.size(); ++i) {
                AlignmentLine * line = results->m_lines[i];
                output += "    " + line->getDescriptiveString() + "\n";
                if (verbosity > 4) {
                    output += "  Seed chain after bridging\n";
                    output += getSeedChainTable(line->m_bridgedSeedChain);
                }
            }
        }
    }

    return results;
}






// Given a point, a slope and rectangle bounds, this function returns the length of the line
// segment which goes through that point with that slope and is in those bounds.
double getLineLength(double x, double y, double slope, double xSize, double ySize) {
    double xStart, yStart, xEnd, yEnd;

    double yIntercept = y - (slope * x);
    if (yIntercept >= 0.0) {
        xStart = 0.0;
        yStart = yIntercept;
    }
    else {
        xStart = -yIntercept / slope;
        yStart = 0.0;
    }

    double yAtXEnd = (slope * xSize) + yIntercept;
    if (yAtXEnd <= ySize) {
        xEnd = xSize;
        yEnd = yAtXEnd;
    }
    else {
        xEnd = (ySize - yIntercept) / slope;
        yEnd = ySize;
    }

    double xLength = xEnd - xStart;
    double yLength = yEnd - yStart;
    return sqrt((xLength * xLength) + (yLength * yLength));
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
    std::string table = "\tSeq 1 pos\tSeq 2 pos\tRotated seq 1 pos\tRotated seq 2 pos\tBand area\tScore\n";
    for (size_t i = 0; i < commonKmers.size(); ++i) {
        table += "\t" + std::to_string(commonKmers[i].m_hPosition) +
                 "\t" + std::to_string(commonKmers[i].m_vPosition) + 
                 "\t" + std::to_string(commonKmers[i].m_rotatedHPosition) +
                 "\t" + std::to_string(commonKmers[i].m_rotatedVPosition) +
                 "\t" + std::to_string(commonKmers[i].m_bandArea) +
                 "\t" + std::to_string(commonKmers[i].m_score) + "\n";
    }
    return table;
}


