#include "semi_global_align.h"

#include <seqan/align.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <math.h>

#include "settings.h"
#include <seqan/basic.h>
#include <seqan/seeds.h>


char * semiGlobalAlignment(char * readNameC, char * readSeqC, int verbosity,
                           char * minimapAlignmentsStr, SeqMap * refSeqs,
                           int matchScore, int mismatchScore, int gapOpenScore,
                           int gapExtensionScore, double lowScoreThreshold, bool returnBad,
                           int sensitivityLevel) {
    int kSize = LEVEL_0_KMER_SIZE;
    if (sensitivityLevel == 1)
        kSize = LEVEL_1_KMER_SIZE;
    else if (sensitivityLevel == 2)
        kSize = LEVEL_2_KMER_SIZE;
    else if (sensitivityLevel == 3)
        kSize = LEVEL_3_KMER_SIZE;

    std::string output;
    std::string returnString;
    std::vector<ScoredAlignment *> returnedAlignments;

//    std::cout << std::endl;  // TEMP

    // Change the read name and sequence to C++ strings.
    std::string readName(readNameC);
    std::string posReadName = readName + "+";
    std::string negReadName = readName + "-";
    std::string posReadSeq(readSeqC);
    std::string negReadSeq;  // Will make later, if necessary.
    int readLength = posReadSeq.length();


    std::vector<std::string> minimapAlignments = splitString(minimapAlignmentsStr, ';');
    if (verbosity > 2) {
        output += "minimap alignments:\n";
        for (auto minimapAlignment : minimapAlignments)
            output += "    " + minimapAlignment + "\n";
    }


    // Debugging information for use in R.
    if (verbosity > 3) {  // only at very high verbosities
        output += "R:library(ggplot2)\n";
        output += "R:dot.plot.1 <- function(all_points) {ggplot() + geom_point(data=all_points,  aes(x=X1, y=X2), size=0.1, alpha=0.1, shape=19) + theme_bw() + coord_equal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))}\n";
        output += "R:dot.plot.2 <- function(all_points, trace_dots) {ggplot() + geom_point(data=all_points,  aes(x=X1, y=X2), size=0.1, alpha=0.02, shape=19) + geom_point(data=trace_dots,  aes(x=X1, y=X2), size=0.1, alpha=1, shape=19, colour=\"red\") + theme_bw() + coord_equal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))}\n";
        output += "R:dot.plot.3 <- function(all_points, filtered_data, trace_dots) {ggplot() + geom_point(data=all_points,  aes(x=X1, y=X2), size=0.1, alpha=0.02, shape=19) + geom_point(data=filtered_data,  aes(x=X1, y=X2), size=0.1, alpha=1, shape=19, colour=\"green\") + geom_point(data=trace_dots,  aes(x=X1, y=X2), size=0.1, alpha=1, shape=19, colour=\"red\") + theme_bw() + coord_equal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))}\n";
    }

    // For each minimap alignment we find the appropriate part of the reference sequence.
    RefRangeMap refRanges;
    for (size_t i = 0; i < minimapAlignments.size(); ++i) {
        std::string minimapStr = minimapAlignments[i];
        std::vector<std::string> minimapStrParts = splitString(minimapStr, ',');

        int readStart = std::stoi(minimapStrParts[0]);
        int readEnd = std::stoi(minimapStrParts[1]);
        char readStrand = minimapStrParts[2][0];
        bool posStrand = readStrand == '+';

        std::string refName = minimapStrParts[3];
        int refStart = std::stoi(minimapStrParts[4]);
        int refEnd = std::stoi(minimapStrParts[5]);
        std::string & refSeq = refSeqs->at(refName);

//        std::cout << readStart << "  " << readEnd << "  " << readStrand << "  ";  // TEMP
//        std::cout << refName << "  " << refStart << "  " << refEnd << std::endl;  // TEMP

        StartEndRange refRange = getRefRange(refStart, refEnd, refSeq.length(),
                                             readStart, readEnd, readLength, posStrand);
//        std::cout << refRange.first << "  " << refRange.second << std::endl;  // TEMP

        // The first time we see the reference/strand, initialise it with an empty vector.
        std::string refNameAndStrand = refName + readStrand;
        if (refRanges.find(refNameAndStrand) == refRanges.end())
            refRanges[refNameAndStrand] = std::vector<StartEndRange>();

        refRanges[refNameAndStrand].push_back(refRange);
    }

    // Simplify the reference ranges by combining overlapping ranges.
    RefRangeMap simplifiedRefRanges;
    for(auto const& r : refRanges) {
        std::string refNameAndStrand = r.first;
//        std::cout << refNameAndStrand << std::endl;  // TEMP
        std::vector<StartEndRange> ranges = r.second;
//        for (int i = 0; i < ranges.size(); ++i)  // TEMP
//            std::cout << "(" << ranges[i].first << ","<< ranges[i].second << ") ";  // TEMP
//        std::cout << std::endl;  // TEMP
        std::vector<StartEndRange> simplifiedRanges = simplifyRanges(ranges);
        simplifiedRefRanges[refNameAndStrand] = simplifiedRanges;
//        for (int i = 0; i < simplifiedRanges.size(); ++i)  // TEMP
//            std::cout << "(" << simplifiedRanges[i].first << ","<< simplifiedRanges[i].second << ") ";  // TEMP
//        std::cout << std::endl;  // TEMP
    }

    // Make a new KmerPositions object for the read. We'll actually add positions later as
    // necessary (because we may not need both the positive strand or the negative strand).
    KmerPositions readKmerPositions;
    bool posPositions = false, negPositions = false;

    // Align to each reference range.
    for(auto r : simplifiedRefRanges) {
        std::string refName = r.first;
        char readStrand = refName.back();
        bool posStrand = readStrand == '+';
        refName.pop_back();
        std::string & refSeq = refSeqs->at(refName);
        std::vector<StartEndRange> ranges = r.second;

        // Prepare some stuff for the read.
        std::string * readSeq;
        KmerPosMap * kmerPositions;
        if (posStrand) {
            if (!posPositions) {
                readKmerPositions.addPositions(posReadName, posReadSeq, kSize);
                posPositions = true;
            }
            readSeq = &posReadSeq;
            kmerPositions = readKmerPositions.getKmerPositions(posReadName);
        }
        else {  // negative strand
            if (!negPositions) {
                negReadSeq = getReverseComplement(posReadSeq);
                readKmerPositions.addPositions(negReadName, negReadSeq, kSize);
                negPositions = true;
            }
            readSeq = &negReadSeq;
            kmerPositions = readKmerPositions.getKmerPositions(negReadName);
        }

        // Work on each range (there's probably just one, but there could be more).
        for (auto range : ranges) {
            std::vector<ScoredAlignment *> a =
                alignReadToReferenceRange(refSeqs, refName, range, refSeq.length(), readName,
                                          readStrand, kmerPositions, kSize, readSeq, matchScore,
                                          mismatchScore, gapOpenScore, gapExtensionScore,
                                          sensitivityLevel, verbosity, output);
            returnedAlignments.insert(returnedAlignments.end(), a.begin(), a.end());
        }
    }

    // The returned string is semicolon-delimited. The last part is the console output and the
    // other parts are alignment description strings.
    for (auto alignment : returnedAlignments) {
        if (alignment != 0)
            returnString += alignment->getFullString() + ";";
    }
    returnString += output;

    return cppStringToCString(returnString);
}


std::vector<ScoredAlignment *> alignReadToReferenceRange(SeqMap * refSeqs, std::string refName,
                                                         StartEndRange refRange, int refLen,
                                                         std::string readName, char readStrand,
                                                         KmerPosMap * kmerPositions, int kSize,
                                                         std::string * readSeq, int matchScore,
                                                         int mismatchScore, int gapOpenScore,
                                                         int gapExtensionScore,
                                                         int sensitivityLevel,
                                                         int verbosity, std::string & output) {
    long long startTime = getTime();
    std::vector<ScoredAlignment *> alignments;

    int refStart = refRange.first;
    int refEnd = refRange.second;
    int readLen = readSeq->length();

//    std::cout << "(" << refName << "," << readStrand << "," << refStart << "," << refEnd << ") ";  // TEMP
//    std::cout << std::flush; // TEMP
    std::string trimmedRefSeq = refSeqs->at(refName).substr(refStart, refEnd-refStart);
    int trimmedRefLen = trimmedRefSeq.length();

    // Find all common k-mer positions.
    std::vector<CommonKmer> commonKmers;
    int maxI = trimmedRefLen - kSize + 1;
    for (int i = 0; i < maxI; ++i) {
        std::string refKmer = trimmedRefSeq.substr(i, kSize);
        if (kmerPositions->find(refKmer) != kmerPositions->end() ) {  // if k-mer is in the read
            std::vector<int> & readPositions = kmerPositions->at(refKmer);
            for (size_t j = 0; j < readPositions.size(); ++j)
                commonKmers.emplace_back(readPositions[j], i);
        }
    }

    // Debugging information for use in R.
    if (verbosity > 3) {  // only at very high verbosities
        std::ofstream allPointsFile;
        std::string filename = readName + readStrand + "_" + refName + "_all_points.tsv";
        allPointsFile.open(filename);
        for (auto k : commonKmers)
            allPointsFile << k.m_hPosition << "\t" << k.m_vPosition << "\n";
        allPointsFile.close();
        output += "R:all.points <- read_delim(\"" + filename + "\", \"\t\", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)\n";
    }


//    for (int i = 0; i < commonKmers.size(); ++i)   // TEMP
//        std::cout << commonKmers[i].m_hPosition << "\t" << commonKmers[i].m_vPosition << "\n";  // TEMP
//    std::cout << std::endl;  // TEMP

    PointCloud cloud;
    addKmerPointsToNanoflann(cloud, commonKmers);
    my_kd_tree_t index(2, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
    index.buildIndex();

    double highestDensityScore = 0.0;
    Point highestDensityPoint = getHighestDensityPoint(100, cloud, index, trimmedRefSeq, readSeq,
                                                       &highestDensityScore);

    Point p = highestDensityPoint;
    std::vector<Point> traceDots;
    traceDots.push_back(p);
//    std::cout << std::endl;
//    std::cout << "(" << p.x << "," << p.y << ") " << std::flush;

    int smallLineTracingStepSize = 250;
    int smallSearchRadius = 500;

    // If the line is 'lost' then we will switch to larger steps to increase our change of
    //'finding' it again.
    int largeLineTracingStepSize = 500;
    int largeSearchRadius = 1000;

    // Start the point collection using points around the starting point.
    std::unordered_set<Point> pointSet;
    std::vector<Point> nearbyPoints = radiusSearchAroundPoint(p, smallSearchRadius, cloud, index);
    for (auto p : nearbyPoints)
        pointSet.insert(p);

    int smallestTraceLineX = p.x, largestTraceLineX = p.x;
    int smallestTraceLineY = p.y, largestTraceLineY = p.y;

    // Trace the line forward then backward.
    int directions[2] = {1, -1};
    for (auto direction : directions) {
        p = highestDensityPoint;
        int lineTracingStepSize = smallLineTracingStepSize;
        int searchRadius = smallSearchRadius;
        int maxX = readLen;
        int maxY = trimmedRefSeq.length();
        bool failed;
        while (true) {
            Point newP = p;
            int step = direction * lineTracingStepSize;
            newP.x += step;
            newP.y += step;
            if (newP.x > maxX || newP.y > maxY)
                p = newP;
            else if (newP.x < 0 || newP.y < 0)
                p = newP;
            else
                p = getHighestDensityPointNearPoint(lineTracingStepSize, newP, cloud, index,
                                                    highestDensityScore, &failed);
                if (failed) {
                    lineTracingStepSize = largeLineTracingStepSize;
                    searchRadius = largeSearchRadius;
                }
                else {
                    lineTracingStepSize = smallLineTracingStepSize;
                    searchRadius = smallSearchRadius;
                }

            if (p.x == -1 || p.y == -1)
                p = newP;
    //        std::cout << "(" << p.x << "," << p.y << ") " << std::flush;
            std::vector<Point> nearbyPoints = radiusSearchAroundPoint(p, searchRadius, cloud, index);
            traceDots.push_back(p);
            for (auto p : nearbyPoints)
                pointSet.insert(p);

            smallestTraceLineX = std::min(p.x, smallestTraceLineX);
            smallestTraceLineY = std::min(p.y, smallestTraceLineY);
            largestTraceLineX = std::max(p.x, largestTraceLineX);
            largestTraceLineY = std::max(p.y, largestTraceLineY);

            if (direction == 1 && (p.x > maxX || p.y > maxY))
                break;
            if (direction == -1 && (p.x < 0 || p.y < 0))
                break;
        }
    }

    // This figures out whether the read alignment is contained within a large contig. If so, then
    // Unicycler won't try as hard to align it well, because it won't be informative for bridging.
    bool containedRead = (smallestTraceLineY > 1000 && trimmedRefLen - largestTraceLineY > 1000);

//    std::cout << readName << ", " << refName << std::endl;
//    std::cout << "Trace line X: " << smallestTraceLineX << " to " << largestTraceLineX << std::endl;
//    std::cout << "Trace line Y: " << smallestTraceLineY << " to " << largestTraceLineY << std::endl;
//    std::cout << "Contained: " << containedRead << std::endl;
//    std::cout << std::endl;


    if (verbosity > 3) {  // only at very high verbosities
        std::ofstream traceDotsFile;
        std::string filename = readName + readStrand + "_" + refName + "_trace_dots.tsv";
        traceDotsFile.open(filename);
        for (auto d : traceDots)
            traceDotsFile << d.x << "\t" << d.y << "\n";
        traceDotsFile.close();
        output += "R:trace.dots <- read_delim(\"" + filename + "\", \"\t\", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)\n";
        
        std::ofstream filteredDataFile;
        filename = readName + readStrand + "_" + refName + "_filtered_data.tsv";
        filteredDataFile.open(filename);
        for (auto d : pointSet)
            filteredDataFile << d.x << "\t" << d.y << "\n";
        filteredDataFile.close();
        output += "R:filtered.data <- read_delim(\"" + filename + "\", \"\t\", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)\n";

        output += "R:dot.plot.1(all.points)\n";
        output += "R:dot.plot.2(all.points, trace.dots)\n";
        output += "R:dot.plot.3(all.points, filtered.data, trace.dots)\n";
    }


    // Now we can do a Seqan alignment around the points we've collected!

    typedef Seed<Simple> TSeed;
    typedef SeedSet<TSeed> TSeedSet;

    String<TSeed> seeds;
    for (auto p : pointSet)
        appendValue(seeds, TSeed(p.x, p.y, kSize));

//    std::cout << "before local chaining: " << length(seeds) << std::endl;

    TSeedSet seedSet;
    for (unsigned i = 0; i < length(seeds); ++i) {
        if (!addSeed(seedSet, seeds[i], 2, Merge()))
            addSeed(seedSet, seeds[i], Single());
    }

//    std::cout << "after local chaining: " << length(seedSet) << std::endl;
//    std::cout << "Resulting seeds.\n";
//    typedef Iterator<TSeedSet>::Type TIter;
//    for (TIter it = begin(seedSet, Standard());
//         it != end(seedSet, Standard()); ++it)
//        std::cout << "(" << beginPositionH(*it) << ", " << endPositionH(*it)
//                  << ", " << beginPositionV(*it) << ", " << endPositionV(*it)
//                  << ", " << lowerDiagonal(*it) << ", " << upperDiagonal(*it)
//                  << ")\n";

    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());

//    std::cout << "after global chaining: " << length(seedChain) << std::endl;
//    for (unsigned i = 0; i < length(seedChain); ++i) {
//        std::cout << "(" << beginPositionH(seedChain[i]) << ", " << endPositionH(seedChain[i])
//                  << ", " << beginPositionV(seedChain[i]) << ", " << endPositionV(seedChain[i])
//                  << ", " << lowerDiagonal(seedChain[i]) << ", " << upperDiagonal(seedChain[i])
//                  << ")\n";
//    }


    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), *readSeq);
    assignSource(row(alignment, 1), trimmedRefSeq);
    AlignConfig<true, true, true, true> alignConfig;
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    int bandSize = LEVEL_0_BAND_SIZE;
    if (sensitivityLevel == 1)
        bandSize = LEVEL_1_BAND_SIZE;
    else if (sensitivityLevel == 2)
        bandSize = LEVEL_2_BAND_SIZE;
    else if (sensitivityLevel == 3)
        bandSize = LEVEL_3_BAND_SIZE;

    // First try a fast alignment with a small band.
    ScoredAlignment * sgAlignment;
    try {
        bandedChainAlignment(alignment, seedChain, scoringScheme, alignConfig, bandSize);
        std::string signedReadName = readName + readStrand;
        sgAlignment = new ScoredAlignment(alignment, signedReadName, refName, readLen, refLen,
                                          refStart, startTime, bandSize, false, false, false,
                                          scoringScheme);
    }
    catch (...) {
        sgAlignment = 0;
    }

    alignments.push_back(sgAlignment);
    return alignments;
}


void addKmerPointsToNanoflann(PointCloud & cloud, std::vector<CommonKmer> & commonKmers) {
    cloud.pts.resize(commonKmers.size());
    for (size_t i = 0; i < commonKmers.size(); ++i) {
        cloud.pts[i].x = commonKmers[i].m_hPosition;
        cloud.pts[i].y = commonKmers[i].m_vPosition;
    }
}


std::vector<Point> radiusSearchAroundPoint(Point point, int radius, PointCloud & cloud,
                                           my_kd_tree_t & index) {
    std::vector<Point> points;
    nanoflann::SearchParams params;
    std::vector<std::pair<size_t,int> > ret_matches;
    const int query_pt[2] = {point.x, point.y};
    index.radiusSearch(query_pt, radius, ret_matches, params);
    for (auto i : ret_matches)
        points.push_back(cloud.pts[i.first]);
    return points;
}

std::vector<Point> getPointsInHighestDensityRegion(int searchRadius, std::string & trimmedRefSeq,
                                                   std::string * readSeq, PointCloud & cloud,
                                                   my_kd_tree_t & index) {

    int xStepCount = int(ceil(readSeq->length() / double(searchRadius)));
    int yStepCount = int(ceil(trimmedRefSeq.length() / double(searchRadius)));
    double xStepSize = double(readSeq->length()) / xStepCount;
    double yStepSize = double(trimmedRefSeq.length()) / yStepCount;

//    std::cout << std::endl;
//    std::cout << "xStepCount=" << xStepCount << std::endl;
//    std::cout << "yStepCount=" << yStepCount << std::endl;
//    std::cout << "xStepSize=" << xStepSize << std::endl;
//    std::cout << "yStepSize=" << yStepSize << std::endl << std::flush;

    nanoflann::SearchParams params;
    double highestDensity = 0.0;
    std::vector<Point> pointsInHighestDensity;

    for (int i = 0; i <= xStepCount; ++i) {
        int xCentre = int(0.5 + i * xStepSize);

        for (int j = 0; j <= yStepCount; ++j) {
            int yCentre = int(0.5 + j * yStepSize);

//            std::cout << "xCentre=" << xCentre << ", yCentre=" << yCentre << std::endl << std::flush;

            const int query_pt[2] = {xCentre, yCentre};

            std::vector<std::pair<size_t,int> > ret_matches;
            const size_t nMatches = index.radiusSearch(query_pt, searchRadius, ret_matches, params);
            double density = double(nMatches);

//            std::cout << "unadjusted density=" << density << std::endl << std::flush;

            // If the search region is on the edge, increase the density (because the region has
            // less area). It would be technically correct to double the density on the edges,
            // but this biases the density peak towards edges, so we only use 1.5 instead.
            if (i == 0 || i == xStepCount) density *= 1.5;
            if (j == 0 || j == yStepCount) density *= 1.5;

//            std::cout << "adjusted density=" << density << std::endl << std::flush;

            if (density > highestDensity) {
//                std::cout << "NEW BEST!" << std::endl << std::flush;
                highestDensity = density;
                pointsInHighestDensity.clear();
                for (auto k : ret_matches)
                    pointsInHighestDensity.push_back(cloud.pts[k.first]);
            }
        }
    }
    return pointsInHighestDensity;
}

Point getHighestDensityPoint(int densityRadius, PointCloud & cloud, my_kd_tree_t & index,
                             std::string & trimmedRefSeq, std::string * readSeq,
                             double * highestDensityScore) {

    std::vector<Point> points = getPointsInHighestDensityRegion(densityRadius * 2, trimmedRefSeq,
                                                                readSeq, cloud, index);
    Point highestDensityPoint = points[0];
    *highestDensityScore = 0.0;

    for (auto point : points) {
        double densityScore = getPointDensityScore(densityRadius, point, cloud, index);
        if (densityScore > *highestDensityScore) {
            *highestDensityScore = densityScore;
            highestDensityPoint = point;
        }
    }
    return highestDensityPoint;
}

Point getHighestDensityPointNearPoint(int densityRadius, Point centre, PointCloud & cloud,
                                      my_kd_tree_t & index, double highestDensityScore,
                                      bool * failed) {
    std::vector<Point> points = radiusSearchAroundPoint(centre, densityRadius, cloud, index);
    if (points.size() == 0)
        return {-1, -1};
    Point highestDensityPoint = centre;
    *failed = true;
    double bestDensityScore = highestDensityScore / 10.0;

    for (auto point : points) {
        double densityScore = getPointDensityScore(densityRadius, point, cloud, index);

        // Boost the density score for points near the centre.
        int distanceFromCentre = abs(point.x - centre.x) + abs(point.y - centre.y);
        densityScore *= (1.0 + ((densityRadius - distanceFromCentre) / densityRadius));

        if (densityScore > bestDensityScore) {
            bestDensityScore = densityScore;
            highestDensityPoint = point;
            *failed = false;
        }
    }

//    std::cout << "Starting point: " << centre.x << "," << centre.y << "\n";
//    std::cout << "Highest density point: " << highestDensityPoint.x << "," << highestDensityPoint.y << "\n";
//    std::cout << "Density score: " << bestDensityScore << "\n";
//    std::cout << "\n";

    return highestDensityPoint;
}


double getPointDensityScore(int densityRadius, Point p, PointCloud & cloud, my_kd_tree_t & index) {
    std::vector<Point> neighbourPoints = radiusSearchAroundPoint(p, densityRadius, cloud, index);
    double densityScore = 0.0;
    for (auto neighbourPoint : neighbourPoints) {
        int xDiff = neighbourPoint.x - p.x;
        int yDiff = neighbourPoint.y - p.y;
        if (xDiff + yDiff > 0)
            densityScore += 1.0 / (abs(xDiff-yDiff) + 1.0);
    }
    return densityScore;
}

//
//
//    // Choose settings based on the sensitivity.
//    int badLineCountSingleAlignment;
//    int badLineCountEntireRead;
//    int badLineCountPartialRead;
//    if (extraSensitive) {
//        badLineCountSingleAlignment = BAD_LINE_COUNT_SINGLE_ALIGNMENT_EXTRA_SENSITIVE;
//        badLineCountEntireRead = BAD_LINE_COUNT_ENTIRE_READ_EXTRA_SENSITIVE;
//        badLineCountPartialRead = BAD_LINE_COUNT_PARTIAL_READ_EXTRA_SENSITIVE;
//    }
//    else {
//        badLineCountSingleAlignment = BAD_LINE_COUNT_SINGLE_ALIGNMENT;
//        badLineCountEntireRead = BAD_LINE_COUNT_ENTIRE_READ;
//        badLineCountPartialRead = BAD_LINE_COUNT_PARTIAL_READ;
//    }
//
//
//    // Change the read name and sequence to C++ strings.
//    std::string readName(readNameC);
//    std::string posReadName = readName + "+";
//    std::string negReadName = readName + "-";
//    std::string posReadSeq(readSeqC);
//    std::string negReadSeq = getReverseComplement(posReadSeq);
//    int readLength = posReadSeq.length();
//
//
//    // Create a CommonKmerSet for the read (both forward and reverse complement) and every reference.
//    std::vector<std::string> referenceNames = refKmerPositions->getAllNames();
//    std::vector<CommonKmerSet *> commonKmerSets;
//    for (size_t i = 0; i < referenceNames.size(); ++i) {
//        std::string refName = referenceNames[i];
//        int refLength = refKmerPositions->getLength(refName);
//        CommonKmerSet * forwardCommonKmerSet = new CommonKmerSet(posReadName, refName, readLength, refLength, expectedSlope, &readKmerPositions, refKmerPositions, kSize);
//        commonKmerSets.push_back(forwardCommonKmerSet);
//        CommonKmerSet * reverseCommonKmerSet = new CommonKmerSet(negReadName, refName, readLength, refLength, expectedSlope, &readKmerPositions, refKmerPositions, kSize);
//        commonKmerSets.push_back(reverseCommonKmerSet);
//    }
//
//    if (verbosity > 2)
//        output += "Seqan alignment attempts (using expected slope of " + std::to_string(expectedSlope) + ")\n";
//
//    // We now extract alignment lines and perform alignments until we have had too many failures.
//    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);
//    std::vector<AlignmentLine *> alignmentLines;
//    std::vector<ScoredAlignment *> allAlignments;
//    std::vector<ScoredAlignment *> goodAlignments;
//    int badAlignmentCount = 0;
//    bool oneAlignmentWholeRead = false;
//    bool entireReadAligned = false;
//    bool needMoreAlignments = true;
//
//    while (needMoreAlignments) {
//        // Extract an alignment line from around the highest scoring point.
//        CommonKmerSet * highestScoringSet = getHighestScoringSet(commonKmerSets);
//        if (highestScoringSet == 0)
//            break;
//        std::string readName = highestScoringSet->m_readName;
//        std::string refName = highestScoringSet->m_refName;
//        AlignmentLine * line = highestScoringSet->extractAlignmentLine();
//
//        // Check to see if the alignment line is good.
//        if (line == 0) {
//            ++badAlignmentCount;
//            if (verbosity > 2)
//                output += "  line: " + refName + ", none, BAD\n";
//        }
//        else if (!line->buildSeedChain(MIN_POINT_COUNT, MIN_ALIGNMENT_LENGTH, kSize)) {
//            ++badAlignmentCount;
//            if (verbosity > 2)
//                output += "  line: " + line->getDescriptiveString() + ", no seed chain, BAD\n";
//            delete line;
//        }
//        else {
//            // If the code got here, then we should have an alignment line with a seed chain ready to
//            // go, so let's perform the alignment!
//            if (verbosity > 2)
//                output += "  line: " + line->getDescriptiveString() + ", GOOD\n";
//            std::string * readSeq = readKmerPositions.getSequence(readName);
//            std::string * refSeq = refKmerPositions->getSequence(refName);
//            ScoredAlignment * alignment = semiGlobalAlignmentOneLine(readName, refName, readSeq, refSeq,
//                                                                     line, verbosity, output, scoringScheme);
//            alignmentLines.push_back(line);
//            allAlignments.push_back(alignment);
//
//            // Check to see if the alignment failed or if it scored too low.
//            if (alignment == 0) {
//                ++badAlignmentCount;
//                if (verbosity > 2)
//                    output += "    alignment: failed, BAD\n";
//            }
//            else if (alignment->m_scaledScore < lowScoreThreshold) {
//                ++badAlignmentCount;
//                if (verbosity > 2)
//                    output += "    alignment: " + alignment->getShortDisplayString() + ", BAD\n";
//            }
//
//            // If the alignment is good, we add it to the results.
//            else {
//                goodAlignments.push_back(alignment);
//                if (verbosity > 2)
//                    output += "    alignment: " + alignment->getShortDisplayString() + ", GOOD\n";
//                if (alignment->getReadAlignmentLength() == readLength)
//                    oneAlignmentWholeRead = true;
//            }
//
//            // Finally, we want to check whether this alignment line is near any of the previous
//            // alignment lines.
//            AlignmentLine * bestNearbyLine = 0;
//            ScoredAlignment * bestNearbyLineAlignment = 0;
//            for (size_t i = 0; i < alignmentLines.size() - 1; ++i) {
//                AlignmentLine * previousLine = alignmentLines[i];
//                ScoredAlignment * previousLineAlignment = allAlignments[i];
//
//                if (previousLine->isNear(line)) {
//
//                    // If this is the first nearby line found, or...
//                    if (bestNearbyLine == 0 ||
//
//                        // this is the first nearby line found with an alignment, or...
//                        (bestNearbyLineAlignment == 0 && previousLineAlignment != 0) ||
//
//                        // this nearby line has a better score than our previous best...
//                        (bestNearbyLineAlignment != 0 && previousLineAlignment != 0 &&
//                         previousLineAlignment->m_scaledScore > bestNearbyLineAlignment->m_scaledScore) )
//                    {
//                        bestNearbyLine = previousLine;
//                        bestNearbyLineAlignment = previousLineAlignment;
//                    }
//                }
//            }
//
//            // If the alignment line was indeed near another one, we try merging the two lines to
//            // make a new line and aligning to that. This is good for cases of long alignments that
//            // aren't well captured by a single line.
//            if (bestNearbyLine != 0) {
//                if (verbosity > 2) {
//                    output += "  merging lines: 1) " + line->getDescriptiveString() + "\n";
//                    output += "                 2) " + bestNearbyLine->getDescriptiveString() + "\n";
//                }
//                AlignmentLine * mergedLine = new AlignmentLine(line, bestNearbyLine);
//                if (!mergedLine->buildSeedChain(MIN_POINT_COUNT, MIN_ALIGNMENT_LENGTH, kSize)) {
//                    if (verbosity > 2)
//                        output += "                 failed\n";
//                    delete mergedLine;
//                }
//                else {
//                    if (verbosity > 2)
//                        output += "                 result: " + mergedLine->getDescriptiveString() + "\n";
//                    ScoredAlignment * mergedLineAlignment = semiGlobalAlignmentOneLine(readName, refName, readSeq, refSeq,
//                                                                                       mergedLine, verbosity, output, scoringScheme);
//                    alignmentLines.push_back(mergedLine);
//                    allAlignments.push_back(mergedLineAlignment);
//                    if (mergedLineAlignment == 0) {
//                        if (verbosity > 2)
//                            output += "                 merged line alignment failed\n";
//                    }
//                    else if (mergedLineAlignment->m_scaledScore < lowScoreThreshold) {
//                        if (verbosity > 2)
//                            output += "                 merged line alignment: " + mergedLineAlignment->getShortDisplayString() + ", BAD\n";
//                    }
//                    else {
//                        goodAlignments.push_back(mergedLineAlignment);
//                        if (verbosity > 2)
//                            output += "                 merged line alignment: " + mergedLineAlignment->getShortDisplayString() + ", GOOD\n";
//                    }
//                }
//            }
//        }
//
//        if (oneAlignmentWholeRead)
//            needMoreAlignments = (badAlignmentCount < badLineCountSingleAlignment);
//        else {
//            entireReadAligned = (fractionOfReadAligned(goodAlignments) == 1.0);
//            if (entireReadAligned)
//                needMoreAlignments = (badAlignmentCount < badLineCountEntireRead);
//            else
//                needMoreAlignments = (badAlignmentCount < badLineCountPartialRead);
//        }
//    }
//
//    // Either all alignments or only good alignments are returned, depending on a parameter.
//
//    if (returnBad)
//        returnedAlignments = &allAlignments;
//    else
//        returnedAlignments = &goodAlignments;
//
//
//    // Clean up.
//    for (size_t i = 0; i < commonKmerSets.size(); ++i)
//        delete commonKmerSets[i];
//    for (size_t i = 0; i < allAlignments.size(); ++i) {
//        if (allAlignments[i] != 0)
//            delete allAlignments[i];
//        delete alignmentLines[i];
//    }
//
//}


//
//
// // Runs an alignment using Seqan between one read and one reference along one line.
// // It starts with a smallish band size (fast) and works up to larger ones to see if they improve
// // the alignment.
//ScoredAlignment * semiGlobalAlignmentOneLine(std::string & readName, std::string & refName,
//                                             std::string * readSeq, std::string * refSeq,
//                                             AlignmentLine * line, int verbosity, std::string & output,
//                                             Score<int, Simple> & scoringScheme) {
//    long long startTime = getTime();
//
//    int trimmedRefLength = line->m_trimmedRefEnd - line->m_trimmedRefStart;
//    std::string trimmedRefSeq = refSeq->substr(line->m_trimmedRefStart, trimmedRefLength);
//
//    Dna5String readSeqSeqan(*readSeq);
//    Dna5String refSeqSeqan(trimmedRefSeq);
//    int readLength = readSeq->length();
//
//    int bandSize = STARTING_BAND_SIZE;
//    ScoredAlignment * bestAlignment = 0;
//    double bestAlignmentScore = std::numeric_limits<double>::min();
//
//    // We perform the alignment with increasing band sizes until the score stops improving or we
//    // reach the max band size.
//    while (true) {
//        ScoredAlignment * alignment = semiGlobalAlignmentOneLineOneBand(readName, refName,
//                                                                        readSeqSeqan, readLength,
//                                                                        refSeqSeqan, trimmedRefLength,
//                                                                        line, bandSize, verbosity,
//                                                                        output, scoringScheme);
//        if (alignment != 0) {
//            double alignmentScore = alignment->m_scaledScore;
//            bool badScore = (alignmentScore <= bestAlignmentScore);
//            bool tooShort = (alignment->getReadAlignmentLength() < MIN_ALIGNMENT_LENGTH) ||
//                            (alignment->getRefAlignmentLength() < MIN_ALIGNMENT_LENGTH);
//            if (badScore || tooShort) {
//                delete alignment;
//                break;
//            }
//            else {
//                if (bestAlignment != 0)
//                    delete bestAlignment;
//                bestAlignment = alignment;
//                bestAlignmentScore = alignmentScore;
//            }
//        }
//        bandSize *= 2;
//        if (bandSize > MAX_BAND_SIZE)
//            break;
//    }
//
//    if (bestAlignment != 0)
//        bestAlignment->m_milliseconds = getTime() - startTime;
//
//    return bestAlignment;
//}
//
//
//
//
//
//
//// This function, given a line, will search for semi-global alignments around that line. The
//// bandSize parameter specifies how far of an area around the line is searched.
//ScoredAlignment * semiGlobalAlignmentOneLineOneBand(std::string & readName, std::string & refName,
//                                                        Dna5String & readSeq, int readLen,
//                                                        Dna5String & refSeq, int refLen,
//                                                        AlignmentLine * line, int bandSize,
//                                                        int verbosity, std::string & output,
//                                                        Score<int, Simple> & scoringScheme) {
//    long long startTime = getTime();
//
//    // I encountered a Seqan crash when the band size exceeded the sequence length, so don't let
//    // that happen.
//    int shortestSeqLen = std::min(readLen, refLen);
//    if (bandSize > shortestSeqLen)
//        bandSize = shortestSeqLen;
//
//    // The reference sequence here is the trimmed reference sequence, not the whole reference
//    // sequence. But the seed chain was made using the same offset as the trimming, so everything
//    // should line up nicely (no offset adjustment needed).
//
//    Align<Dna5String, ArrayGaps> alignment;
//    resize(rows(alignment), 2);
//    assignSource(row(alignment, 0), readSeq);
//    assignSource(row(alignment, 1), refSeq);
//    AlignConfig<true, true, true, true> alignConfig;
//
//    ScoredAlignment * sgAlignment;
//    try {
//        bandedChainAlignment(alignment, line->m_bridgedSeedChain, scoringScheme, alignConfig,
//                             bandSize);
//        sgAlignment = new ScoredAlignment(alignment, readName, refName, readLen, refLen,
//                                          line->m_trimmedRefStart, startTime, bandSize, false,
//                                          false, false, scoringScheme);
//    }
//    catch (...) {
//        if (verbosity > 2)
//            output += "  Alignment failed, bandwidth = " + std::to_string(bandSize) + "\n";
//        sgAlignment = 0;
//    }
//
//    return sgAlignment;
//}
//
//
//
//
//// This function is used to conduct a short alignment for the sake of extending a GraphMap
//// alignment.
//char * startExtensionAlignment(char * read, char * ref,
//                               int matchScore, int mismatchScore, int gapOpenScore,
//                               int gapExtensionScore) {
//    long long startTime = getTime();
//    std::string output;
//
//    Dna5String sequenceH = read;
//    Dna5String sequenceV = ref;
//    std::string readName = "";
//    std::string refName = "";
//
//    Align<Dna5String, ArrayGaps> alignment;
//    resize(rows(alignment), 2);
//    assignSource(row(alignment, 0), sequenceH);
//    assignSource(row(alignment, 1), sequenceV);
//    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);
//
//    // The only free gaps are at the start of ref (the reference sequence).
//    AlignConfig<false, true, false, false> alignConfig;
//    globalAlignment(alignment, scoringScheme, alignConfig);
//
//    ScoredAlignment startAlignment(alignment, readName, refName, length(read), length(ref),
//                                   0, startTime, 0, false, true, true, scoringScheme);
//    return cppStringToCString(startAlignment.getFullString());
//}
//
//
//
//// This function is used to conduct a short alignment for the sake of extending a GraphMap
//// alignment.
//char * endExtensionAlignment(char * read, char * ref,
//                             int matchScore, int mismatchScore, int gapOpenScore,
//                             int gapExtensionScore) {
//    long long startTime = getTime();
//    std::string output;
//
//    Dna5String sequenceH = read;
//    Dna5String sequenceV = ref;
//    std::string readName = "";
//    std::string refName = "";
//
//    Align<Dna5String, ArrayGaps> alignment;
//    resize(rows(alignment), 2);
//    assignSource(row(alignment, 0), sequenceH);
//    assignSource(row(alignment, 1), sequenceV);
//    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);
//
//    // The only free gaps are at the end of ref (the reference sequence).
//    AlignConfig<false, false, true, false> alignConfig;
//    globalAlignment(alignment, scoringScheme, alignConfig);
//
//    ScoredAlignment endAlignment(alignment, readName, refName, length(read), length(ref),
//                                 0, startTime, 0, true, false, false, scoringScheme);
//    return cppStringToCString(endAlignment.getFullString());
//}
//


double fractionOfReadAligned(std::vector<ScoredAlignment *> & alignments) {
    if (alignments.size() == 0)
        return true;
    std::vector<std::pair<int, int> > ranges;
    for (size_t i = 0; i < alignments.size(); ++i) {
        ScoredAlignment * alignment = alignments[i];
        int start, end;
        if (alignment->isRevComp()) {
            start = alignment->m_readLength - alignment->m_readEndPos;
            end = alignment->m_readLength - alignment->m_readStartPos;
        }
        else {
            start = alignment->m_readStartPos;
            end = alignment->m_readEndPos;
        }
        ranges.push_back(std::pair<int, int>(start, end));
    }
    std::vector<std::pair<int, int> > simplifiedRanges = simplifyRanges(ranges);
    int alignedLength = 0;
    for (size_t i = 0; i < simplifiedRanges.size(); ++i)
        alignedLength += simplifiedRanges[i].second - simplifiedRanges[i].first;
    return double(alignedLength) / alignments[0]->m_readLength;
}


std::pair<int, int> getRefRange(int refStart, int refEnd, int refLen,
                                int readStart, int readEnd, int readLen, bool posStrand) {
    int halfReadLen = 1 + readLen / 2;
    int readBasesBeforeStart = readStart;
    int readBasesAfterEnd = readLen - readEnd;
    if (!posStrand)
        std::swap(readBasesBeforeStart, readBasesAfterEnd);
    int newRefStart = refStart - readBasesBeforeStart - halfReadLen;
    int newRefEnd = refEnd + readBasesAfterEnd + halfReadLen;
    newRefStart = std::max(0, newRefStart);
    newRefEnd = std::min(refLen, newRefEnd);
    return std::pair<int, int>(newRefStart, newRefEnd);
}


std::vector<StartEndRange> simplifyRanges(std::vector<StartEndRange> & ranges) {
    std::sort(ranges.begin(),ranges.end());
    std::vector<std::pair<int, int> > simplifiedRanges;
    std::vector<std::pair<int, int> >::iterator it = ranges.begin();
    std::pair<int,int> current = *(it)++;
    while (it != ranges.end()){
       if (current.second >= it->first){
           current.second = std::max(current.second, it->second); 
       } else {
           simplifiedRanges.push_back(current);
           current = *(it);
       }
       it++;
    }
    simplifiedRanges.push_back(current);
    return simplifiedRanges;
}

//
//// Given a vector of CommonKmerSet pointers, this function returns the one with the highest score.
//CommonKmerSet * getHighestScoringSet(std::vector<CommonKmerSet *> & commonKmerSets) {
//    double allSetsMaxScore = std::numeric_limits<double>::min();
//    CommonKmerSet * bestSet = 0;
//    for (size_t i = 0; i < commonKmerSets.size(); ++i) {
//        double setMaxScore = commonKmerSets[i]->m_maxScore;
//        if (setMaxScore > allSetsMaxScore) {
//            allSetsMaxScore = setMaxScore;
//            bestSet = commonKmerSets[i];
//        }
//    }
//
//    // The best set has to at least reach the low threshold.
//    if (allSetsMaxScore < MIN_LINE_SCORE)
//        return 0;
//
//    return bestSet;
//}



