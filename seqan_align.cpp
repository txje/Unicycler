#include <stdio.h>
#include <chrono>
#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seeds.h>
#include <string>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <utility>
#include <iterator>
#include <mutex>
#include <limits>

#define KMER_SIZE 5

using namespace seqan;
using namespace std::chrono;

extern "C" {

typedef Dna5String TSequence;
typedef Align<TSequence, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;
typedef Seed<Simple> TSeed;
typedef SeedSet<TSeed> TSeedSet;
typedef Row<TAlign>::Type TRow;
typedef Iterator<TRow>::Type TRowIterator;
typedef std::unordered_map<std::string, std::vector<int> > KmerPosMap;

enum CigarType {MATCH, INSERTION, DELETION, CLIP, NOTHING};

class CommonKmer
{
public:
    CommonKmer(std::string sequence, int hPosition, int vPosition, double angle);
    static double getRotationAngle(double slope) {return -atan(slope);}

    std::string m_sequence;
    int m_hPosition;
    int m_vPosition;
    double m_rotatedHPosition;
    double m_rotatedVPosition;
    double m_bandArea;
    double m_score; // Scaled kmer density
};

// KmerPositions is a class that holds maps of k-mer positions for named sequences. It exists so we
// don't have to repeatedly find the same k-mer sets over and over. Whenever it makes a k-mer map,
// it stores it for later. On destruction it deletes all of its k-mer maps.
class KmerPositions
{
public:
    KmerPositions() {}
    ~KmerPositions();
    void addPositions(char * nameC, char * sequenceC);
    void deletePositions(std::string & name);
    KmerPosMap * getKmerPositions(std::string & name);

    std::unordered_map<std::string, KmerPosMap *> m_kmerPositions;
};


// These are the functions that will be called by the Python script.
char * bandedSemiGlobalAlignment(char * read, int readLen,
                                 char * ref, int refLen, int refOffset,
                                 double slope, double intercept, int bandSize,
                                 int verbosity,
                                 int matchScore, int mismatchScore,
                                 int gapOpenScore, int gapExtensionScore,
                                 char * kmerLocations);
char * findAlignmentLines(char * readNameC, char * refNameC, int readLength, int refLength,
                          double expectedSlope, int verbosity, KmerPositions * kmerPositions);
char * startExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
                               int matchScore, int mismatchScore, int gapOpenScore,
                               int gapExtensionScore);
char * endExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
                             int matchScore, int mismatchScore, int gapOpenScore,
                             int gapExtensionScore);
void free_c_string(char * p) {free(p);}
KmerPositions * newKmerPositions() {return new KmerPositions();}
void addKmerPositions(KmerPositions * kmerPositions, char * nameC, char * sequenceC);
void deleteKmerPositions(KmerPositions * kmerPositions, char * nameC);
void deleteAllKmerPositions(KmerPositions * kmerPositions) {delete kmerPositions;}


// These functions are internal to this C++ code.
void addBridgingSeed(TSeedSet & seedSet, int hStart, int vStart, int hEnd, int vEnd);
std::vector<CommonKmer> getCommonKmers(std::string & readName, std::string & refName,
                                       double expectedSlope, int verbosity, std::string & output,
                                       KmerPositions * kmerPositions);
char * turnAlignmentIntoDescriptiveString(Align<Dna5String, ArrayGaps> * alignment,
                                          int refOffset, long long startTime, std::string output,
                                          bool startImmediately, bool goToEnd);
CigarType getCigarType(char b1, char b2, bool alignmentStarted);
std::string getCigarPart(CigarType type, int length);
char * cppStringToCString(std::string cpp_string);
std::string vectorToString(std::vector<int> * v);
void fixSeedChainToLine(String<TSeed> * seedChain, double bandSize);
void getSeedChainSlopeAndIntercept(String<TSeed> * seedChain, int firstI, int lastI,
                                   double * slope, double * intercept);
void getSlopeAndIntercept(int hStart, int hEnd, int vStart, int vEnd,
                                   double * slope, double * intercept);
double getMedian(std::vector<double> & v);
void printKmerSize(int kmerSize, int locationCount, std::string & output);
double getLineLength(double x, double y, double slope, double xSize, double ySize);
void linearRegression(std::vector<CommonKmer> & pts, double * slope, double * intercept);
void parseKmerLocationsFromString(std::string & str, std::vector<int> & v1, std::vector<int> & v2,
                                  int refOffset);
std::string getKmerTable(std::vector<CommonKmer> & commonKmers);
std::string getSeedChainTable(String<TSeed> & seedChain);
void addSeedMerge(TSeedSet & seedSet, TSeed & seed);
void getMeanAndStDev(std::vector<double> & v, double & mean, double & stdev);



// This function, given a line, will search for semi-global alignments around that line. The
// bandSize parameter specifies how far of an area around the line is searched. It is generally
// much faster than exhaustiveSemiGlobalAlignment, though it may not find the optimal alignment.
// A lower bandSize is faster with a larger chance of missing the optimal alignment.
char * bandedSemiGlobalAlignment(char * read, int readLen,
                                 char * ref, int refLen, int refOffset,
                                 double slope, double intercept, int bandSize,
                                 int verbosity,
                                 int matchScore, int mismatchScore,
                                 int gapOpenScore, int gapExtensionScore,
                                 char * kmerLocations)
{
    std::string output = "";
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    // Extreme slope values will not work.
    if (slope > 1.5)
        return cppStringToCString(output + ";Failed: slope too large");
    if (slope < 0.5)
        return cppStringToCString(output + ";Failed: slope too small");

    TSequence sequenceH = read;
    TSequence sequenceV = ref;

    std::string kmerLocationsStr(kmerLocations);
    std::vector<int> readKmerLocations;
    std::vector<int> refKmerLocations;
    parseKmerLocationsFromString(kmerLocationsStr, readKmerLocations, refKmerLocations, refOffset);

    // Build a Seqan seed set using our common k-mers.
    TSeedSet seedSet;
    for (int i = 0; i < readKmerLocations.size(); ++i)
    {
        TSeed seed(readKmerLocations[i], refKmerLocations[i], KMER_SIZE);
        addSeedMerge(seedSet, seed);
    }

    // We now get a Seqan global chain of the seeds.
    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());
    int seedsInChain = length(seedChain);
    if (seedsInChain == 0)
        return cppStringToCString(output + ";Failed: no global seed chain");

    if (verbosity > 4)
    {
        output += "  Globally chained seeds before bridging\n";
        output += getSeedChainTable(seedChain);
    }

    // Now we create a new seed chain will all of the gaps bridged. This will help keep alignment
    // in a narrow band, even when the seeds are spaced apart.
    TSeedSet bridgedSeedSet;

    // Create a seed bridge for the start of the chain by following the slope backwards from the
    // first seed.
    TSeed firstSeed = seedChain[0];
    int firstSeedHStart = beginPositionH(firstSeed);
    int firstSeedVStart = beginPositionV(firstSeed);
    if (firstSeedHStart > 0 && firstSeedVStart > 0)
    {
        double vPosAtStartOfH = firstSeedVStart - (slope * firstSeedHStart);
        if (vPosAtStartOfH >= 0.0)
            addBridgingSeed(bridgedSeedSet, 0, std::round(vPosAtStartOfH), firstSeedHStart, firstSeedVStart);
        else
            addBridgingSeed(bridgedSeedSet, std::round(-vPosAtStartOfH / slope), 0, firstSeedHStart, firstSeedVStart);
    }
    addSeedMerge(bridgedSeedSet, firstSeed);
    
    // Fill in any gaps in the middle of the seed chain.
    for (int i = 1; i < seedsInChain; ++i)
    {
        TSeed seed1 = seedChain[i-1];
        TSeed seed2 = seedChain[i];
        int seed1HEnd = endPositionH(seedChain[i-1]);
        int seed1VEnd = endPositionV(seedChain[i-1]);
        int seed2HStart = beginPositionH(seedChain[i]);
        int seed2VStart = beginPositionV(seedChain[i]);
        addBridgingSeed(bridgedSeedSet, seed1HEnd, seed1VEnd, seed2HStart, seed2VStart);
        addSeedMerge(bridgedSeedSet, seed2);
    }

    // Create a seed bridge for the end of the chain by following the slope forwards from the
    // last seed.
    TSeed lastSeed = seedChain[seedsInChain - 1];
    int lastSeedHEnd = endPositionH(lastSeed);
    int lastSeedVEnd = endPositionV(lastSeed);
    if (lastSeedHEnd < readLen && lastSeedVEnd < refLen)
    {
        double vPosAtStartOfH = lastSeedVEnd - (slope * lastSeedHEnd);
        double vPosAtEndOfH = (slope * readLen) + vPosAtStartOfH;
        if (vPosAtEndOfH <= refLen)
            addBridgingSeed(bridgedSeedSet, lastSeedHEnd, lastSeedVEnd, readLen,
                            std::round(vPosAtEndOfH));
        else
            addBridgingSeed(bridgedSeedSet, lastSeedHEnd, lastSeedVEnd,
                            std::round((refLen - vPosAtStartOfH) / slope), refLen);
    }

    String<TSeed> bridgedSeedChain;
    chainSeedsGlobally(bridgedSeedChain, bridgedSeedSet, SparseChaining());

    if (verbosity > 4)
    {
        output += "  Seed chain after bridging\n";
        output += getSeedChainTable(bridgedSeedChain);
    }

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);
    AlignConfig<true, true, true, true> alignConfig;

    int score = bandedChainAlignment(alignment, bridgedSeedChain, scoringScheme, alignConfig,
                                     bandSize);

    return turnAlignmentIntoDescriptiveString(&alignment, refOffset, startTime, output,
                                              false, false);
}

// This function takes the seed chain which should end at the point hStart, vStart. New seeds will be
// added to the chain to reach the point hEnd, vEnd.
void addBridgingSeed(TSeedSet & seedSet, int hStart, int vStart, int hEnd, int vEnd)
{
    int hDiff = hEnd - hStart;
    int vDiff = vEnd - vStart;

    // If the start and end are against each other, then there's nothing to bridge.
    if (hDiff == 0 || vDiff == 0)
        return;

    TSeed seed(hStart, vStart, hEnd, vEnd);
    addSeedMerge(seedSet, seed);
}

// This function searches for lines in the 2D read-ref space that represent likely semi-global
// alignments.
char * findAlignmentLines(char * readNameC, char * refNameC, int readLength, int refLength,
                          double expectedSlope, int verbosity, KmerPositions * kmerPositions)
{
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    std::string readName(readNameC);
    std::string refName(refNameC);
    std::string output;

    std::vector<CommonKmer> commonKmers = getCommonKmers(readName, refName, expectedSlope,
                                                         verbosity, output, kmerPositions);

    if (commonKmers.size() < 2)
    {
        long long endTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();
        return cppStringToCString(output + ";" + std::to_string(endTime - startTime) + ";Failed: too few common kmers");
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
    int bandSize = 16; // TO DO: MAKE THIS A PARAMETER?
    int halfBandSize = bandSize / 2;
    double maxScore = 0.0;

    // Get the full band length using the middle point of the alignment rectangle.
    double fullBandLength = getLineLength(readLength / 2.0, refLength / 2.0,
                                          expectedSlope, readLength, refLength);

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

    int commonKmerCount = commonKmers.size();
    int startKmerIndex = 0 + (halfBandSize - 1);
    int endKmerIndex = commonKmerCount - (halfBandSize - 1); //One past last

    // Now we loop through the CommonKmer points, calculating their k-mer density (score) along the way!
    for (int i = 0; i < commonKmerCount; ++i)
    {
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
        else
        {
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
            else if (bandStartY <= c1Y && bandEndY >= c1Y && bandEndY <= c3Y)
            {
                double trapezoidArea = (c1Y - bandStartY) * ((bandStartLength + c1BandLength) / 2.0);
                double parallelogramArea = (bandEndY - c1Y) * bandEndLength;
                bandArea = trapezoidArea + parallelogramArea;
            }
            else if (bandStartY >= c1Y && bandStartY <= c3Y && bandEndY >= c3Y)
            {
                double trapezoidArea = (bandEndY - c3Y) * ((c3BandLength + bandEndLength) / 2.0);
                double parallelogramArea = (c3Y - bandStartY) * bandStartLength;
                bandArea = trapezoidArea + parallelogramArea;
            }

            // The final, most complex scenario is when the band start and end span both C1 and C3.
            // This would be unusual, as it would require either a very large band or a very sparse
            // set of CommonKmers.
            else
            {
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

    if (verbosity > 4)
    {
        output += "  Common k-mer positions:\n";
        output += getKmerTable(commonKmers);
        output += "  Max score: " + std::to_string(maxScore) + "\n";
    }

    // Now group all of the line points. For a line group to form, a point must score above the
    // high threshold. The group will continue (in both directions) until the score drops below
    // the low threshold.
    std::vector<std::vector<CommonKmer> > lineGroups;
    double lowThreshold = 2.0; //TO DO: MAKE THIS A PARAMETER?
    double highThreshold = 20.0; //TO DO: MAKE THIS A PARAMETER?
    bool lineInProgress = false;
    for (int i = 0; i < commonKmers.size(); ++i)
    {
        if (lineInProgress)
        {
            if (commonKmers[i].m_score >= lowThreshold)
                lineGroups.back().push_back(commonKmers[i]);
            else // This line group is done.
                lineInProgress = false;
        }
        else if (commonKmers[i].m_score >= highThreshold)
        {
            // It's time to start a new line group!
            lineGroups.push_back(std::vector<CommonKmer>());
            lineInProgress = true;

            // Step backwards to find where the group should start (the first point over the low
            // threshold).
            int groupStartPoint = i;
            while (groupStartPoint >= 0 && commonKmers[groupStartPoint].m_score >= lowThreshold)
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
    double mergeDistance = 100.0; // TO DO: MAKE THIS A PARAMETER?
    std::vector<std::vector<CommonKmer> > mergedLineGroups;
    if (lineGroups.size() > 0)
        mergedLineGroups.push_back(lineGroups[0]);
    for (int i = 1; i < lineGroups.size(); ++i)
    {
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

    // We are only interested in line groups which seem to span their full possible length (i.e.
    // not line groups caused by short, local alignments) and for which the band is reasonably 
    // long (to avoid things like 3 bp alignments).
    double minimumAlignmentLength = 20.0; //TO DO: MAKE THIS A PARAMETER?
    std::vector<std::vector<CommonKmer> > lengthFilteredLineGroups;
    for (int i = 0; i < lineGroups.size(); ++i)
    {
        // Determine the mean position for the line group and use it to calculate the band length.
        int groupCount = lineGroups[i].size();
        double vSum = 0.0, hSum = 0.0;
        for (int j = 0; j < groupCount; ++j)
        {
            hSum += lineGroups[i][j].m_hPosition;
            vSum += lineGroups[i][j].m_vPosition;
        }
        double meanH = hSum / groupCount;
        double meanV = vSum / groupCount;
        double bandLength = getLineLength(meanH, meanV, expectedSlope, readLength, refLength);

        // Exclude alignments which are too short.
        if (bandLength < minimumAlignmentLength)
        {
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
        if (rotatedHStdDev < 0.5 * uniformStdDev)
        {
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

    for (int i = 0; i < lineGroups.size(); ++i)
    {
        std::vector<CommonKmer> * lineGroup = &(lineGroups[i]);

        // Sort by rotated horizontal position so we proceed along the line from start to end. 
        std::sort(lineGroup->begin(), lineGroup->end(), [](const CommonKmer & a, const CommonKmer & b) {
            return a.m_rotatedHPosition < b.m_rotatedHPosition;   
        });

        // Store the max slope for each point in the line group.
        std::vector<double> maxSlopes;
        maxSlopes.reserve(lineGroup->size());
        int groupSize = lineGroup->size();
        for (int j = 0; j < groupSize; ++j)
        {
            CommonKmer * thisPoint = &((*lineGroup)[j]);
            double maxSlope = 0.0;
            for (int k = -steps; k <= steps; ++k)
            {
                int neighbourIndex = j + k;
                if (neighbourIndex >= 0 && neighbourIndex < groupSize)
                {
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
        for (int j = 0; j < groupSize; ++j)
        {
            double maxSlope = maxSlopes[j];
            if (maxSlope <= slopeCutoff)
                fixedLineGroup.push_back((*lineGroup)[j]);
        }
        lineGroups[i] = fixedLineGroup;
    }

    // Remove any line groups with too few points.
    int minPointCount = 4; // TO DO: MAKE A PARAMETER?
    lineGroups.erase(std::remove_if(lineGroups.begin(), lineGroups.end(), 
                                    [&minPointCount](std::vector<CommonKmer> i) {return i.size() < minPointCount;}),
                     lineGroups.end());
    if (verbosity > 4)
        output += "  Number of lines after point count filtering: " + std::to_string(lineGroups.size()) + "\n";

    // Perform a simple least-squares linear regression on each line group.
    // If the slope is too far from 1.0, we throw the line out.
    double minSlope = 0.5; // TO DO: MAKE THIS A PARAMETER?
    double maxSlope = 1.5; // TO DO: MAKE THIS A PARAMETER?
    std::vector<double> slopes;
    std::vector<double> intercepts;
    std::vector<std::vector<CommonKmer> > slopeFilteredLineGroups;
    for (int i = 0; i < lineGroups.size(); ++i)
    {
        double slope, intercept;
        linearRegression(lineGroups[i], &slope, &intercept);
        if (slope >= minSlope && slope <= maxSlope)
        {
            slopeFilteredLineGroups.push_back(lineGroups[i]);
            slopes.push_back(slope);
            intercepts.push_back(intercept);
        }
    }
    lineGroups = slopeFilteredLineGroups;
    if (verbosity > 4)
        output += "  Number of lines after slope filtering: " + std::to_string(lineGroups.size()) + "\n";

    if (lineGroups.size() == 0)
    {
        long long endTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();
        return cppStringToCString(output + ";" + std::to_string(endTime - startTime) + ";Failed: no lines found");
    }
    
    if (verbosity > 3)
        output += "  Lines found:\n";

    std::string linesString;
    for (int i = 0; i < lineGroups.size(); ++i)
    {
        if (linesString.length() > 0)
            linesString += ";";
        linesString += std::to_string(slopes[i]) + "," + std::to_string(intercepts[i]);

        // Add the k-mer locations to the returned string.
        for (int j = 0; j < lineGroups[i].size(); ++j)
            linesString += "," + std::to_string(lineGroups[i][j].m_hPosition) + "," + std::to_string(lineGroups[i][j].m_vPosition);

        if (verbosity > 3)
            output += "    slope = " + std::to_string(slopes[i]) + ", intercept = " + std::to_string(intercepts[i]) + "\n";
        if (verbosity > 4)
            output += getKmerTable(lineGroups[i]);
    }

    long long endTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();
    return cppStringToCString(output + ";" + std::to_string(endTime - startTime) + ";" + linesString);
}


// This function is used to conduct a short alignment for the sake of extending a GraphMap
// alignment.
char * startExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
                               int matchScore, int mismatchScore, int gapOpenScore,
                               int gapExtensionScore)
{
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();
    std::string output;

    TSequence sequenceH = read;
    TSequence sequenceV = ref;

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    // The only free gaps are at the start of ref (the reference sequence).
    AlignConfig<false, true, false, false> alignConfig;
    int score = globalAlignment(alignment, scoringScheme, alignConfig);

    return turnAlignmentIntoDescriptiveString(&alignment, 0, startTime, output, false, true);
}

// This function is used to conduct a short alignment for the sake of extending a GraphMap
// alignment.
char * endExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
                             int matchScore, int mismatchScore, int gapOpenScore,
                             int gapExtensionScore)
{
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();
    std::string output;

    TSequence sequenceH = read;
    TSequence sequenceV = ref;

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    // The only free gaps are at the end of ref (the reference sequence).
    AlignConfig<false, false, true, false> alignConfig;
    int score = globalAlignment(alignment, scoringScheme, alignConfig);

    return turnAlignmentIntoDescriptiveString(&alignment, 0, startTime, output, true, false);
}


// This function returns a list of the k-mers common to the two sequences.
std::vector<CommonKmer> getCommonKmers(std::string & readName, std::string & refName,
                                       double expectedSlope, int verbosity, std::string & output,
                                       KmerPositions * kmerPositions)
{
    std::vector<CommonKmer> commonKmers;
    double rotationAngle = CommonKmer::getRotationAngle(expectedSlope);

    KmerPosMap * readKmerPositions = kmerPositions->getKmerPositions(readName);
    KmerPosMap * refKmerPositions = kmerPositions->getKmerPositions(refName);

    KmerPosMap * smaller = readKmerPositions;
    KmerPosMap * larger = refKmerPositions;
    bool refKmersSmaller = false;
    if (smaller->size() > larger->size())
    {
        std::swap(smaller, larger);
        refKmersSmaller = true;
    }

    for (KmerPosMap::iterator i = smaller->begin(); i != smaller->end(); ++i)
    {
        std::string kmer = i->first;
        KmerPosMap::iterator j = larger->find(kmer);
        if (j != larger->end())
        {
            // If the code got here, then a common k-mer was found!
            std::vector<int> * readPositions = &(i->second);
            std::vector<int> * refPositions = &(j->second);
            if (refKmersSmaller)
                std::swap(readPositions, refPositions);

            for (int k = 0; k < readPositions->size(); ++k)
            {
                for (int l = 0; l < refPositions->size(); ++l)
                    commonKmers.push_back(CommonKmer(kmer, (*readPositions)[k], (*refPositions)[l], rotationAngle));
            }
        }
    }

    return commonKmers;



















































    // // Now that we've chose a good k-mer size, we can build the vector of CommonKmer objects.
    // std::unordered_set<std::string> * readKmers = readKmerSets->getKmerPositions(readName, readSeq, kSize);
    // std::unordered_set<std::string> * refKmers = refKmerSets->getKmerPositions(refName, refSeq, kSize);

    // // Loop through the reference sequence, and for all kmers also present in the read sequence,
    // // add to a map.
    // std::map<std::string, std::vector<int> > refPositions;
    // int kCount = refSeq.size() - kSize;
    // for (int i = 0; i < kCount; ++i)
    // {
    //     std::string kmer = refSeq.substr(i, kSize);
    //     if (readKmers->find(kmer) != readKmers->end()) // If kmer is in readKmers
    //     {
    //         if (refPositions.find(kmer) == refPositions.end()) // If kmer is not in refPositions
    //             refPositions[kmer] = std::vector<int>();
    //         refPositions[kmer].push_back(i);
    //     }
    // }
    
    // // Loop through the read sequence, and for all kmers also present in the reference, create a
    // // CommonKmer object for each reference position.
    // kCount = readSeq.size() - kSize;
    // for (int i = 0; i < kCount; ++i)
    // {
    //     std::string kmer = readSeq.substr(i, kSize);
    //     if (refKmers->find(kmer) != refKmers->end()) // If kmer is in refKmers
    //     {
    //         for (int j = 0; j < refPositions[kmer].size(); ++j)
    //             commonKmers.push_back(CommonKmer(kmer, i, refPositions[kmer][j], rotationAngle));
    //     }
    // }
    // return commonKmers;
}



// This function turns an alignment into a string which has the start/end positions, the CIGAR and
// the time it took to align.
char * turnAlignmentIntoDescriptiveString(Align<Dna5String, ArrayGaps> * alignment,
                                          int refOffset, long long startTime, std::string output,
                                          bool startImmediately, bool goToEnd)
{
    // Extract the alignment sequences into C++ strings, as the TRow type doesn't seem to have
    // constant time random access.
    std::ostringstream stream1;
    stream1 << row(*alignment, 0);
    std::string readAlignment =  stream1.str();
    std::ostringstream stream2;
    stream2 << row(*alignment, 1);
    std::string refAlignment =  stream2.str();

    int alignmentLength = std::max(readAlignment.size(), refAlignment.size());
    if (alignmentLength == 0)
        return cppStringToCString(output + ";Failed: alignment length zero");

    // Build a CIGAR string of the alignment.
    std::string cigarString;
    CigarType currentCigarType;
    int currentCigarLength = 0;
    int readStart = -1, refStart = -1;
    int readBases = 0, refBases = 0;

    bool alignmentStarted = false;
    bool readStarted = false, refStarted = false;

    if (startImmediately)
    {
        alignmentStarted = true;
        readStarted = true;
        refStarted = true;
        readStart = 0;
        refStart = 0;
    }

    for (int i = 0; i < alignmentLength; ++i)
    {
        char base1 = readAlignment[i];
        char base2 = refAlignment[i];

        // We consider the alignment to have started when we've encountered a base in both
        // sequences (though not necessarily at the same time).
        if (base1 != '-')
            readStarted = true;
        if (base2 != '-')
            refStarted = true;
        if (readStarted && refStarted && !alignmentStarted)
        {
            readStart = readBases;
            refStart = refBases;
            alignmentStarted = true;
        }

        CigarType cigarType = getCigarType(base1, base2, alignmentStarted);
        if (i == 0)
            currentCigarType = cigarType;
        if (cigarType == currentCigarType)
            ++currentCigarLength;
        else
        {
            cigarString.append(getCigarPart(currentCigarType, currentCigarLength));
            currentCigarType = cigarType;
            currentCigarLength = 1;
        }

        if (base1 != '-')
            ++readBases;
        if (base2 != '-')
            ++refBases;
    }

    int readEnd = readBases;
    int refEnd = refBases;
    if (currentCigarType == INSERTION && !goToEnd)
    {
        currentCigarType = CLIP;
        readEnd -= currentCigarLength;
    }
    else if (currentCigarType == DELETION && !goToEnd)
    {
        currentCigarType = NOTHING;
        refEnd -= currentCigarLength;
    }
    cigarString.append(getCigarPart(currentCigarType, currentCigarLength));

    long long endTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count(); 
    int millisecs = endTime - startTime;

    std::string finalString = output + ";" + 
                              cigarString + ";" +
                              std::to_string(readStart) + ";" + 
                              std::to_string(readEnd) + ";" + 
                              std::to_string(refStart + refOffset) + ";" + 
                              std::to_string(refEnd + refOffset) + ";" + 
                              std::to_string(millisecs);

    return cppStringToCString(finalString);
}


char * cppStringToCString(std::string cpp_string)
{
    char * c_string = (char*)malloc(sizeof(char) * (cpp_string.size() + 1));
    std::copy(cpp_string.begin(), cpp_string.end(), c_string);
    c_string[cpp_string.size()] = '\0';
    return c_string;
}

std::string vectorToString(std::vector<int> * v)
{
    std::stringstream ss;
    for(size_t i = 0; i < v->size(); ++i)
    {
        if (i != 0)
            ss << ",";
        ss << (*v)[i];
    }
    return ss.str();
}

CigarType getCigarType(char b1, char b2, bool alignmentStarted)
{
    if (b1 == '-')
    {
        if (alignmentStarted)
            return DELETION;
        else
            return NOTHING;
    }
    else if (b2 == '-')
    {
        if (alignmentStarted)
            return INSERTION;
        else
            return CLIP;
    }
    else
        return MATCH;
}

std::string getCigarPart(CigarType type, int length)
{
    std::string cigarPart = std::to_string(length);
    if (type == DELETION)
        cigarPart.append("D");
    else if (type == INSERTION)
        cigarPart.append("I");
    else if (type == CLIP)
        cigarPart.append("S");
    else if (type == MATCH)
        cigarPart.append("M");
    else //type == NOTHING
        return "";
    return cigarPart;
}


// This function gets the slope and intercept of the seed chain line as defined by the seeds at the
// given indices. Slope is defined as the change in V position over the change in H position.
// Intercept is defined as the V position at the H position of 0.
void getSeedChainSlopeAndIntercept(String<TSeed> * seedChain, int i1, int i2,
                                   double * slope, double * intercept)
{
    int hStart = beginPositionH((*seedChain)[i1]);
    int hEnd = endPositionH((*seedChain)[i2]);
    int vStart = beginPositionV((*seedChain)[i1]);
    int vEnd = endPositionV((*seedChain)[i2]);

    getSlopeAndIntercept(hStart, hEnd, vStart, vEnd, slope, intercept);
}


void getSlopeAndIntercept(int hStart, int hEnd, int vStart, int vEnd,
                          double * slope, double * intercept)
{
    *slope = -1.0;
    *intercept = -1.0;

    int hDiff = hEnd - hStart;
    int vDiff = vEnd - vStart;

    if (hDiff > 0)
    {
        *slope = double(vDiff) / hDiff;
        *intercept = vStart - (*slope * hStart);
    }
}


// Gets the median from an already-sorted vector.
double getMedian(std::vector<double> & v)
{
    size_t size = v.size();
    if (size % 2 == 0)
        return (v[size / 2 - 1] + v[size / 2]) / 2;
    else 
        return v[size / 2];
}


// Creates the CommonKmer object using the position in the two sequences.
// It then rotates the point relative to the origin by the given angle (in radians).
CommonKmer::CommonKmer(std::string sequence, int hPosition, int vPosition, double angle) :
    m_sequence(sequence),
    m_hPosition(hPosition),
    m_vPosition(vPosition),
    m_bandArea(0.0),
    m_score(1.0)
{
    double s = sin(angle);
    double c = cos(angle);
    m_rotatedHPosition = (m_hPosition * c) - (m_vPosition * s);
    m_rotatedVPosition = (m_hPosition * s) + (m_vPosition * c);
}

void printKmerSize(int kmerSize, int locationCount, std::string & output)
{
    output += "  " + std::to_string(locationCount) + " " + std::to_string(kmerSize) + "-mers in common\n";
}

// Given a point, a slope and rectangle bounds, this function returns the length of the line
// segment which goes through that point with that slope and is in those bounds.
double getLineLength(double x, double y, double slope, double xSize, double ySize)
{
    double xStart, yStart, xEnd, yEnd;

    double yIntercept = y - (slope * x);
    if (yIntercept >= 0.0)
    {
        xStart = 0.0;
        yStart = yIntercept;
    }
    else
    {
        xStart = -yIntercept / slope;
        yStart = 0.0;
    }

    double yAtXEnd = (slope * xSize) + yIntercept;
    if (yAtXEnd <= ySize)
    {
        xEnd = xSize;
        yEnd = yAtXEnd;
    }
    else
    {
        xEnd = (ySize - yIntercept) / slope;
        yEnd = ySize;
    }

    double xLength = xEnd - xStart;
    double yLength = yEnd - yStart;
    return sqrt((xLength * xLength) + (yLength * yLength));
}


// Adapted from:
// http://stackoverflow.com/questions/11449617/how-to-fit-the-2d-scatter-data-with-a-line-with-c
void linearRegression(std::vector<CommonKmer> & pts, double * slope, double * intercept)
{
    int n = pts.size();
    double sumH = 0.0, sumV = 0.0, sumHV = 0.0, sumHH = 0.0;
    for (int i = 0; i < n; ++i)
    {
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


void parseKmerLocationsFromString(std::string & str, std::vector<int> & v1, std::vector<int> & v2,
                                  int refOffset)
{
    std::stringstream ss(str);
    while(ss.good())
    {
        std::string readLocation;
        std::string refLocation;
        std::getline(ss, readLocation, ',');
        std::getline(ss, refLocation, ',');
        v1.push_back(std::stoi(readLocation));
        v2.push_back(std::stoi(refLocation) - refOffset);
    }
}

std::string getKmerTable(std::vector<CommonKmer> & commonKmers)
{
    std::string table = "\tSeq 1 pos\tSeq 2 pos\tRotated seq 1 pos\tRotated seq 2 pos\tBand area\tScore\n";
    for (int i = 0; i < commonKmers.size(); ++i)
    {
        table += "\t" + std::to_string(commonKmers[i].m_hPosition) +
                 "\t" + std::to_string(commonKmers[i].m_vPosition) + 
                 "\t" + std::to_string(commonKmers[i].m_rotatedHPosition) +
                 "\t" + std::to_string(commonKmers[i].m_rotatedVPosition) +
                 "\t" + std::to_string(commonKmers[i].m_bandArea) +
                 "\t" + std::to_string(commonKmers[i].m_score) + "\n";
    }
    return table;
}

std::string getSeedChainTable(String<TSeed> & seedChain)
{
    std::string table = "\tH start\tH end\tV start\tV end\n";
    int seedsInChain = length(seedChain);
    for (int i = 0; i < seedsInChain; ++i)
    {
        table += "\t" + std::to_string(beginPositionH(seedChain[i])) +
                 "\t" + std::to_string(endPositionH(seedChain[i])) +
                 "\t" + std::to_string(beginPositionV(seedChain[i])) +
                 "\t" + std::to_string(endPositionV(seedChain[i])) + "\n";
    }
    return table;
}

// This function adds a seed to a seed set. First it tries to merge the seed, and if that doesn't
// work it adds it as a separate seed.
void addSeedMerge(TSeedSet & seedSet, TSeed & seed)
{
    if (!addSeed(seedSet, seed, 1, Merge()))
        addSeed(seedSet, seed, Single());
}

void getMeanAndStDev(std::vector<double> & v, double & mean, double & stdDev)
{
    mean = 0.0;
    stdDev = 0.0;
    int count = v.size();
    if (count < 1)
        return;
    double devSum = 0.0;
    for (int i = 0; i < count; ++i)
        mean += v[i];
    mean /= count;
    for (int i = 0; i < count; ++i)
    {
        double dev = v[i] - mean;
        devSum += dev * dev;
    }
    stdDev = sqrt(devSum / v.size());
}

// This is the destructor for KmerPositions. It cleans up all the KmerPosMaps which were allocated
// on the heap.
KmerPositions::~KmerPositions()
{
    for (std::unordered_map<std::string, KmerPosMap *>::iterator i = m_kmerPositions.begin(); i != m_kmerPositions.end(); ++i)
        delete i->second;
}

// This function adds a sequence to the KmerPositions object. It creates a new KmerPosMap on the
// heap (will be deleted in destructor), fills it up and adds it to m_kmerPositions.
void KmerPositions::addPositions(char * nameC, char * sequenceC)
{
    std::string name(nameC);
    std::string sequence(sequenceC);

    KmerPosMap * posMap = new KmerPosMap();
    int kCount = sequence.size() - KMER_SIZE + 1;
    for (int i = 0; i < kCount; ++i)
    {
        std::string kmer = sequence.substr(i, KMER_SIZE);
        if (posMap->find(kmer) == posMap->end())
            (*posMap)[kmer] = std::vector<int>();
        (*posMap)[kmer].push_back(i);
    }

    m_kmerPositions[name] = posMap;
}

void addKmerPositions(KmerPositions * kmerPositions, char * nameC, char * sequenceC)
{
    kmerPositions->addPositions(nameC, sequenceC);
}

void KmerPositions::deletePositions(std::string & name)
{
    KmerPosMap * kmerPosMap = getKmerPositions(name);
    if (kmerPosMap != 0)
    {
        m_kmerPositions.erase(name);
        delete kmerPosMap;
    }
}

void deleteKmerPositions(KmerPositions * kmerPositions, char * nameC)
{
    std::string name(nameC);
    kmerPositions->deletePositions(name);
}

// This function retrieves a KmerPosMap from the object using the name as a key. If the name isn't
// in the map, it returns 0.
KmerPosMap * KmerPositions::getKmerPositions(std::string & name)
{
    if (m_kmerPositions.find(name) == m_kmerPositions.end())
        return 0;
    else
        return m_kmerPositions[name];
}

}
