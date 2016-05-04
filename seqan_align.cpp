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
    double m_bandLength; // Length of band at expected slope
    int m_bandCount; // Number of neighbour points in band
    double m_score; // Score calculated from bandCount, bandLength and overall kmer density
};


// KmerSets is a class that holds sets of k-mers for named sequences. It exists so we don't have to
// repeatedly find the same k-mer sets over and over. Whenever it makes a k-mer set, it stores it
// for later. On destruction it deletes all of its k-mer sets.
class KmerSets
{
public:
    KmerSets() {}
     ~KmerSets();
    std::unordered_set<std::string> * getKmerSet(std::string & name, std::string & sequence, int kSize);

private:
    std::unordered_map<std::string, std::unordered_set<std::string> *> m_kmerSets;
    std::mutex m_mutex;
};


// These are the functions that will be called by the Python script.
char * bandedSemiGlobalAlignment(char * read, int readLen,
                                 char * ref, int refLen, int refOffset,
                                 double slope, double intercept, int kSize, int bandSize,
                                 int verbosity,
                                 int matchScore, int mismatchScore,
                                 int gapOpenScore, int gapExtensionScore,
                                 char * kmerLocations);
char * findAlignmentLines(char * readSeqC, char * readNameC, char * refSeqC, char * refNameC,
                          double expectedSlope, int verbosity,
                          KmerSets * readKmerSets, KmerSets * refKmerSets);
char * startExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
                               int matchScore, int mismatchScore, int gapOpenScore,
                               int gapExtensionScore);
char * endExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
                             int matchScore, int mismatchScore, int gapOpenScore,
                             int gapExtensionScore);
void free_c_string(char * p) {free(p);}
KmerSets * newKmerSets() {return new KmerSets();}
void deleteKmerSets(KmerSets * kmerSets) {delete kmerSets;}


// These functions are internal to this C++ code.
void addBridgingSeed(TSeedSet & seedSet, int hStart, int vStart, int hEnd, int vEnd);
std::vector<CommonKmer> getCommonKmers(std::string & readSeq, std::string & readName,
                                       std::string & refSeq, std::string & refName,
                                       double expectedSlope, int verbosity, std::string & output,
                                       KmerSets * readKmerSets, KmerSets * refKmerSets);
long long getCommonKmerCount(std::string & readSeq, std::string & readName,
                             std::string & refSeq, std::string & refName,
                             int kSize, KmerSets * readKmerSets, KmerSets * refKmerSets);
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
std::unordered_set<std::string> * makeKmerSet(std::string & sequence, int kSize);
double getScoreThreshold(std::vector<CommonKmer> & commonKmers, int verbosity, std::string & output);
void getMeanAndStDev(std::vector<double> & v, double & mean, double & stdev);



// This function, given a line, will search for semi-global alignments around that line. The
// bandSize parameter specifies how far of an area around the line is searched. It is generally
// much faster than exhaustiveSemiGlobalAlignment, though it may not find the optimal alignment.
// A lower bandSize is faster with a larger chance of missing the optimal alignment.
char * bandedSemiGlobalAlignment(char * read, int readLen,
                                 char * ref, int refLen, int refOffset,
                                 double slope, double intercept, int kSize, int bandSize,
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
        TSeed seed(readKmerLocations[i], refKmerLocations[i], kSize);
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

// This function searches for lines in the 2D ref-read space that represent likely semi-global
// alignments.
char * findAlignmentLines(char * readSeqC, char * readNameC, char * refSeqC, char * refNameC,
                          double expectedSlope, int verbosity,
                          KmerSets * readKmerSets, KmerSets * refKmerSets)
{
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    std::string readSeq(readSeqC);
    std::string readName(readNameC);
    std::string refSeq(refSeqC);
    std::string refName(refNameC);
    std::string output;

    std::vector<CommonKmer> commonKmers = getCommonKmers(readSeq, readName, refSeq, refName,
                                                         expectedSlope, verbosity, output,
                                                         readKmerSets, refKmerSets);
    if (commonKmers.size() < 2)
    {
        long long endTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();
        return cppStringToCString(output + ";" + std::to_string(endTime - startTime) + ";Failed: too few common kmers");
    }

    int kSize = commonKmers[0].m_sequence.length();

    double commonKmerDensity = double(commonKmers.size()) / (double(readSeq.length()) * double(refSeq.length()));
    if (verbosity > 4)
        output += "  Common k-mer density: " + std::to_string(commonKmerDensity) + "\n";

    // Sort by rotated vertical position so lines should be roughly horizontal.
    std::sort(commonKmers.begin(), commonKmers.end(), [](const CommonKmer & a, const CommonKmer & b) {
        return a.m_rotatedVPosition < b.m_rotatedVPosition;   
    });

    // Score each point based on the number of other points in its band. The score is scaled by
    // the length of the band so short bands aren't penalised.
    double bandSize = 200.0; // TO DO: MAKE THIS A PARAMETER?
    int windowStart = 0;
    int windowEnd = 0;
    double maxScore = 0.0;
    for (int i = 0; i < commonKmers.size(); ++i)
    {
        double y = commonKmers[i].m_rotatedVPosition;
        while (y - commonKmers[windowStart].m_rotatedVPosition > bandSize)
            ++windowStart;
        while (windowEnd < commonKmers.size() &&
               commonKmers[windowEnd].m_rotatedVPosition - y < bandSize)
            ++windowEnd;
        double bandLength = getLineLength(commonKmers[i].m_hPosition, commonKmers[i].m_vPosition,
                                          expectedSlope, readSeq.length(), refSeq.length());
        int bandCount = windowEnd - windowStart;
        double bandKmerDensity = bandCount / (bandLength * 2.0 * bandSize);
        double score = bandKmerDensity / commonKmerDensity;
        if (score > maxScore)
            maxScore = score;

        commonKmers[i].m_bandLength = bandLength;
        commonKmers[i].m_bandCount = bandCount;
        commonKmers[i].m_score = score;
    }

    if (verbosity > 4)
    {
        output += "  Common k-mer positions:\n";
        output += getKmerTable(commonKmers);
        output += "  Max score: " + std::to_string(maxScore) + "\n";
    }

    // Now group all of the line points. A line group begins when the score exceeds a threshold and
    // it ends when the score drops below the threshold.
    std::vector<std::vector<CommonKmer> > lineGroups;
    double scoreThreshold = getScoreThreshold(commonKmers, verbosity, output);
    bool lineInProgress = false;
    for (int i = 0; i < commonKmers.size(); ++i)
    {
        if (commonKmers[i].m_score > scoreThreshold)
        {
            if (!lineInProgress)
            {
                lineGroups.push_back(std::vector<CommonKmer>());
                lineInProgress = true;
            }
            lineGroups.back().push_back(commonKmers[i]);
        }
        else // score is below threshold
            lineInProgress = false;
    }

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

    // Perform a simple least-squares linear regression on each line group.
    // If the slope is too far from 1.0, we throw the line out.
    double minSlope = 0.5; // TO DO: MAKE THIS A PARAMETER?
    double maxSlope = 1.5; // TO DO: MAKE THIS A PARAMETER?
    std::vector<double> slopes;
    std::vector<double> intercepts;
    std::vector<std::vector<CommonKmer> > filteredLineGroups;
    for (int i = 0; i < lineGroups.size(); ++i)
    {
        double slope, intercept;
        linearRegression(lineGroups[i], &slope, &intercept);
        if (slope >= minSlope && slope <= maxSlope)
        {
            filteredLineGroups.push_back(lineGroups[i]);
            slopes.push_back(slope);
            intercepts.push_back(intercept);
        }
    }
    lineGroups = filteredLineGroups;

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
        linesString += std::to_string(slopes[i]) + "," + std::to_string(intercepts[i]) + "," + std::to_string(kSize);

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
std::vector<CommonKmer> getCommonKmers(std::string & readSeq, std::string & readName,
                                       std::string & refSeq, std::string & refName,
                                       double expectedSlope, int verbosity, std::string & output,
                                       KmerSets * readKmerSets, KmerSets * refKmerSets)
{
    std::vector<CommonKmer> commonKmers;
    double rotationAngle = CommonKmer::getRotationAngle(expectedSlope);

    // We will dynamically choose a k-mer size that gives the maximum number of common k-mers.
    // Sizes too large will give fewer k-mers because it is less likely that sequences will share a
    // large k-mer. Sizes too small will give fewer k-mers because there are fewer possible k-mers.
    // Though we scale the value a little by multiplying by the k-mer size. This slightly pushs the
    // function towards larger k-mers.
    int startingKSize = 8;
    int kSize = startingKSize;
    long long commonKmerCount = getCommonKmerCount(readSeq, readName, refSeq, refName,
                                                   kSize, readKmerSets, refKmerSets);
    long long bestScore = commonKmerCount * kSize;
    int bestK = kSize;
    if (verbosity > 3) printKmerSize(kSize, commonKmerCount, output);

    // Try larger k sizes, to see if that helps.
    while (true)
    {
        ++kSize;
        commonKmerCount = getCommonKmerCount(readSeq, readName, refSeq, refName,
                                             kSize, readKmerSets, refKmerSets);
        long long score = commonKmerCount * kSize;
        if (verbosity > 3) printKmerSize(kSize, commonKmerCount, output);
        if (score > bestScore)
        {
            bestK = kSize;
            bestScore = score;
        }
        else
            break;
    }

    // If larger k sizes didn't help, try smaller k sizes.
    if (bestK == startingKSize)
    {
        kSize = startingKSize;
        while (true)
        {
            --kSize;
            commonKmerCount = getCommonKmerCount(readSeq, readName, refSeq, refName,
                                                 kSize, readKmerSets, refKmerSets);
            long long score = commonKmerCount * kSize;
            if (verbosity > 3) printKmerSize(kSize, commonKmerCount, output);
            if (score > bestScore)
            {
                bestK = kSize;
                bestScore = score;
            }
            else
                break;
        }
    }

    if (verbosity > 3)
        output += "  Best k size: " + std::to_string(bestK) + "\n";

    // Now that we've chose a good k-mer size, we can build the vector of CommonKmer objects.
    kSize = bestK;
    std::unordered_set<std::string> * readKmers = readKmerSets->getKmerSet(readName, readSeq, kSize);
    std::unordered_set<std::string> * refKmers = refKmerSets->getKmerSet(refName, refSeq, kSize);

    // Loop through the reference sequence, and for all kmers also present in the read sequence,
    // add to a map.
    std::map<std::string, std::vector<int> > refPositions;
    int kCount = refSeq.size() - kSize;
    for (int i = 0; i < kCount; ++i)
    {
        std::string kmer = refSeq.substr(i, kSize);
        if (readKmers->find(kmer) != readKmers->end()) // If kmer is in readKmers
        {
            if (refPositions.find(kmer) == refPositions.end()) // If kmer is not in refPositions
                refPositions[kmer] = std::vector<int>();
            refPositions[kmer].push_back(i);
        }
    }
    
    // Loop through the read sequence, and for all kmers also present in the reference, create a
    // CommonKmer object for each reference position.
    kCount = readSeq.size() - kSize;
    for (int i = 0; i < kCount; ++i)
    {
        std::string kmer = readSeq.substr(i, kSize);
        if (refKmers->find(kmer) != refKmers->end()) // If kmer is in refKmers
        {
            for (int j = 0; j < refPositions[kmer].size(); ++j)
                commonKmers.push_back(CommonKmer(kmer, i, refPositions[kmer][j], rotationAngle));
        }
    }
    return commonKmers;
}

// This function returns the number of k-mers the two sequences have in common.
long long getCommonKmerCount(std::string & readSeq, std::string & readName,
                             std::string & refSeq, std::string & refName,
                             int kSize, KmerSets * readKmerSets, KmerSets * refKmerSets)
{
    // Get the k-mer sets for both read and reference. If they have been accessed before, then this
    // should be very fast as they'll already exist. But if they haven't, the set(s) will be built
    // and this will take a bit longer.
    std::unordered_set<std::string> * readKmers = readKmerSets->getKmerSet(readName, readSeq, kSize);
    std::unordered_set<std::string> * refKmers = refKmerSets->getKmerSet(refName, refSeq, kSize);

    // Count the k-mers in both sets.
    std::unordered_set<std::string> * smallerSet = readKmers;
    std::unordered_set<std::string> * largerSet = refKmers;
    if (readKmers->size() > refKmers->size())
        std::swap(smallerSet, largerSet);
    int intersectionSize = 0;
    for (std::unordered_set<std::string>::iterator i = smallerSet->begin(); i != smallerSet->end(); ++i)
    {
        if (largerSet->find(*i) != largerSet->end())
            ++intersectionSize;
    }
    return intersectionSize;
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
    m_score(0.0)
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
    std::string table = "\tSeq 1 pos\tSeq 2 pos\tRotated seq 1 pos\tRotated seq 2 pos\tBand length\tBand count\tScore\n";
    for (int i = 0; i < commonKmers.size(); ++i)
    {
        table += "\t" + std::to_string(commonKmers[i].m_hPosition) +
                 "\t" + std::to_string(commonKmers[i].m_vPosition) + 
                 "\t" + std::to_string(commonKmers[i].m_rotatedHPosition) +
                 "\t" + std::to_string(commonKmers[i].m_rotatedVPosition) +
                 "\t" + std::to_string(commonKmers[i].m_bandLength) +
                 "\t" + std::to_string(commonKmers[i].m_bandCount) +
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


// This is the destructor for KmerSets. It cleans up all stuff allocated on the heap.
KmerSets::~KmerSets()
{
    for (std::unordered_map<std::string, std::unordered_set<std::string> *>::iterator i = m_kmerSets.begin(); i != m_kmerSets.end(); ++i)
        delete i->second;
}

// This function returns a k-mer set for a given sequence and k size. If the k-mer set already
// exists, it will just return it immediately. If not, it will build it using the sequence.
// The function is locked by a mutex so multiple threads don't try to build the same k-mer set at
// the same time.
std::unordered_set<std::string> * KmerSets::getKmerSet(std::string & name, std::string & sequence, int kSize)
{
    m_mutex.lock();
    std::string mapKey = name + "_" + std::to_string(kSize);
    if (m_kmerSets.find(mapKey) == m_kmerSets.end())
        m_kmerSets[mapKey] = makeKmerSet(sequence, kSize);
    m_mutex.unlock();
    return m_kmerSets[mapKey];
}


// This function produces a k-mer set for the given sequence. It allocates it on the heap and so
// will need to be deleted later.
std::unordered_set<std::string> * makeKmerSet(std::string & sequence, int kSize)
{
    std::unordered_set<std::string> * kmerSet = new std::unordered_set<std::string>();
    int kCount = sequence.size() - kSize;
    for (int i = 0; i < kCount; ++i)
        kmerSet->insert(sequence.substr(i, kSize));
    return kmerSet;
}


// This function determines a good score threshold where points exceeding the threshold are
// considered part of a line and points below are not.
double getScoreThreshold(std::vector<CommonKmer> & commonKmers, int verbosity, std::string & output)
{
    int count = commonKmers.size();

    // First we set an initial threshold by taking the mean of the 1% and 99% percentiles.
    std::vector<double> scores;
    scores.reserve(count);
    for (int i = 0; i < count; ++i)
        scores.push_back(commonKmers[i].m_score);
    std::sort(scores.begin(), scores.end());
    int firstPercentileIndex = int((count / 100.0) + 0.5) - 1;
    int ninetyNinthPercentileIndex = int((99.0 * count / 100.0) + 0.5) - 1;
    double firstPercentileScore = scores[firstPercentileIndex];
    double ninetyNinthPercentileScore = scores[ninetyNinthPercentileIndex];
    if (verbosity > 4)
    {
        output += "  1st percentile score: " + std::to_string(firstPercentileScore) + "\n";
        output += "  99th percentile score: " + std::to_string(ninetyNinthPercentileScore) + "\n";
    }

    // We expect the scores to be distributed in three different ways:
    //  1) Unimodal distribution near 1.0 (no line)
    //  2) Bimodal distribution: the first near 1.0 (points not in a line) and the second at a
    //     considerably higher value (points in a line).
    //  3) Unimodal distribution well over 1.0 (all point are in a line).

    // We first distinguish between the unimodal possibilities and the bimodal possibility.
    // This is done by getting the mean and stdev for the scores below and above the threshold.
    // If the distributions appear to overlap, we decide the scores are unimodal.
    if (ninetyNinthPercentileScore > firstPercentileScore)
    {
        double threshold = (firstPercentileScore + ninetyNinthPercentileScore) / 2.0;

        std::vector<double> scoresBelow;
        std::vector<double> scoresAbove;
        for (int i = 0; i < count; ++i)
        {
            if (scores[i] < threshold)
                scoresBelow.push_back(scores[i]);
            else
                scoresAbove.push_back(scores[i]);
        }
        double lowGroupMean, lowGroupStdDev;
        getMeanAndStDev(scoresBelow, lowGroupMean, lowGroupStdDev);
        double highGroupMean, highGroupStdDev;
        getMeanAndStDev(scoresAbove, highGroupMean, highGroupStdDev);

        // The zValue is how many standard deviations away from the means we must go before the
        // upper and lower distributions meet. A low zValue means the two groups are close and
        // likely to be two halves of a unimodal distribution. A high zValue means that the
        // two groups are separate and probably from a bimodal distribution.
        double zValue = (highGroupMean - lowGroupMean) / (highGroupStdDev + lowGroupStdDev);
        double zValueCutoff = 2.0; //TO DO: MAKE A PARAMETER?
        if (verbosity > 4)
            output += "  threshold of " + std::to_string(threshold) + " has zValue of " + std::to_string(zValue) + "\n";
        if (zValue > zValueCutoff)
        {
            if (verbosity > 4)
                output += "  Exceeded threshold, using bimodal distribution\n";
            return threshold;
        }
    }

    // If the code got here, then we seem to have a unimodal distribution. This means that either
    // none of the points are in a line (more likely) or all of the points are in a line (less
    // likely). We will distinguish between the two by looking at the variance of rotated point
    // positions. If it's very small, we assume the points all form a line. If it's large, we
    // assume no line.
    double stdDevThreshold = 100.0; // TO DO: come up with a better empirical value for this.
    std::vector<double> rotatedVPositions;
    rotatedVPositions.reserve(count);
    for (int i = 0; i < count; ++i)
        rotatedVPositions.push_back(commonKmers[i].m_rotatedVPosition);
    double mean, stdDev;
    getMeanAndStDev(rotatedVPositions, mean, stdDev);

    // A low standard deviation means all points are in a line. Set the threshold to 0 so all
    // points are included.
    if (stdDev < stdDevThreshold)
    {
        if (verbosity > 4)
            output += "  Unimodal distribution, all points in line\n";
        return 0.0;
    }

    //A high standard deviation means no points are in a line. Set the threshold to a high value
    // so all points are excluded.
    else
    {
        if (verbosity > 4)
            output += "  Unimodal distribution, no line detected\n";
        return std::numeric_limits<double>::max();
    }
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

}
