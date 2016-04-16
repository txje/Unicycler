#include <stdio.h>
#include <chrono>
#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seeds.h>
#include <string>
#include <set>
#include <tuple>
#include <vector>
#include <map>
#include <algorithm>

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
typedef std::tuple<std::string, int, int> Kmer;
typedef std::map<std::string, std::tuple<int, int>> KmerDict;
typedef std::tuple<int, int, int, int> CommonLocation;

enum CigarType {MATCH, INSERTION, DELETION, CLIP, NOTHING};

// These are the functions that will be called by the Python script.
char * semiGlobalAlignmentAroundLine(char * s1, char * s2, int s1Len, int s2Len, double slope,
                                     double intercept, int bandSize, int debugOutput);
char * findAlignmentLine(char * s1, char * s2, int s1Len, int s2Len, double expectedSlope,
                         int debugOutput);
char * exhaustiveSemiGlobalAlignment(char * s1, char * s2, int s1Len, int s2Len);
void free_c_string(char * p);

// These functions are internal to this C++ code.
char * turnAlignmentIntoDescriptiveString(Align<Dna5String, ArrayGaps> * alignment,
                                          long long startTime);
std::vector<Kmer> getSeqKmers(std::string seq, int kSize);
std::vector<CommonLocation> getCommonLocations(std::string s1, std::string s2, int kSize);
CigarType getCigarType(char b1, char b2, bool alignmentStarted);
std::string getCigarPart(CigarType type, int length);
char * cppStringToCString(std::string cpp_string);
std::string vectorToString(std::vector<int> * v);
void fixSeedChainToLine(String<TSeed> * seedChain, double bandSize);
void getSeedChainSlopeAndIntercept(String<TSeed> * seedChain, int firstI, int lastI,
                                   double * slope, double * intercept);
void getSlopeAndIntercept(int hStart, int hEnd, int vStart, int vEnd,
                                   double * slope, double * intercept);
char * getSlopeAndInterceptString(double slope, double intercept);
double getMedian(std::vector<double> * v);




char * semiGlobalAlignmentAroundLine(char * s1, char * s2, int s1Len, int s2Len, double slope,
                                     double intercept, int bandSize, int debugOutput)
{
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    // Extreme slope values will not work, so check for them now.
    if (slope > 1.5)
        return strdup("Failed: slope too large");
    if (slope < 0.5)
        return strdup("Failed: slope too small");

    TSequence sequenceH = s1;
    TSequence sequenceV = s2;

    // Convert the slope and intercept into start/end coordinates for the two sequences.
    int hStart, hEnd, vStart, vEnd;
    if (intercept >= 0.0)
    {
        hStart = 0;
        vStart = std::round(intercept);
    }
    else
    {
        hStart = std::round(-intercept / slope);
        vStart = 0;
    }
    double vPosAtEndOfH = (slope * s1Len) + intercept;
    if (vPosAtEndOfH <= s2Len)
    {
        hEnd = s1Len;
        vEnd = std::round(vPosAtEndOfH);
    }
    else
    {
        hEnd = std::round((s2Len - intercept) / slope);
        vEnd = s2Len;
    }
    int hDiff = hEnd - hStart;
    int vDiff = vEnd - vStart;

    if (debugOutput > 0)
    {
        std::cout << std::endl;
        std::cout << "ALIGNMENT RANGE FROM SLOPE" << std::endl;
        std::cout << "--------------------------" << std::endl;
        std::cout << "hStart: " << hStart << ", " << "hEnd: " << hEnd << ", " << "vStart: " << vStart << ", " << "vEnd: " << vEnd << std::endl << std::endl;;
    }

    // Now that we have an estimate for the start/end of both sequences, we build a single seed
    // chain around that line.
    String<TSeed> seedChain;

    // If the two sequences are the same length, then a single seed will do.
    if (hDiff == vDiff)
        appendValue(seedChain, TSeed(hStart, vStart, hEnd, vEnd));

    // If one of the distances is longer, then we need multiple seeds to step across at a slope
    // not equal to 1.
    else
    {
        int additionalSteps = std::abs(vDiff - hDiff);
        bool additionalHSteps = hDiff > vDiff;

        int seedCount = additionalSteps + 1;
        double pieceSize = double(std::min(hDiff, vDiff)) / seedCount;        
        double targetPosition = 0;
        int distanceSoFar = 0;
        int hPos = hStart;
        int vPos = vStart;
        for (int i = 0; i < seedCount; ++i)
        {
            targetPosition += pieceSize;
            int thisPieceSize = std::round(targetPosition - distanceSoFar);
            distanceSoFar += thisPieceSize;
            appendValue(seedChain, TSeed(hPos, vPos, hPos + thisPieceSize, vPos + thisPieceSize));
            hPos += thisPieceSize;
            vPos += thisPieceSize;
            if (additionalHSteps)
                hPos += 1;
            else
                vPos += 1;
        }
    }

    if (debugOutput > 0)
    {
        int seedsInChain = length(seedChain);
        std::cout << std::endl;
        std::cout << "SEED POSITIONS" << std::endl;
        std::cout << "--------------" << std::endl;
        std::cout << "H start\tH end\tV start\tV end" << std::endl;
        for (int i = 0; i < seedsInChain; ++i)
        {
            std::cout << beginPositionH(seedChain[i]) << "\t" << endPositionH(seedChain[i]) << "\t";
            std::cout << beginPositionV(seedChain[i]) << "\t" << endPositionV(seedChain[i]) << std::endl;
        }
        std::cout << std::endl;
    }

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(1, -1, -1);
    AlignConfig<true, true, true, true> alignConfig;
    bandedChainAlignment(alignment, seedChain, scoringScheme, alignConfig, bandSize);

    return turnAlignmentIntoDescriptiveString(&alignment, startTime);
}

char * findAlignmentLine(char * s1, char * s2, int s1Len, int s2Len, double expectedSlope,
                         int debugOutput)
{
    std::string s1Str(s1);
    std::string s2Str(s2);

    // We will dynamically choose a k-mer size that gives a useful density of common locations.
    int targetKCount = double(s1Len) * double(s2Len) * 0.000002;
    if (targetKCount < 10)
        return strdup("Failed: sequence too short");
    int minimumKCount = targetKCount / 4;
    int kSize = 10;
    std::vector<CommonLocation> commonLocations = getCommonLocations(s1Str, s2Str, kSize);

    if (debugOutput > 0)
        std::cout << "K-MER: " << kSize << ", " << commonLocations.size() << " SITES IN COMMON" << std::endl;

    // If the starting k-mer gave too many locations, we increase it until we're in the correct
    // range.
    if (commonLocations.size() > targetKCount)
    {
        while (true)
        {
            ++kSize;
            commonLocations = getCommonLocations(s1Str, s2Str, kSize);
            if (debugOutput > 0)
                std::cout << "K-MER: " << kSize << ", " << commonLocations.size() << " SITES IN COMMON" << std::endl;

            // If we're reached the target range, that's good...
            if (commonLocations.size() <= targetKCount)
            {
                // But if we went under the minimum, then we need to back down one k-mer step, even
                // if it takes us over our target.
                if (commonLocations.size() < minimumKCount)
                {
                    --kSize;
                    commonLocations = getCommonLocations(s1Str, s2Str, kSize);
                    if (debugOutput > 0)
                        std::cout << "K-MER: " << kSize << ", " << commonLocations.size() << " SITES IN COMMON" << std::endl;
                }
                break;
            }
        }
    }

    // If the starting k-mer gave too few locations, we decrease it until we're in the correct
    // range.
    if (commonLocations.size() < minimumKCount)
    {
        while (true)
        {
            --kSize;

            if (kSize < 2)
                return strdup("Failed: kmer size too small");

            commonLocations = getCommonLocations(s1Str, s2Str, kSize);
            if (debugOutput > 0)
                std::cout << "K-MER: " << kSize << ", " << commonLocations.size() << " SITES IN COMMON" << std::endl;

            if (commonLocations.size() >= minimumKCount)
                break;
        }
    }

    if (debugOutput > 0)
        std::cout << std::endl;

    //We should now have a reasonably-sized set of common kmers.
    if (commonLocations.size() < 2)
        return strdup("Failed: too few common kmers");

    if (debugOutput > 1)
    {
        std::cout << std::endl;
        std::cout << "COMMON K-MER POSITIONS" << std::endl;
        std::cout << "----------------------" << std::endl;
        std::cout << "Seq 1 pos\tSeq 2 pos" << std::endl;
        for (int i = 0; i < commonLocations.size(); ++i)
        {
            CommonLocation l = commonLocations[i];
            std::cout << std::get<0>(l) << "\t" << std::get<2>(l) << std::endl;
        }
        std::cout << std::endl;
    }

    // Build a Seqan seed set using our common k-mers.
    TSeedSet seedSet;
    for (int i = 0; i < commonLocations.size(); ++i)
    {
        CommonLocation l = commonLocations[i];
        TSeed seed(std::get<0>(l), std::get<2>(l), std::get<1>(l), std::get<3>(l));
        if (!addSeed(seedSet, seed, 1, Merge()))
            addSeed(seedSet, seed, Single());
    }

    if (debugOutput > 1)
    {
        std::cout << std::endl;
        std::cout << "SEED POSITIONS BEFORE GLOBAL CHAINING" << std::endl;
        std::cout << "-------------------------------------" << std::endl;
        std::cout << "H start\tH end\tV start\tV end" << std::endl;
        for (Iterator<TSeedSet>::Type it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it)
        {
            std::cout << beginPositionH(*it) << "\t" << endPositionH(*it) << "\t";
            std::cout << beginPositionV(*it) << "\t" << endPositionV(*it) << std::endl;
        }
        std::cout << std::endl;
    }

    // We now get a Seqan global chain of the seeds.
    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());
    int seedsInChain = length(seedChain);
    if (seedsInChain == 0)
        return strdup("Failed: no global seed chain");

    if (debugOutput > 1)
    {
        std::cout << std::endl;
        std::cout << "SEED POSITIONS AFTER GLOBAL CHAINING" << std::endl;
        std::cout << "------------------------------------" << std::endl;
        std::cout << "H start\tH end\tV start\tV end" << std::endl;
        for (int i = 0; i < seedsInChain; ++i)
        {
            std::cout << beginPositionH(seedChain[i]) << "\t" << endPositionH(seedChain[i]) << "\t";
            std::cout << beginPositionV(seedChain[i]) << "\t" << endPositionV(seedChain[i]) << std::endl;
        }
        std::cout << std::endl;
    }

    // In the unlikely case that our kmers have all merged into a single seed, then we can use
    // that one seed's start/end to get our line.
    if (seedsInChain == 1)
    {
        int hStart = beginPositionH(seedChain[0]);
        int hEnd = endPositionH(seedChain[0]);
        int vStart = beginPositionV(seedChain[0]);
        int vEnd = endPositionV(seedChain[0]);
        double slope, intercept;
        getSlopeAndIntercept(hStart, hEnd, vStart, vEnd, &slope, &intercept);
        if (slope == -1.0 && intercept == -1.0)
            return strdup("Failed: single seed bad slope and intercept");
        return getSlopeAndInterceptString(slope, intercept);
    }

    // We now want to extract the most likely line which represents that global chain.
    // Sample pairs of seeds in the chain and determine their slope and intercept. The median of
    // these values will be our returned answer.
    int pairDistance = seedsInChain / 5; // TO DO: make this a parameter?
    std::vector<double> slopes;
    std::vector<double> intercepts;
    for (int i = 0; i < seedsInChain - pairDistance; ++i)
    {
        int hStart = beginPositionH(seedChain[i]);
        int vStart = beginPositionV(seedChain[i]);
        int hEnd = beginPositionH(seedChain[i + pairDistance]);
        int vEnd = beginPositionV(seedChain[i + pairDistance]);
        double slope, intercept;
        getSlopeAndIntercept(hStart, hEnd, vStart, vEnd, &slope, &intercept);
        if (slope == -1.0 && intercept == -1.0)
            return strdup("Failed: bad slope and intercept");
        slopes.push_back(slope);
        intercepts.push_back(intercept);
    }

    double medianSlope = getMedian(&slopes);
    double medianIntercept = getMedian(&intercepts);

    if (debugOutput > 1)
    {
        std::cout << std::endl;
        std::cout << "SLOPES AND INTERCEPTS (SORTED)" << std::endl;
        std::cout << "------------------------------" << std::endl;
        for (int i = 0; i < slopes.size(); ++i)
            std::cout << "Slope: " << slopes[i] << ", Intercept: " << intercepts[i] << std::endl;
        std::cout << std::endl;
    }
    if (debugOutput > 0)
        std::cout << "MEDIAN SLOPE: " << medianSlope << ", MEDIAN INTERCEPT: " << medianIntercept << std::endl << std::endl;

    return getSlopeAndInterceptString(medianSlope, medianIntercept);
}

// This function does the full semi-global alignment using the entirety of both sequences. It will
// be slow but will always find the ideal alignment.
char * exhaustiveSemiGlobalAlignment(char * s1, char * s2, int s1Len, int s2Len)
{
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    TSequence sequenceH = s1;
    TSequence sequenceV = s2;

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(1, -1, -1);
    AlignConfig<true, true, true, true> alignConfig;
    globalAlignment(alignment, scoringScheme, alignConfig);

    return turnAlignmentIntoDescriptiveString(&alignment, startTime);
}


char * turnAlignmentIntoDescriptiveString(Align<Dna5String, ArrayGaps> * alignment,
                                          long long startTime)
{
    // Extract the alignment sequences into C++ strings, as the TRow type doesn't seem to have
    // constant time random access.
    std::ostringstream stream1;
    stream1 << row(*alignment, 0);;
    std::string s1Alignment =  stream1.str();
    std::ostringstream stream2;
    stream2 << row(*alignment, 1);;
    std::string s2Alignment =  stream2.str();

    int alignmentLength = std::max(s1Alignment.size(), s2Alignment.size());
    if (alignmentLength == 0)
        return strdup("Failed: alignment length zero");

    // Build a CIGAR string of the alignment.
    std::string cigarString;
    CigarType currentCigarType;
    int currentCigarLength = 0;
    int matchCount = 0;
    int mismatchCount = 0;
    int insertionCount = 0;
    int deletionCount = 0;
    std::vector<int> s2MismatchPositions;
    std::vector<int> s2InsertionPositions;
    std::vector<int> s2DeletionPositions;
    int s1Start = -1, s2Start = -1;
    int s1Bases = 0, s2Bases = 0;
    bool alignmentStarted = false;
    for (int i = 0; i < alignmentLength; ++i)
    {
        char base1 = s1Alignment[i];
        char base2 = s2Alignment[i];

        if (base1 != '-' && base2 != '-' && !alignmentStarted)
        {
            s1Start = s1Bases;
            s2Start = s2Bases;
            alignmentStarted = true;
        }

        CigarType cigarType = getCigarType(base1, base2, alignmentStarted);
        if (i == 0)
            currentCigarType = cigarType;

        if (cigarType == MATCH)
        {
            if (base1 == base2)
                ++matchCount;
            else
            {
                ++mismatchCount;
                s2MismatchPositions.push_back(s2Bases);
            }
        }
        else if (cigarType == DELETION)
        {
            ++deletionCount;
            s2DeletionPositions.push_back(s2Bases);
        }
        else if (cigarType == INSERTION)
        {
            ++insertionCount;
            s2InsertionPositions.push_back(s2Bases);
        }

        if (cigarType == currentCigarType)
            ++currentCigarLength;
        else
        {
            cigarString.append(getCigarPart(currentCigarType, currentCigarLength));
            currentCigarType = cigarType;
            currentCigarLength = 1;
        }

        if (base1 != '-')
            ++s1Bases;
        if (base2 != '-')
            ++s2Bases;
    }

    int s1End = s1Bases;
    int s2End = s2Bases;
    if (currentCigarType == INSERTION)
    {
        currentCigarType = CLIP;
        insertionCount -= currentCigarLength;
        s2InsertionPositions.resize(insertionCount);
        s1End -= currentCigarLength;
    }
    else if (currentCigarType == DELETION)
    {
        currentCigarType = NOTHING;
        deletionCount -= currentCigarLength;
        s2DeletionPositions.resize(deletionCount);
        s2End -= currentCigarLength;
    }    
    cigarString.append(getCigarPart(currentCigarType, currentCigarLength));

    int editDistance = mismatchCount + insertionCount + deletionCount;
    int alignedLength = matchCount + mismatchCount + insertionCount + deletionCount;
    double percentIdentity = 100.0 * matchCount / alignedLength;

    long long endTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count(); 
    int milliseconds = endTime - startTime;

    std::string finalString = cigarString + "," +
                              std::to_string(s1Start) + "," + 
                              std::to_string(s1End) + "," + 
                              std::to_string(s2Start) + "," + 
                              std::to_string(s2End) + "," + 
                              std::to_string(alignedLength) + "," + 
                              std::to_string(matchCount) + "," + 
                              std::to_string(mismatchCount) + "," + 
                              vectorToString(&s2MismatchPositions) + "," + 
                              std::to_string(insertionCount) + "," + 
                              vectorToString(&s2InsertionPositions) + "," + 
                              std::to_string(deletionCount) + "," + 
                              vectorToString(&s2DeletionPositions) + "," + 
                              std::to_string(editDistance) + "," + 
                              std::to_string(percentIdentity) + "," + 
                              std::to_string(milliseconds);

    return cppStringToCString(finalString);
}

char * getSlopeAndInterceptString(double slope, double intercept)
{
    std::string slopeAndInterceptString = std::to_string(slope) + "," + std::to_string(intercept);
    return cppStringToCString(slopeAndInterceptString);

}

// Frees dynamically allocated memory for a c string. Called by Python after the string has been
// received.
void free_c_string(char * p)
{
    free(p);
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
            ss << ";";
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

// Returns a list of all Kmers in a sequence.
std::vector<Kmer> getSeqKmers(std::string seq, int kSize)
{
    std::vector<Kmer> kmers;
    int kCount = seq.size() - kSize;
    kmers.reserve(kCount);
    for (int i = 0; i < kCount; ++i)
    {
        Kmer kmer(seq.substr(i, kSize), i, i + kSize);
        kmers.push_back(kmer);
    }
    return kmers;
}

// Returns a list of all Kmers common to both lists.
std::vector<CommonLocation> getCommonLocations(std::string s1, std::string s2, int kSize)
{
    std::vector<Kmer> s1Kmers = getSeqKmers(s1, kSize);
    std::vector<Kmer> s2Kmers = getSeqKmers(s2, kSize);

    // Store all s1 kmers in a map of seq -> positions.
    KmerDict s1KmerPositions;
    for (int i = 0; i < s1Kmers.size(); ++i)
    {
        std::string s1KmerSeq = std::get<0>(s1Kmers[i]);
        std::tuple<int, int> s1Positions(std::get<1>(s1Kmers[i]), std::get<2>(s1Kmers[i]));
        s1KmerPositions[s1KmerSeq] = s1Positions;
    }

    // For all s2 kmers, see if they are in the s1 map. If so, they are common.
    std::vector<CommonLocation> commonLocations;
    for (int i = 0; i < s2Kmers.size(); ++i)
    {
        std::string s2KmerSeq = std::get<0>(s2Kmers[i]);
        if (s1KmerPositions.count(s2KmerSeq))
        {
            std::tuple<int, int> s1Position = s1KmerPositions[s2KmerSeq];
            int s1Start = std::get<0>(s1Position);
            int s1End = std::get<1>(s1Position);
            int s2Start = std::get<1>(s2Kmers[i]);
            int s2End = std::get<2>(s2Kmers[i]);
            CommonLocation commonLocation(s1Start, s1End, s2Start, s2End);
            commonLocations.push_back(commonLocation);
        }
    }

    return commonLocations;
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

double getMedian(std::vector<double> * v)
{
    size_t size = v->size();
    std::sort(v->begin(), v->end());
    if (size % 2 == 0)
        return ((*v)[size / 2 - 1] + (*v)[size / 2]) / 2;
    else 
        return (*v)[size / 2];
}

}
