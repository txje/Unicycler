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
#include <random>

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
                                     double intercept, int bandSize);
char * findAlignmentLine(char * s1, char * s2, int s1Len, int s2Len);
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
bool areSeedsInBand(String<TSeed> * seedChain, int i1, int i2, double slope, double lowIntercept,
                    double highIntercept);
bool isSeedInBand(TSeed * seed, double slope, double lowIntercept, double highIntercept);
char * getSlopeAndInterceptString(double slope, double intercept);
double getMedian(std::vector<double> * v);




char * semiGlobalAlignmentAroundLine(char * s1, char * s2, int s1Len, int s2Len, double slope,
                                     double intercept, int bandSize)
{
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    TSequence sequenceH = s1;
    TSequence sequenceV = s2;

    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE
    // TO DO: BUILD SEED CHAIN AROUND LINE

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(1, -1, -1);
    AlignConfig<true, true, true, true> alignConfig;
    bandedChainAlignment(alignment, seedChain, scoringScheme, alignConfig, bandSize);

    return turnAlignmentIntoDescriptiveString(&alignment, startTime);
}

char * findAlignmentLine(char * s1, char * s2, int s1Len, int s2Len)
{
    std::string s1Str(s1);
    std::string s2Str(s2);

    int targetKCount = std::min(s1Len, s2Len);
    if (targetKCount < 10)
        return strdup("Failed: sequence too short");
    int minimumKCount = targetKCount / 4;

    // We will dynamically choose a k-mer size that gives a useful number of common locations.
    int kSize = 10;
    std::vector<CommonLocation> commonLocations = getCommonLocations(s1Kmers, s2Kmers, kSize);

    // If the starting k-mer gave too many locations, we increase it until we're in the correct
    // range.
    if (commonLocations.size() > targetKCount)
    {
        while (true)
        {
            ++kSize;
            commonLocations = getCommonLocations(s1Kmers, s2Kmers, kSize);

            // If we're reached the target range, that's good...
            if (commonLocations.size() <= targetKCount)
            {
                // But if we went under the minimum, then we need to back down one k-mer step, even
                // if it takes us over our target.
                if (commonLocations.size() < minimumKCount)
                {
                    --kSize;
                    commonLocations = getCommonLocations(s1Kmers, s2Kmers, kSize);
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

            commonLocations = getCommonLocations(s1Kmers, s2Kmers, kSize);
            if (commonLocations.size() >= minimumKCount)
                break;
        }
    }

    //We should now have a reasonably-sized set of common kmers.
    if (commonLocations.size() < 2)
        return strdup("Failed: too few common kmers");

    // // Debugging code: output common kmer positions
    // std::cout << std::endl << "Common kmer positions" << std::endl;
    // for (int i = 0; i < commonLocations.size(); ++i)
    // {
    //     CommonLocation l = commonLocations[i];
    //     std::cout << std::get<0>(l) << "\t" << std::get<2>(l) << std::endl;
    // }

    // Build a Seqan seed set using our common k-mers.
    TSeedSet seedSet;
    for (int i = 0; i < commonLocations.size(); ++i)
    {
        CommonLocation l = commonLocations[i];
        TSeed seed(std::get<0>(l), std::get<2>(l), std::get<1>(l), std::get<3>(l));
        if (!addSeed(seedSet, seed, 1, Merge()))
            addSeed(seedSet, seed, Single());
    }

    // // Debugging code: output seed positions before chaining
    // std::cout << std::endl << "Seed positions before global chaining" << std::endl;
    // for (Iterator<TSeedSet>::Type it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it)
    // {
    //     std::cout << beginPositionH(*it) << "\t" << beginPositionV(*it) << std::endl;
    //     std::cout << endPositionH(*it) << "\t" << endPositionV(*it) << std::endl;
    // }

    // We now get a Seqan global chain of the seeds.
    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());
    int seedsInChain = length(seedChain);
    if (seedsInChain == 0)
        return strdup("Failed: no global seed chain");

    // // Debugging code: output seed positions after chaining
    // std::cout << std::endl << "Seed positions after global chaining" << std::endl;
    // for (int i = 0; i < seedsInChain; ++i)
    // {
    //     std::cout << beginPositionH(seedChain[i]) << "\t" << beginPositionV(seedChain[i]) << std::endl; //TEMP
    //     std::cout << endPositionH(seedChain[i]) << "\t" << endPositionV(seedChain[i]) << std::endl; //TEMP
    // }

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
            return strdup("Failed: bad slope and intercept");
        return getSlopeAndInterceptString(slope, intercept);
    }

    // We now want to extract the most likely line which represents that global chain.
    // Randomly sample pairs of seeds in the chain and determine their slope and intercept. The
    // median of these values will be our returned answer.
    int randomSampleCount = 100;
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, seedsInChain - 1);
    std::vector<double> slopes;
    std::vector<double> intercepts;
    for (int i = 0; i < randomSampleCount; ++i)
    {
        int i1, i2;
        i1 = random_integer = uni(rng);
        do
        {
            i2 = random_integer = uni(rng);
        } while (i1 == i2);
        if (i1 > i2)
            std::swap(i1, i2);

        int hStart = beginPositionH(seedChain[i1]);
        int hEnd = beginPositionH(seedChain[i1]);
        int vStart = beginPositionV(seedChain[i2]);
        int vEnd = beginPositionV(seedChain[i2]);
        double slope, intercept;
        getSlopeAndIntercept(hStart, hEnd, vStart, vEnd, &slope, &intercept);
        if (slope == -1.0 && intercept == -1.0)
            return strdup("Failed: bad slope and intercept");
        slopes.push_back(slope);
        intercepts.push_back(intercept);
    }

    double medianSlope = getMedian(&slopes);
    double medianIntercept = getMedian(&intercepts);
    return getSlopeAndInterceptString(slope, intercept);
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

// This function checks to see if all of the seeds within the given index range in the chain fall
// between two lines.
bool areSeedsInBand(String<TSeed> * seedChain, int i1, int i2, double slope, double lowIntercept,
                    double highIntercept)
{
    int seedsInChain = length(*seedChain);
    int h, v, lowerV, upperV;
    for (int i = i1; i <= i2; ++i)
    {
        if (!isSeedInBand(&((*seedChain)[i]), slope, lowIntercept, highIntercept))
            return false;
    }
    return true;
}

bool isSeedInBand(TSeed * seed, double slope, double lowIntercept, double highIntercept)
{
    int h = beginPositionH(*seed);
    int v = beginPositionV(*seed);
    if (v < slope * h + lowIntercept)
        return false;
    if (v > slope * h + highIntercept)
        return false;
    return true;
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
