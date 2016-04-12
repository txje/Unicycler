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
typedef Align<Dna5String, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;
typedef Seed<Simple> TSeed;
typedef SeedSet<TSeed> TSeedSet;

typedef std::tuple<std::string, int, int> Kmer;
typedef std::map<std::string, std::tuple<int, int>> KmerDict;
typedef std::tuple<int, int, int, int> CommonLocation;

std::vector<Kmer> getSeqKmers(std::string seq, int strLen, int kSize)
{
	std::vector<Kmer> kmers;
	int kCount = strLen - kSize;
	kmers.reserve(kCount);
	for (int i = 0; i < kCount; ++i)
	{
		Kmer kmer(seq.substr(i, kSize), i, i + kSize);
		kmers.push_back(kmer);
	}
	return kmers;
}

std::vector<CommonLocation> getCommonLocations(std::vector<Kmer> s1Kmers, std::vector<Kmer> s2Kmers)
{
	std::vector<CommonLocation> commonLocations;

	// Store all s1 kmers in a map of seq -> positions.
	KmerDict s1KmerPositions;
	for (int i = 0; i < s1Kmers.size(); ++i)
	{
		std::string s1KmerSeq = std::get<0>(s1Kmers[i]);
		std::tuple<int, int> s1Positions(std::get<1>(s1Kmers[i]), std::get<2>(s1Kmers[i]));
		s1KmerPositions[s1KmerSeq] = s1Positions;
	}

	// For all s2 kmers, see if they are in the s1 map. If so, they are common.
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

struct LocationCompare
{
    bool operator()(CommonLocation const &l1, CommonLocation const &l2)
    {
    	if (std::get<0>(l1) < std::get<0>(l2))
    		return true;
	    if (std::get<0>(l1) > std::get<0>(l2))
	    	return false;
	    return std::get<2>(l1) < std::get<2>(l2);
    }
};

std::vector<CommonLocation> mergeAdjacentLocations(std::vector<CommonLocation> locations)
{
	std::sort(std::begin(locations), std::end(locations), LocationCompare());

	std::vector<CommonLocation> mergedLocations;
	CommonLocation currentLocation = locations[0];
	int currentS1Start = std::get<0>(currentLocation);
	for (int i = 1; i < locations.size(); ++i)
	{
		CommonLocation location = locations[i];
		if (std::get<0>(location) == currentS1Start + 1)
		{
			std::get<1>(currentLocation) += 1;
			std::get<3>(currentLocation) += 1;
		}
		else
		{
			mergedLocations.push_back(currentLocation);
			currentLocation = location;
			currentS1Start = std::get<0>(location);
		}
	}
	mergedLocations.push_back(currentLocation);

	return mergedLocations;
}


void semiGlobalAlign(char * s1, char * s2, int s1Len, int s2Len)
{
	int kSize = 7;

	std::string s1Str(s1);
	std::string s2Str(s2);

	std::vector<Kmer> s1Kmers = getSeqKmers(s1Str, s1Len, kSize);
	std::vector<Kmer> s2Kmers = getSeqKmers(s2Str, s2Len, kSize);
	std::vector<CommonLocation> commonLocations = getCommonLocations(s1Kmers, s2Kmers);
	std::vector<CommonLocation> mergedLocations = mergeAdjacentLocations(commonLocations);

	for (int i = 0; i < mergedLocations.size(); ++i)
	{
		CommonLocation l = mergedLocations[i];
		std::cout << std::get<0>(l) << " " << std::get<1>(l) << " " << std::get<2>(l) << " " << std::get<3>(l) << std::endl;
	}




    long long ms_before = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    TSequence sequenceH = s1;
    TSequence sequenceV = s2;

    Score<int, Simple> scoringScheme(1, -1, -1);
    AlignConfig<true, true, true, true> alignConfig;

	TSeedSet seedSet;
	for (int i = 0; i < mergedLocations.size(); ++i)
	{
		CommonLocation l = mergedLocations[i];
		addSeed(seedSet, TSeed(std::get<0>(l), std::get<2>(l), std::get<1>(l), std::get<3>(l)), Single());
	}


	String<TSeed> seedChain;
	chainSeedsGlobally(seedChain, seedSet, SparseChaining());

	Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);

    int result = bandedChainAlignment(alignment, seedChain, scoringScheme, alignConfig, 20);

    std::cout << "Score: " << result << std::endl;
    std::cout << alignment << std::endl;






    // TSeed seed(1681, 1073, 1691, 1083);


    // extendSeed(seed, seq1, seq2, EXTEND_BOTH, scoringScheme, 3, GappedXDrop());

    // std::cout << "seed1: " << infix(seq1, beginPositionH(seed), endPositionH(seed)) << "\n"
    //           << "seed2: " << infix(seq2, beginPositionV(seed), endPositionV(seed)) << "\n";


    // TAlign align;
    // resize(rows(align), 2);
    // assignSource(row(align, 0), infix(seq1, beginPositionH(seed),
    //                                   endPositionH(seed)));
    // assignSource(row(align, 1), infix(seq2, beginPositionV(seed),
    //                                   endPositionV(seed)));

    // globalAlignment(align, scoringScheme, alignConfig);
    // std::cout << align << "\n";

    long long ms_after = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();
    long long elapsedTime = ms_after - ms_before;
    std::cout << "Milliseconds: " << elapsedTime << std::endl;
}


}
