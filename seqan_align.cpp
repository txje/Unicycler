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

enum CigarType {MATCH, INSERTION, DELETION};

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


void semiGlobalAlign(char * s1, char * s2, int s1Len, int s2Len, int kSize, int bandSize)
{
	long long time1 = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

	std::string s1Str(s1);
	std::string s2Str(s2);
	std::vector<Kmer> s1Kmers = getSeqKmers(s1Str, s1Len, kSize);
	std::vector<Kmer> s2Kmers = getSeqKmers(s2Str, s2Len, kSize);
	std::vector<CommonLocation> commonLocations = getCommonLocations(s1Kmers, s2Kmers);

    long long time2 = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

	TSeedSet seedSet;
	for (int i = 0; i < commonLocations.size(); ++i)
	{
		CommonLocation l = commonLocations[i];
		TSeed seed(std::get<0>(l), std::get<2>(l), std::get<1>(l), std::get<3>(l));
    	if (!addSeed(seedSet, seed, 1, Merge()))
        	addSeed(seedSet, seed, Single());
	}

    long long time3 = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    TSequence sequenceH = s1;
    TSequence sequenceV = s2;

	String<TSeed> seedChain;
	chainSeedsGlobally(seedChain, seedSet, SparseChaining());

    long long time4 = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

	Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(1, -1, -1);
    AlignConfig<true, true, true, true> alignConfig;
    int result = bandedChainAlignment(alignment, seedChain, scoringScheme, alignConfig, bandSize);

    long long time5 = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    // std::cout << alignment << std::endl;

	int alignmentLength = std::max(length(row(alignment, 0)), length(row(alignment, 1)));

	for (int i = 0; i < alignmentLength; ++i)
	{
		std::cout << row(alignment, 0)[i] << row(alignment, 1)[i] << std::endl;

	}

	// for (int i = 0; i < length(rows(alignment)); ++i)
 //    {
 //    	TRowIterator it = iter(row(alignment, i), 0);
 //        TRowIterator itEnd = iter(row(alignment, i), aliLength);
 //        int pos = 0;
 //        std::cout << "Row " << i << " contains gaps at positions: ";
 //        std::cout << std::endl;
 //        while (it != itEnd)
 //        {
 //            if (isGap(it))
 //                std::cout << pos << std::endl;
 //            ++it;
 //            ++pos;
 //        }
 //    }


    std::cout << "Score: " << result << std::endl << std::endl;

    std::cout << "Milliseconds to find common kmers: " << time2 - time1 << std::endl;
    std::cout << "Milliseconds to find add seeds:    " << time3 - time2 << std::endl;
    std::cout << "Milliseconds to find chain seeds:  " << time4 - time3 << std::endl;
    std::cout << "Milliseconds to perform alignment: " << time5 - time4 << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Total milliseconds:                " << time5 - time1 << std::endl;
}


}
