
#ifndef KMERS_H
#define KMERS_H

#include <string>
#include <cmath>
#include <unordered_map>
#include <vector>
#include "settings.h"

typedef std::unordered_map<std::string, std::vector<int> > KmerPosMap;

// CommonKmer is a class to hold the positions of common k-mers between a read and a reference.
// It holds the position in both original and rotated coordinates, and will ultimately hold the
// point's score as well.
class CommonKmer {
public:
    CommonKmer(int hPosition, int vPosition, float angle);
    static float getRotationAngle(float slope) {return -atan(slope);}

    int m_hPosition;
    int m_vPosition;
    float m_rotatedHPosition;
    float m_rotatedVPosition;
    float m_score; // Scaled kmer density
};


// KmerPositions is a class that holds maps of k-mer positions for named sequences. It exists so we
// don't have to repeatedly find the same k-mer sets over and over when finding alignment lines.
class KmerPositions {
public:
    KmerPositions() {}
    ~KmerPositions();
    void addPositions(char * nameC, char * sequenceC);
    void deletePositions(char * nameC);
    KmerPosMap * getKmerPositions(std::string & name);

private:
    std::unordered_map<std::string, KmerPosMap *> m_kmerPositions;
    std::unordered_map<std::string, std::string> m_sequences;
};


// This class holds a set of common k-mers for a particular read-ref combination.
class CommonKmerSet {
public:
    CommonKmerSet(std::string & readName, std::string & refName,
                  int readLength, int refLength, int bandSize,
                  float expectedSlope, KmerPositions * kmerPositions);

    std::vector<CommonKmer> m_commonKmers;
    float m_maxScore;
};


float getLineLength(float x, float y, float slope, float xSize, float ySize);


// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {

    KmerPositions * newKmerPositions();

    void addKmerPositions(KmerPositions * kmerPositions, char * name, char * sequence);

    void deleteKmerPositions(KmerPositions * kmerPositions, char * name);

    void deleteAllKmerPositions(KmerPositions * kmerPositions);

}

#endif // KMERS_H

