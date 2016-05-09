
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
    CommonKmer(std::string sequence, int hPosition, int vPosition, double angle);
    static double getRotationAngle(double slope) {return -atan(slope);}

    std::string m_sequence;
    int m_hPosition;
    int m_vPosition;
    double m_rotatedHPosition;
    double m_rotatedVPosition;
    double m_bandArea; // TO DO: DELETE LATER TO SAVE MEMORY
    double m_score; // Scaled kmer density
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
};


std::vector<CommonKmer> getCommonKmers(std::string & readName, std::string & refName,
                                       double expectedSlope, int verbosity, std::string & output,
                                       KmerPositions * kmerPositions);


KmerPositions * newKmerPositions();
void addKmerPositions(KmerPositions * kmerPositions, char * name, char * sequence);
void deleteKmerPositions(KmerPositions * kmerPositions, char * name);
void deleteAllKmerPositions(KmerPositions * kmerPositions);


#endif // KMERS_H

