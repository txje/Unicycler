
#ifndef KMERS_H
#define KMERS_H

#include <string>
#include <cmath>
#include <unordered_map>
#include <vector>
#include "settings.h"
#include <mutex>

typedef std::unordered_map<std::string, std::vector<int> > KmerPosMap;


class CommonKmer {
public:
    CommonKmer(int hPosition, int vPosition);
    int m_hPosition;
    int m_vPosition;
};


// KmerPositions is a class that holds maps of k-mer positions for named sequences. It exists so we
// don't have to repeatedly find the same k-mer sets over and over.
class KmerPositions {
public:
    KmerPositions() {}
    ~KmerPositions();
    void addPositions(std::string & name, std::string & sequence, int kSize);
    KmerPosMap * getKmerPositions(std::string & name);
    std::string * getSequence(std::string & name);
    std::vector<std::string> getAllNames();
    int getLength(std::string & name);

private:
    std::unordered_map<std::string, KmerPosMap *> m_kmerPositions;
    std::unordered_map<std::string, std::string> m_sequences;
    std::mutex m_mutex;
};

KmerPositions * newKmerPositions();

void addKmerPositions(KmerPositions * kmerPositions, char * nameC, char * sequenceC, int kSize);

void deleteAllKmerPositions(KmerPositions * kmerPositions);


#endif // KMERS_H

