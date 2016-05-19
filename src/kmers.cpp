
#include "kmers.h"




// Creates the CommonKmer object using the position in the two sequences.
// It then rotates the point relative to the origin by the given angle (in radians).
CommonKmer::CommonKmer(int hPosition, int vPosition, float angle) :
    m_hPosition(hPosition),
    m_vPosition(vPosition),
    m_score(1.0) {
    float s = sin(angle);
    float c = cos(angle);
    m_rotatedHPosition = (m_hPosition * c) - (m_vPosition * s);
    m_rotatedVPosition = (m_hPosition * s) + (m_vPosition * c);
}


// This is the destructor for KmerPositions. It cleans up all the KmerPosMaps which were allocated
// on the heap.
KmerPositions::~KmerPositions() {
    for (std::unordered_map<std::string, KmerPosMap *>::iterator i = m_kmerPositions.begin(); i != m_kmerPositions.end(); ++i)
        delete i->second;
}


// Returns a vector all of k-mer position names (should be exact the same as the the sequence names).
std::vector<std::string> KmerPositions::getAllNames() {
    std::vector<std::string> returnVector;
    m_mutex.lock();
    for (std::unordered_map<std::string, KmerPosMap *>::iterator i = m_kmerPositions.begin(); i != m_kmerPositions.end(); ++i)
        returnVector.push_back(i->first);
    m_mutex.unlock();
    return returnVector;
}

// Returns the length of the sequence with the given name.
int KmerPositions::getLength(std::string & name) {
    int returnedLength = 0;
    m_mutex.lock();
    if (m_sequences.find(name) != m_sequences.end())
        returnedLength = m_sequences[name].length();
    m_mutex.unlock();
    return returnedLength;
}

// This function adds a sequence to the KmerPositions object. It creates a new KmerPosMap on the
// heap (will be deleted in destructor), fills it up and adds it to m_kmerPositions.
void KmerPositions::addPositions(std::string & name, std::string & sequence) {
    KmerPosMap * posMap = new KmerPosMap();
    int kCount = sequence.size() - KMER_SIZE + 1;
    for (int i = 0; i < kCount; ++i) {
        std::string kmer = sequence.substr(i, KMER_SIZE);
        if (posMap->find(kmer) == posMap->end())
            (*posMap)[kmer] = std::vector<int>();
        (*posMap)[kmer].push_back(i);
    }

    m_mutex.lock();
    m_sequences[name] = sequence;
    m_kmerPositions[name] = posMap;
    m_mutex.unlock();
}

void KmerPositions::deletePositions(std::string & name) {
    m_mutex.lock();
    if (m_sequences.find(name) != m_sequences.end())
        m_sequences.erase(name);
    m_mutex.unlock();

    KmerPosMap * kmerPosMap = getKmerPositions(name);
    
    m_mutex.lock();
    if (kmerPosMap != 0) {
        m_kmerPositions.erase(name);
        delete kmerPosMap;
    }
    m_mutex.unlock();
}

// This function retrieves a KmerPosMap from the object using the name as a key. If the name isn't
// in the map, it returns 0.
KmerPosMap * KmerPositions::getKmerPositions(std::string & name) {
    KmerPosMap * returnedMap = 0;
    m_mutex.lock();
    if (m_kmerPositions.find(name) != m_kmerPositions.end())
        returnedMap =  m_kmerPositions[name];
    m_mutex.unlock();
    return returnedMap;
}


std::string * KmerPositions::getSequence(std::string & name) {
    std::string * returnedSequence = 0;
    m_mutex.lock();
    if (m_kmerPositions.find(name) != m_kmerPositions.end())
        returnedSequence = &(m_sequences[name]);
    m_mutex.unlock();
    return returnedSequence;
}


KmerPositions * newKmerPositions() {
    return new KmerPositions();
}

void addKmerPositions(KmerPositions * kmerPositions, char * nameC, char * sequenceC) {
    std::string name(nameC);
    std::string sequence(sequenceC);
    kmerPositions->addPositions(name, sequence);
}

// void deleteKmerPositions(KmerPositions * kmerPositions, char * name) {
//     kmerPositions->deletePositions(name);
// }

void deleteAllKmerPositions(KmerPositions * kmerPositions) {
    delete kmerPositions;
}


