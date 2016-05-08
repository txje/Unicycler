

enum CigarType {MATCH, INSERTION, DELETION, CLIP, NOTHING};


class Alignment {
public:
    Alignment(Align<Dna5String, ArrayGaps> & alignment, int refOf8fset, long long startTime,
              bool startImmediately, bool goToEnd, ScoringScheme & scoringScheme);
    std::string getFullString();
    std::string getShortDisplayString();

    int m_readStartPos;
    int m_readEndPos;
    int m_refStartPos;
    int m_refEndPos;
    std::string m_cigar;
    int m_rawScore;
    double m_scaledScore;
    int m_milliseconds;

private:
    CigarType getCigarType(char b1, char b2, bool alignmentStarted);
    std::string getCigarPart(CigarType type, int length);
    int getCigarScore(CigarType type, int length, ScoringScheme & scoringScheme);
};

