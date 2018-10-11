#ifndef ALIGNNOGAP_H_
#define ALIGNNOGAP_H_


class MatchMismatch {
public:
	MatchMismatch();
	~MatchMismatch();

	int match;
	int mismatch;

	void init();
};


class AlignNoGap {
public:
	AlignNoGap();
	AlignNoGap(int th_overlap);
	AlignNoGap(std::string seq1, std::string seq2, int th_overlap);
	~AlignNoGap();

	
	bool isAligned;
	int posStart1;
	int posStart2;
	int match;
	int max_allowed;
	std::string s1,s2;
	int minOverlap;


	void init();
	void setSequences(std::string seq1, std::string seq2);
	void setMinimumOverlap(int v);
	int compTwo(int start1, int start2);
	MatchMismatch compTwo_mistch_allow(int start1, int start2, int max_allowed);
	int doAlign_mistch_allow(int max_allowed);
	bool doAlign();
	std::string *getAlignedSequence();
	void printAlignment();
	std::string getSumSeq();
	std::string getSumSeqMod();
	std::string getSumSeqMod2();

};



#endif /* ALIGNOGAP_H_ */
