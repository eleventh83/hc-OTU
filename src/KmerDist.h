#ifndef KMERDIST_H_
#define KMERDIST_H_

class KmerDist{
public:
	KmerDist();
	KmerDist(int k);
	KmerDist(std::string seqA, std::string seqB, int k);
	~KmerDist();


	std::vector<int> numeric_seqA;
	std::vector<int> numeric_seqB;
	std::vector<int> freqList_A;
	std::vector<int> freqList_B;
	int kmer;

	void setSeqA(std::string seq);
	void setSeqB(std::string seq);
	void setFreqA(std::string seq);
	void setFreqB(std::string seq);
	static std::vector<int> convseq(std::string seq);
	std::vector<int> getFreq(std::vector<int> numeric_seq);
	static std::vector<int> getFreqFromNT(std::string seq, int kmer);
	double getsim();
	static double getsim(std::vector<int> freqList_A, std::vector<int> freqList_B, int kmer);
};

#endif /* KMERDIST_H_ */







