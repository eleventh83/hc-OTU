#ifndef ASSEMBLEGLOBAL_H_
#define ASSEMBLEGLOBAL_H_

class AssembleGlobal{
public:
	AssembleGlobal(int initialCapacity);
	AssembleGlobal(std::vector<std::string> inSeq);
	~AssembleGlobal();

	std::string contig;
	std::vector<std::string> seq;
	std::vector<Nucleotide> nt;
	std::vector<int> Outlier;
	
	std::string getContig();
	std::string getContigWithoutGap();
	void addSeq(std::string oneSeq);
	static std::string SumTwoSeq(std::string seq1, std::string seq2);
	void doAssemble();
	void doConsensus();
	void print();
};

#endif /* ASSEMBLEGLOBAL_H_ */


