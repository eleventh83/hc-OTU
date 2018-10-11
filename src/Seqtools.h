#ifndef SEQTOOLS_H_
#define SEQTOOLS_H_

class Seqtools{
public:
	Seqtools();
	~Seqtools();

	static int isMatchedDegenerate(char base1, char base2);
	static double calSimilarity(std::string seq1, std::string seq2);
	static std::string condense(std::string seq);
	static std::string unAlign(std::string seq);
	static std::string padEndGaps(std::string, char fill);

};

#endif /* SEQTOOLS_H_ */
