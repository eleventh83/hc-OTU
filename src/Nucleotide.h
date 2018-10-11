#ifndef NUCLEOTIDE_H_
#define NUCLEOTIDE_H_

class Nucleotide{
public:
	Nucleotide();
	~Nucleotide();

	int A;
	int C;
	int G;
	int T;
	int Gap;
	int endGap;
	int N;
	int total;

	char getMajority_pyro();
	char getMajority_normal();
};

#endif /* NUCLEOTIDE_H_ */
