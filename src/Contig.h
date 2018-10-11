#ifndef CONTIG_H_
#define CONTIG_H_

class Contig{

public:
	Contig();
	~Contig();

	std::string seq;
	int noClones;
	std::string title;
	std::vector<FastaSeq> seqList;
	std::vector<int> idxList;

	//int compareTo(Contig ctg);

};

#endif /* CONTIG_H_ */
