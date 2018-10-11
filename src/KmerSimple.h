#ifndef KMERSIMPLE_H_
#define KMERSIMPLE_H_

class KmerSimple{
public:
	KmerSimple();
	KmerSimple(int qIdx, std::vector<ContigCluster> ctgList);
	~KmerSimple();

	int ksize;
	int numSeq;
	int queryIdx;
	std::vector<ContigCluster> contigList;
	std::vector<int> hitLine;
	
	std::vector<int> runkmer(int k);
};

#endif /* KMERSIMPLE_H_ */



