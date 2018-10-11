#ifndef KMERHITLINE_H_
#define KMERHITLINE_H_

class KmerHitLine{
public:
	KmerHitLine();
	KmerHitLine(int idx, double score);
	~KmerHitLine();

	int contigidx;
	double kmerscore;

	void setIdx(int idx);
	void setScore(double score);

	//int compareTo(kmerHitLine o);	// not done !!
};

#endif /* KMERHITLINE_H_ */
