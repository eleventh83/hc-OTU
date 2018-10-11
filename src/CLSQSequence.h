#ifndef CLSQSEQUENCE_H_
#define CLSQSEQUENCE_H_


class CLSQSequence:public Fastq {
public:
	CLSQSequence();
	~CLSQSequence();

	int noClones;
	std::string contig;
	std::string condensed_seq;
	bool solo;

	std::vector<int> clonesFor;
	std::vector<int> clonesRev;
	std::vector<int> profile;
	std::vector<int> kmercount;


	int compareTo(CLSQSequence o);
	std::string getSequence();
	std::string getCondensedSeq();
	void setCondensedSeq(std::string condensedseq);
	int getCloneNum();
	void setCloneNum(int cloneNo);
	void addClone();
	void addCloneFor(int idx);
	std::string getContig();
	void setContig(std::string contig);

	void setProfile(std::vector<int> profile);
	void setKmercount(std::vector<int> kmercount);
	void setSolo(bool single);
};

#endif /* CLSQSEQUENCE_H_ */


