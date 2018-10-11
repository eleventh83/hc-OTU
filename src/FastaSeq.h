#ifndef FASTASEQ_H_
#define FASTASEQ_H_

class FastaSeq{
public:
	FastaSeq();
	FastaSeq(const FastaSeq &fa);
	FastaSeq(std::string title, std::string sequence);
	FastaSeq(std::string title, std::string sequence, int noAddedSeq);

	~FastaSeq();

	std::string title;
	std::string acc;
	int gi;
	std::string sequence;
	std::string quality;
	int noAddedSeq;
	int length;
	int cluster;
	int index;
	double coverage;
	int n_reads;
	std::string scaffold;

	//int compareTo(FastaSeq o);	// not done!!
	void setSeq(std::string title, std::string sequence);
	void setSeq(std::string title, std::string sequence, int noAddedSeq);
	std::string getFastaFormat();
	void print();
	void writeAsBlastInfile(std::string fileName);

	void parseNcbi(std::string title);
};

	

#endif /* FASTASEQ_H_ */
