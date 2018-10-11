#ifndef FASTQ_H_
#define FASTQ_H_


class Fastq:public Fasta {
private:
	std::string quality;

public:
	Fastq();
	~Fastq();

	std::string getQuality();
	void setQuality( std::string quality);
	void print();
};

#endif /* FASTQ_H_ */

