#ifndef FASTA_H_
#define FASTA_H_

class Fasta {
protected:
	std::string title;
	std::string sequence;

public:
	Fasta();
	~Fasta();

	std::string getTitle();
	std::string getSequence();
	void setTitle(std::string title);
	void setSequence(std::string sequence);
	//FILE toFile(std::string fileName);
};

#endif /*FASTA_H_ */


