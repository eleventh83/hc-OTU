#ifndef CLSQFASTAREADER_H_
#define CLSQFASTAREADER_H_

class CLSQFastaReader {
public:
	CLSQFastaReader(std::string filename);
	~CLSQFastaReader();

	std::string fileName;
	std::vector<CLSQSequence> clSQList;

	void read();
	std::vector<CLSQSequence> getList();
};


#endif /* CLSQFASTAREADER_H_ */
