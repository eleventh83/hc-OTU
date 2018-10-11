#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "generalFunctions.h"
#include "Seqtools.h"
#include "Fasta.h"
#include "Fastq.h"
#include "CLSQSequence.h"
#include "CLSQFastaReader.h"


#define MAX_SIZE 1000

CLSQFastaReader::CLSQFastaReader(std::string filename){
	this->fileName = filename;
}

CLSQFastaReader::~CLSQFastaReader(){}

void CLSQFastaReader::read(){
	std::string tmpline;
	char inputString[MAX_SIZE];
	CLSQSequence tmpSQ;

	std::ifstream inFile(fileName.c_str());
	while(!inFile.eof()){
		inFile.getline(inputString,MAX_SIZE);
		if(inputString[0] == '>'){
			tmpline = inputString;
			tmpline = tmpline.substr(1);
			
			tmpSQ.setTitle(trim(tmpline));

			inFile.getline(inputString,MAX_SIZE);
			tmpSQ.setSequence(trim(inputString));

			tmpSQ.setCondensedSeq(Seqtools::condense(tmpSQ.getSequence()));

			clSQList.push_back(tmpSQ);
		}
		//std::cout << inputString << "\n";
	}
	inFile.close();
}

std::vector<CLSQSequence> CLSQFastaReader::getList(){
	return clSQList;
}


