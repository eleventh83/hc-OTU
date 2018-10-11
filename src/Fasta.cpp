#include <string>

#include "Fasta.h"

Fasta::Fasta(){}
Fasta::~Fasta(){}

std::string Fasta::getTitle(){
	return title;
}

std::string Fasta::getSequence(){
	return sequence;
}

void Fasta::setTitle(std::string title){
	this->title = title;
}

void Fasta::setSequence(std::string sequence){
	this->sequence = sequence;
}

//FILE Fasta::toFile(std::string fileName){	// not done
//	FILE tmpFile;
//	return tmpFile;
//}
