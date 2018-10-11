#include <iostream>
#include <string>

#include "Fasta.h"
#include "Fastq.h"


Fastq::Fastq(){}
Fastq::~Fastq(){}


std::string Fastq::getQuality(){
	return quality;
}

void Fastq::setQuality(std::string quality){
	this->quality = quality;
}

void Fastq::print(){
	std::cout << "Name " + title << "\n";
	std::cout << "Seqs " + sequence << "\n";
	std::cout << "Qual " + quality << "\n";
}


