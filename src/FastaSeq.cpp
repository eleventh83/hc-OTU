#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>

#include "generalFunctions.h"
#include "FastaSeq.h"


FastaSeq::FastaSeq(){
	title = "";
	acc = "";
	gi = -1;
	sequence = "";
	quality = "";
	noAddedSeq = 0;
	length = -1;
	cluster = -1;
	index = -1;
	coverage = -1.0;
	n_reads = -1;
	scaffold = "";
}

FastaSeq::FastaSeq(const FastaSeq &fa){
	title = fa.title;
	acc = fa.acc;
	gi = -1;
	sequence = fa.sequence;
	quality = "";
	noAddedSeq = fa.noAddedSeq;
	length = -1;
	cluster = -1;
	index = -1;
	coverage = -1.0;
	n_reads = -1;
	scaffold = "";
}

FastaSeq::FastaSeq(std::string title, std::string sequence){
	this->title = title;
	this->acc = "";
	gi = -1;
	this->sequence = sequence;
	quality = "";
	this->noAddedSeq = 0;
	length = -1;
	cluster = -1;
	index = -1;
	coverage = -1.0;
	n_reads = -1;
	scaffold = "";
}

FastaSeq::FastaSeq(std::string title, std::string sequence, int noAddedSeq){
	this->title = title;
	this->acc = "";
	gi = -1;
	this->sequence = sequence;
	quality = "";
	this->noAddedSeq = noAddedSeq;
	length = -1;
	cluster = -1;
	index = -1;
	coverage = -1.0;
	n_reads = -1;
	scaffold = "";
}

FastaSeq::~FastaSeq(){}


//int FastaSeq::compareTo(FastaSeq o){}


void FastaSeq::setSeq(std::string title, std::string sequence){
	this->title = title;
	this->sequence = sequence;
	this->acc = "";
}

void FastaSeq::setSeq(std::string title, std::string sequence, int noAddedSeq){
	this->title = title;
	this->sequence = sequence;
	this->noAddedSeq = noAddedSeq;
	this->acc = "";
}


std::string FastaSeq::getFastaFormat(){
	std::string str = ">" + title + "\n" + sequence;
	return  str;
}


void FastaSeq::print(){
	std::cout << getFastaFormat() + "\n";
}


void FastaSeq::writeAsBlastInfile(std::string fileName){
	std::ofstream fw;
	fw.open(fileName.c_str());
	fw << ">" + title + "\n" + sequence;
	fw.close();
}


void FastaSeq::parseNcbi(std::string title){
	this->title = title;
	if(title == "") return;
	std::vector<std::string> s = StringSplit(title,"\\|",5);
	if(s.size() > 1){
		if(s[0] == "gi")	gi = atoi(s[1].c_str());
		if(s[0] == ">gi")	gi = atoi(s[1].c_str());
	}
	if(s.size() > 3){
		if(toUpperCase(s[2]) == "GB")	acc = s[3];
		if(toUpperCase(s[2]) == "EMB")	acc = s[3];
		if(toUpperCase(s[2]) == "DBJ")	acc = s[3];
		if(toUpperCase(s[2]) == "REF")	acc = s[3];
		if(toUpperCase(s[2]) == "TPG")	acc = s[3];
	}
}


