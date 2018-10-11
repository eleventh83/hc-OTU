//#include <boost/numeric/ublas/vector.hpp>
#include <vector>
#include <string>

#include "Fasta.h"
#include "Fastq.h"
#include "CLSQSequence.h"


CLSQSequence::CLSQSequence(){
	noClones = 1;
	//contig = "";
	//condensed_seq = "";
	solo = false;
}

CLSQSequence::~CLSQSequence(){};


std::string CLSQSequence::getCondensedSeq(){
	return condensed_seq;
}

std::string CLSQSequence::getSequence(){
	return sequence;
}

void CLSQSequence::setCondensedSeq(std::string condensedseq){
	this->condensed_seq = condensedseq;
}

int CLSQSequence::getCloneNum(){
	return noClones;
}

void CLSQSequence::setCloneNum(int cloneNo){
	noClones = cloneNo;
}

void CLSQSequence::addClone(){
	noClones++;
}

void CLSQSequence::addCloneFor(int idx){
	clonesFor.push_back(idx);
}

std::string CLSQSequence::getContig(){
	return contig;
}

void CLSQSequence::setContig(std::string contig){
	this->contig = contig;
}

void CLSQSequence::setProfile(std::vector<int> profile){
		this->profile = profile;
}

void CLSQSequence::setKmercount(std::vector<int> kmercount){
		this->kmercount = kmercount;
}

void CLSQSequence::setSolo(bool single){
	solo = single;
}



//int CLSQSequence::compareTo(CLSQSequence o){	//not done
//	if(this->noClones < o.noClones) return 1;
//	if(this->noClones == o.noClones){
//}
