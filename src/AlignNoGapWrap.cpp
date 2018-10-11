#include <string>
#include <vector>


#include "Fasta.h"
#include "Fastq.h"
#include "CLSQSequence.h"
#include "FastaSeq.h"
#include "Contig.h"
#include "Tools.h"
#include "AlignNoGap.h"
#include "AlignNoGapWrap.h"


AlignNoGapWrap::AlignNoGapWrap(){
	s1 = "";
	minOverlap = 200;
	idxstart = 0;
	idxend = 0;
}

AlignNoGapWrap::AlignNoGapWrap(std::string qseq, std::vector<CLSQSequence> seqsList, int startidx, int endidx, int th_overlap){
	s1 = qseq;
	minOverlap = th_overlap;
	seqList = seqsList;
	idxstart = startidx;
	idxend = endidx;
}

AlignNoGapWrap::~AlignNoGapWrap(){}


void AlignNoGapWrap::setQseq(std::string qseq){
	s1 = qseq;
}

void AlignNoGapWrap::setRange(int startidx, int endidx){
	idxstart = startidx;
	idxend = endidx;
}

void AlignNoGapWrap::clear(){
	match_idx.clear();
}

void AlignNoGapWrap::setMinOverlap(int th_overlap){
	minOverlap = th_overlap;
}

void AlignNoGapWrap::setSeqList(std::vector<CLSQSequence> seqsList){
	seqList = seqsList;
}

std::vector<int> AlignNoGapWrap::call(){
	match_idx.reserve(idxend-idxstart+1);
	for(int i=idxstart; i<=idxend; i++){
		init();
		s2 = seqList.at(i).getCondensedSeq();
		int flag = doAlign_mistch_allow(max_allowed);
		if(flag != -1) match_idx.push_back(i);
	}
	return match_idx;
}
