#include <iostream>
#include <string>
#include <vector>

#include "FastaSeq.h"
#include "Contig.h"
#include "Tools.h"
#include "AlignNoGap.h"




// MatchMismatch

MatchMismatch::MatchMismatch(){
	match = 0;
	mismatch = 0;
}

MatchMismatch::~MatchMismatch(){}

void MatchMismatch::init(){
	match = 0;
	mismatch = 0;
}





// AlignNoGap

AlignNoGap::AlignNoGap(){
	init();
	minOverlap = 50;
}

AlignNoGap::AlignNoGap(int th_overlap){
	init();
	minOverlap = th_overlap;
}

AlignNoGap::AlignNoGap(std::string seq1, std::string seq2, int th_overlap){
	init();
	minOverlap = th_overlap;
	s1 = seq1;
	s2 = seq2;
}

AlignNoGap::~AlignNoGap(){}

void AlignNoGap::init(){
	posStart1 = -1;
	posStart2 = -1;
	match = 0;
	max_allowed = 0;
	isAligned = false;
}

void AlignNoGap::setSequences(std::string seq1, std::string seq2){
	init();
	s1 = seq1;
	s2 = seq2;
}

void AlignNoGap::setMinimumOverlap(int v){
	minOverlap = v;
}

int AlignNoGap::compTwo(int start1, int start2){
	int match = 0;
	for (;;){
		if(s1[start1++]!=s2[start2++])	return 0;
		match++;
		if(start1 >= (signed)s1.size())	break;
		if(start2 >= (signed)s2.size())	break;
	}
	return match;
}


MatchMismatch AlignNoGap::compTwo_mistch_allow(int start1, int start2, int max_allowed){
	MatchMismatch mm;
	for(;;){
		if(s1[start1] != s2[start2])	mm.mismatch++;
		else	mm.match++;

		if(mm.mismatch > max_allowed) return mm;
		start1++;
		start2++;
		if(start1 >= (signed)s1.size())	break;
		if(start2 >= (signed)s2.size())	break;
	}
	return mm;
}

int AlignNoGap::doAlign_mistch_allow(int max_allowed){
	MatchMismatch mm;
	int i=0;
	while(i<(signed)s2.size()-minOverlap){
		mm = compTwo_mistch_allow(0, i, max_allowed);
		if(mm.mismatch>max_allowed){
			if(mm.match>0)	i+= mm.match;
			else	i++;
			continue;
		}
		if(mm.match>minOverlap){
			posStart1=0;
			posStart2=i;
			isAligned = true;
			return mm.mismatch;
		}
		else	break;
	}

	i=0;
	while(i<(signed)s1.size()-minOverlap){
		mm = compTwo_mistch_allow(i, 0, max_allowed);
		if(mm.mismatch>max_allowed){
			if(mm.match>0)	i += mm.match;
			else	i++;
			continue;
		}
		if(mm.match>minOverlap){
			posStart1=i;
			posStart2=0;
			isAligned = true;
			return mm.mismatch;
		}
		else	break;
	}
	return -1;
}


bool AlignNoGap::doAlign(){
	for(int i=0; i<(signed)s2.size(); i++){
		match = compTwo(0,i);
		if(match>minOverlap){
			posStart1=0;
			posStart2=i;
			return	true;
		}
	}
	for(int i=0; i<(signed)s1.size(); i++){
		match = compTwo(i,0);
		if(match>minOverlap){
			posStart1=i;
			posStart2=0;
			return	true;
		}
	}
	match = 0;
	return	false;
}

std::string *AlignNoGap::getAlignedSequence(){
	std::string* result = new std::string[2];
	result[0] = s1;
	result[1] = s2;
	if(posStart1==0)	result[0] = Tools::padString(result[0],posStart2,'~',true); 
	else	result[1] = Tools::padString(result[1],posStart1,'~',true);;
	return result;
}


void AlignNoGap::printAlignment(){
	std::string* seq;
	seq = getAlignedSequence();
	std::cout << "Seq1	" << seq[0] << "\n";
	std::cout << "Seq2	" << seq[1] << "\n";
	delete [] seq;
}

std::string AlignNoGap::getSumSeq(){

	//if(!isAligned)	return "";

	std::string sum;
	int len1 = s1.size();
	int len2 = s2.size();
	int initOverlapPos = posStart1 + posStart2;

	if(posStart1 == 0 )	len1 += posStart2;
	else	len2 += posStart1;
	if(len1>len2)	sum.resize(len1);
	else	sum.resize(len2);

	for(int i=0; i<initOverlapPos; i++){
		if(posStart1==0)	sum[i] = s2[i];
		else	sum[i] = s1[i];
	}
	for(int i=initOverlapPos; i<match+initOverlapPos; i++){
		if(posStart1==0)	sum[i] = s2[i];
		else	sum[i] = s1[i];
	}
	for(int i=match+initOverlapPos; i<(signed)sum.size(); i++){
		if(len1>len2){
			if(posStart1==0)	sum[i]=s1[i-initOverlapPos];
			else	sum[i]=s1[i];
		}else{
			if(posStart2==0)	sum[i]=s2[i-initOverlapPos];
			else	sum[i]=s2[i];
		}
	}
	return sum;
}

std::string AlignNoGap::getSumSeqMod(){
	std::string sum;
	int len1 = s1.size();
	int len2 = s2.size();
	int tlen1 = len1;
	int tlen2 = len2;

	if(posStart1==0)	tlen1 += posStart2;
	else	tlen2 += posStart2;
	if(tlen1>tlen2)	sum.resize(tlen1);
	else	sum.resize(tlen2);

	if(posStart1 == 0 && posStart2 == 0){
		if(len1 > len2){
			for(int i=0; i<len1; i++)	sum[i] = s1[i];
			return sum;
		}else{
			for(int i=0; i<len2; i++)	sum[i] = s2[i];
			return sum;
		}
	}

	if(posStart1 == 0){
		for(int i=0; i<len2; i++)	sum[i] = s2[i];
		int cnt = len2;
		if(len2 - posStart2 < len1)
			for(int i=len2-posStart2; i<len1; i++)	sum[cnt++] = s1[i];
	}
	else{
		for(int i=0; i<len1; i++)	sum[i] = s1[i];
		int cnt = len1;
		if(len1 - posStart1 < len2)
			for(int i=len1-posStart1; i<len2; i++)	sum[cnt++] = s2[i];
	}
	return	sum;
}


std::string AlignNoGap::getSumSeqMod2(){
	std::string sum;
	int len1 = s1.size();
	int len2 = s2.size();
	int tlen1 = len1;
	int tlen2 = len2;

	if(posStart1==0) tlen1 += posStart2;
	else	tlen2 += posStart2;
	if(tlen1>tlen2)	sum.resize(tlen1);
	else	sum.resize(tlen2);

	if(posStart1 == 0 && posStart2 == 0){
		if(len1 > len2){
			sum = s1;
			return sum;
		}else{
			sum = s2;
			return sum;
		}
	}

	if(posStart1 == 0){
		sum = s2;
		if(len2 - posStart2 < len1)
			sum += s1.substr(len2-posStart2,len1-len2+posStart2-1);
	}else{
		sum = s1;
		if(len1 - posStart1 < len2)
			sum += s2.substr(len1-posStart1,len2-len1+posStart1-1);
	}
	return	sum;
}


