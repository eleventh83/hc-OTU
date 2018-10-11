#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <omp.h>

#include "KmerDist.h"


KmerDist::KmerDist(){
	kmer = 0;
}

KmerDist::KmerDist(int k){
	kmer = k;
}

KmerDist::KmerDist(std::string seqA, std::string seqB, int k){
	numeric_seqA = convseq(seqA);
	numeric_seqB = convseq(seqB);
	kmer = k;
}

KmerDist::~KmerDist(){}

void KmerDist::setSeqA(std::string seq){
	numeric_seqA = convseq(seq);
}

void KmerDist::setSeqB(std::string seq){
	numeric_seqB = convseq(seq);
}

void KmerDist::setFreqA(std::string seq){
	numeric_seqA = convseq(seq);
	freqList_A = getFreq(numeric_seqA);
}

void KmerDist::setFreqB(std::string seq){
	numeric_seqB = convseq(seq);
	freqList_B = getFreq(numeric_seqB);
}

std::vector<int> KmerDist::convseq(std::string seq){
	std::vector<int> numeric_seq(seq.size());
//#pragma omp parallel for
	for(int i=0; i<seq.size(); i++){
		switch(seq.at(i)){
		case 'A':
			numeric_seq[i] = 1;
			break;
		case 'C':
			numeric_seq[i] = 2;
			break;
		case 'G':
			numeric_seq[i] = 3;
			break;
		default:
			numeric_seq[i] = 4;
			break;
		}
	}
	return numeric_seq;
}

std::vector<int> KmerDist::getFreq(std::vector<int> numeric_seq){
	std::vector<int> freqList((int)pow((double)4,kmer),0);
	int idx = 0;
	for(int i=0; i<numeric_seq.size()-kmer+1; i++){
		int idx = 0;
		//for(int j=i; j<i+kmer; j++) idx += (numeric_seq[j]-1)*(int)pow((double)4,kmer-j+i-1);
		for(int j=i; j<i+kmer; j++){
			int pow_factor = 1;
			for(int k=0; k<kmer-j+i-1; k++)	pow_factor *= 4;
			
			idx += (numeric_seq[j]-1)*pow_factor;
		}
		freqList.at(idx)++;
	}
	return freqList;
}


std::vector<int> KmerDist::getFreqFromNT(std::string seq, int kmer){
	std::vector<int> numeric_seq = KmerDist::convseq(seq);
	std::vector<int> freqList((int)pow((double)4,kmer),0);
	int idx = 0;
	for(int i=0; i<numeric_seq.size()-kmer+1; i++){
		int idx = 0;
		for(int j=i; j<i+kmer; j++){
			int pow_factor = 1;
			for(int k=0; k<kmer-j+i-1; k++)	pow_factor *= 4;
			
			idx += (numeric_seq[j]-1)*pow_factor;
		}
		freqList.at(idx)++;
	}
	return freqList;
}

double KmerDist::getsim(){
	double sumcnt = 0;
	for(int i=0; i<freqList_A.size(); i++){
		sumcnt += std::min(freqList_A[i], freqList_B[i]);
	}
	double normfactor = std::min(numeric_seqA.size(), numeric_seqB.size())-kmer+1;
	return 1-sumcnt/normfactor;
}



double KmerDist::getsim(std::vector<int> freqA, std::vector<int> freqB, int kmer){
	double sumcnt = 0;
	int N_A = 0;
	int N_B = 0;

	for(int i=0; i<freqA.size(); i++){
		sumcnt += std::min(freqA[i], freqB[i]);
		N_A += freqA[i];
		N_B += freqB[i];
	}

	double normfactor = std::min(N_A, N_B)-kmer+1;
	return 1-sumcnt/normfactor;
}



