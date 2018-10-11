#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "FastaSeq.h"
#include "Contig.h"
#include "ContigCluster.h"
#include "KmerDist.h"
#include "KmerHitLine.h"
#include "KmerSimple.h"



bool compareKmerHitLine(KmerHitLine a, KmerHitLine b){
	return a.kmerscore < b.kmerscore;
}


KmerSimple::KmerSimple(){
	ksize = 5;	// default word size
	//querySeq = "";
	//seqHeader = "";
	numSeq = 0;
	queryIdx = 0;
}

KmerSimple::KmerSimple(int qIdx, std::vector<ContigCluster> ctgList){
	ksize = 5;
	//querySeq = "";
	//seqHeader = "";
	numSeq = 0;

	queryIdx = qIdx;
	contigList = ctgList;
}

KmerSimple::~KmerSimple(){}


/*
std::vector<int> KmerSimple::runkmer(int k){
	std::vector<KmerHitLine> hlk;
	std::vector<KmerHitLine> t_hlk(contigList.size());
	std::string querySeq = contigList[queryIdx].seq;

	//#pragma omp parallel for schedule(auto)
	for(int i=0; i<contigList.size(); i++){
		if((contigList[i].noClones == 0) || contigList[i].solo) continue;

		KmerDist kmd(k);
		kmd.setFreqA(querySeq);
		kmd.setFreqB(contigList[i].seq);
		double dist = kmd.getsim();

		if(dist < 0.30) {
			t_hlk[i].setIdx(i);
			t_hlk[i].setScore(dist);
		}
	}

	for(int i=0; i<t_hlk.size(); i++){		
			if(t_hlk[i].kmerscore < 0.30) hlk.push_back(t_hlk[i]);
	}


	std::sort(hlk.begin(), hlk.end(),compareKmerHitLine);


	hitLine.resize(hlk.size());
	//#pragma omp parallel for
	for(int i=0; i<hlk.size(); i++){
		hitLine[i] = hlk[i].contigidx;
	}

	return hitLine;
}
*/



std::vector<int> KmerSimple::runkmer(int k){
	std::vector<KmerHitLine> hlk;
	std::vector<KmerHitLine> t_hlk(contigList.size());
	std::vector<int> querySeqFreq = contigList[queryIdx].kmercount;

	/*
	std::vector<int> targetCtg;
	for(int i=0; i<contigList.size(); i++){
			if((contigList[i].noClones != 0) && !contigList[i].solo)
					targetCtg.push_back(i);
	}
	#pragma omp parallel for schedule(auto)
	for(int i=0; i<targetCtg.size(); i++){
			double dist = KmerDist::getsim(querySeqFreq, contigList[targetCtg[i]].kmercount, 5);

			if(dist < 0.30){
					t_hlk[i].setIdx(targetCtg[i]);
					t_hlk[i].setScore(dist);
			}
	}
	*/

	
	#pragma omp parallel for schedule(dynamic)
	for(int i=0; i<contigList.size(); i++){
		if((contigList[i].noClones == 0) || contigList[i].solo) continue;

		double dist = KmerDist::getsim(querySeqFreq, contigList[i].kmercount, 5);

		if(dist < 0.30) {
			t_hlk[i].setIdx(i);
			t_hlk[i].setScore(dist);
			/*
				#pragma omp critical (UPDATE)
				{
					hlk.push_back(KmerHitLine(i,dist));
				}
				*/
		}
	}
	
	
	for(int i=0; i<t_hlk.size(); i++){
		if(t_hlk[i].kmerscore < 0.30) hlk.push_back(t_hlk[i]);
	}
	

	std::sort(hlk.begin(), hlk.end(),compareKmerHitLine);


	hitLine.resize(hlk.size());
	//#pragma omp parallel for schedule(static)
	for(int i=0; i<hlk.size(); i++){
		hitLine[i] = hlk[i].contigidx;
	}

	return hitLine;
}



