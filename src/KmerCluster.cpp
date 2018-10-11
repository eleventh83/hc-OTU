#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <sys/time.h>

#include "generalFunctions.h"
#include "FastaSeq.h"
#include "Contig.h"
#include "Tools.h"
#include "ContigCluster.h"
#include "KmerDist.h"
#include "KmerHitLine.h"
#include "KmerSimple.h"
#include "Nucleotide.h"
#include "Seqtools.h"
#include "needle.h"
#include "ssw.h"
#include "ssw_cpp.h"
#include "AssembleGlobal.h"
#include "KmerCluster.h"



bool compareKmer(KmerHitLine a, KmerHitLine b){	return a.kmerscore < b.kmerscore; }
bool compareCC(ContigCluster a, ContigCluster b){ return a.noClones > b.noClones; }

static inline double myclock(){
		struct timeval tval;
		gettimeofday(&tval, NULL);
		return (tval.tv_sec *1000 + tval.tv_usec/1000);
}



KmerCluster::KmerCluster(){
	temp_dir = "./tmp";
	sim_cutoff = 0.97;
	sim_cutoff_margin = 0.03;
	solo_cutoff = sim_cutoff - 0.05;
	totalClones = 0;
	totalOTUs = 0;
	//double stemp = Tools.getTimeStampNow();
	//int blast_type = Common.BLAST_TYPE_BLASTN;		// 1�̴�
	int blast_type = 1;
	std::string blast_path = "";
}

KmerCluster::KmerCluster(double sim_cutoff){
	temp_dir = "./tmp";
	this->sim_cutoff = sim_cutoff;
	sim_cutoff_margin = 0.03;
	solo_cutoff = sim_cutoff - 0.05;
	totalClones = 0;
	totalOTUs = 0;
	//double stemp = Tools.getTimeStampNow();
	//int blast_type = Common.BLAST_TYPE_BLASTN;
	//std::string blast_path = "";
}

KmerCluster::~KmerCluster(){};

void KmerCluster::runContig(std::vector<Contig> ctgList, double sim_cutoff, double sim_cutoff_margin, std::string temp_dir, int blast_type, std::string blast_path){
	this->sim_cutoff = sim_cutoff;
	this->solo_cutoff = this->sim_cutoff-0.05;
	this->sim_cutoff_margin = sim_cutoff_margin;
	this->temp_dir = temp_dir;
	this->blast_type = blast_type;
	this->blast_path = blast_path;
	runContig(ctgList);
}





// runContig (serial) java

void KmerCluster::runContig(std::vector<Contig> ctgList){
	for(int i=0; i<ctgList.size(); i++){
		ContigCluster cc;
		cc.noClones = ctgList.at(i).noClones;
		cc.seq = ctgList.at(i).seq;
		cc.contigList.push_back(ctgList.at(i));
		cc.ContigIdxList = ctgList.at(i).idxList;
		ccList.push_back(cc);
	}




	KmerDist kmerD(5);
	for(int i=0; i<ccList.size();i++){
		if(ccList.at(i).noClones == 0)	continue;

		kmerD.setFreqA(ccList.at(i).seq);
		std::vector<int> hitline;
		std::vector<KmerHitLine> hlk;
		double ksim = 0.0;
		for(int j=0; j<ccList.size();j++){
			if(ccList.at(j).noClones == 0 || ccList.at(i).solo)	continue;
			kmerD.setFreqB(ccList.at(j).seq);
			ksim = kmerD.getsim();
			if(ksim < 0.15){
				KmerHitLine tKHL(j, ksim);
				hlk.push_back(tKHL);
			}
		}
		std::sort(hlk.begin(), hlk.end(),compareKmer);

		for(int j=0; j<hlk.size(); j++){
			hitline.push_back(hlk.at(j).contigidx);
		}

		Needle nwa(8,1);
		int t_len = ccList.at(i).seq.size() * 2;
		char* aligned1 = new char[t_len];
		char* aligned2 = new char[t_len];

		double sim = 0.0;
		int hit = 0;
		bool self_found = false;
		for(int n=0; n<hitline.size(); n++){
			hit = hitline.at(n);

			if(i==hit){
				self_found = true;
				continue;
			}

			nwa.Align(ccList.at(i).seq.c_str(),ccList.at(hit).seq.c_str(),aligned1, aligned2);
			sim = Tools::calculate(aligned1,aligned2);

			if(sim >= sim_cutoff){
				if(hit<i){
					ccList.at(hit).contigList.insert(ccList.at(hit).contigList.end(), ccList.at(i).contigList.begin(),ccList.at(i).contigList.end());
					ccList.at(hit).noClones += ccList.at(i).noClones;
					ccList.at(i).noClones = 0;
					ccList.at(hit).ContigIdxList.insert(ccList.at(hit).ContigIdxList.end(), ccList.at(i).ContigIdxList.begin(), ccList.at(i).ContigIdxList.end());
					break;
				}
				else{
					ccList.at(i).contigList.insert(ccList.at(i).contigList.end(), ccList.at(hit).contigList.begin(), ccList.at(hit).contigList.end());
					ccList.at(i).noClones += ccList.at(hit).noClones;
					ccList.at(hit).noClones = 0;
					ccList.at(i).ContigIdxList.insert(ccList.at(i).ContigIdxList.end(), ccList.at(hit).ContigIdxList.begin(), ccList.at(hit).ContigIdxList.end());
				}
			}
			else{
				if(sim<solo_cutoff){
					if(self_found && n==1){
						ccList.at(i).solo = true;
						break;
					}
				}
				if(sim < (sim_cutoff - sim_cutoff_margin))	break;
			}
		}
	}


	totalClones = 0;
	totalOTUs = 0;

	for(int i=0; i<ccList.size(); i++){
		totalClones += ccList.at(i).noClones;
		if(ccList.at(i).noClones > 0) totalOTUs++;
	}


	std::cout << "\nOTUS = " + int2str(totalOTUs);
	std::cout << "\nClones = " + int2str(totalClones) + "\n";


}



// Parallel

void KmerCluster::runContig_p1(std::vector<Contig> ctgList){

	ccList.resize(ctgList.size());
	for(int i=0; i<ctgList.size(); i++){
		ccList[i].noClones = ctgList[i].noClones;
		ccList[i].seq = ctgList[i].seq;
		ccList[i].contigList.push_back(ctgList[i]);
		ccList[i].ContigIdxList = ctgList[i].idxList;
		ccList[i].kmercount = KmerDist::getFreqFromNT(ccList[i].seq,5);
	}

	std::vector<int> hitLine;
	for(int i=0; i<ccList.size()-1; i++){
		if(ccList.at(i).noClones == 0) continue;
		
		// parallelized
		KmerSimple ks(i,ccList);
		hitLine = ks.runkmer(5);
		
		if(hitLine.size() == 0)	continue;
		
		int lastHitIdx = getValHitIdx(hitLine, ccList, i);
		
		if(lastHitIdx == 0 && hitLine[0] == i){
			ccList[i].solo = true;
			continue;
		}
		
		for(int j=0; j<=lastHitIdx; j++){
			int idx_hit = hitLine[j];
			if(idx_hit < i){
				ccList[idx_hit].contigList.insert(ccList[idx_hit].contigList.end(), ccList[i].contigList.begin(), ccList[i].contigList.end());
				ccList[idx_hit].noClones += ccList[i].noClones;
				ccList[idx_hit].ContigIdxList.insert(ccList[idx_hit].ContigIdxList.end(), ccList[i].ContigIdxList.begin(),ccList[i].ContigIdxList.end());
				ccList[i].noClones = 0;
				ccList[i].contigList.clear();
				ccList[i].ContigIdxList.clear();
			}else if(idx_hit > i){
				ccList[i].contigList.insert(ccList[i].contigList.end(), ccList[idx_hit].contigList.begin(), ccList[idx_hit].contigList.end());
				ccList[i].noClones += ccList[idx_hit].noClones;
				ccList[i].ContigIdxList.insert(ccList[i].ContigIdxList.end(), ccList[idx_hit].ContigIdxList.begin(), ccList[idx_hit].ContigIdxList.end());
				ccList[idx_hit].noClones = 0;
				ccList[idx_hit].contigList.clear();
				ccList[idx_hit].ContigIdxList.clear();
			}
		}
	}

	totalClones = 0;
	totalOTUs = 0;

	for(int i=0; i<ccList.size(); i++){
		totalClones += ccList[i].noClones;
		if(ccList[i].noClones > 0) totalOTUs++;
	}
	
	std::sort(ccList.begin(), ccList.end(),compareCC);
	
	std::cout << "\nOTUS = " + int2str(totalOTUs);
	std::cout << "\nClones = " + int2str(totalClones) + "\n";
}



// Parallel - calculate whole similarities from Kmer filtering

void KmerCluster::runContig_all(std::vector<Contig> ctgList){

	ccList.resize(ctgList.size());
	for(int i=0; i<ctgList.size(); i++){
		ccList[i].noClones = ctgList[i].noClones;
		ccList[i].seq = ctgList[i].seq;
		ccList[i].contigList.push_back(ctgList[i]);
		ccList[i].ContigIdxList = ctgList[i].idxList;
		ccList[i].kmercount = KmerDist::getFreqFromNT(ccList[i].seq,5);
	}


	double BEGIN, END, eTIME;

	std::vector<int> hitLine;
	for(int i=0; i<ccList.size()-1; i++){
		if(ccList.at(i).noClones == 0) continue;

		// parallelized
		
		BEGIN = myclock();

		KmerSimple ks(i,ccList);
		hitLine = ks.runkmer(5);

		END = myclock();
		eTIME += (END-BEGIN)/1000;



		if(hitLine.size() == 0)	continue;

		
		std::vector<int> valHitIdx = getValHitIdxAll(hitLine, ccList, i);


		if(valHitIdx.size() == 1 && hitLine[0] == i){
			ccList[i].solo = true;
			continue;
		}

		for(int j=0; j<valHitIdx.size(); j++){
			int idx_hit = hitLine[valHitIdx.at(j)];
			if(idx_hit < i){
				ccList[idx_hit].contigList.insert(ccList[idx_hit].contigList.end(), ccList[i].contigList.begin(), ccList[i].contigList.end());
				ccList[idx_hit].noClones += ccList[i].noClones;
				ccList[idx_hit].ContigIdxList.insert(ccList[idx_hit].ContigIdxList.end(), ccList[i].ContigIdxList.begin(),ccList[i].ContigIdxList.end());
				ccList[i].noClones = 0;
				ccList[i].contigList.clear();
				ccList[i].ContigIdxList.clear();
			}else if(idx_hit > i){
				ccList[i].contigList.insert(ccList[i].contigList.end(), ccList[idx_hit].contigList.begin(), ccList[idx_hit].contigList.end());
				ccList[i].noClones += ccList[idx_hit].noClones;
				ccList[i].ContigIdxList.insert(ccList[i].ContigIdxList.end(), ccList[idx_hit].ContigIdxList.begin(), ccList[idx_hit].ContigIdxList.end());
				ccList[idx_hit].noClones = 0;
				ccList[idx_hit].contigList.clear();
				ccList[idx_hit].ContigIdxList.clear();
			}
		}
	}

	totalClones = 0;
	totalOTUs = 0;

	for(int i=0; i<ccList.size(); i++){
		totalClones += ccList[i].noClones;
		if(ccList[i].noClones > 0) totalOTUs++;
	}

	std::sort(ccList.begin(), ccList.end(),compareCC);

	std::cout << "\nOTUS = " + int2str(totalOTUs);
	std::cout << "\nClones = " + int2str(totalClones) + "\n";

	std::cout << "\nElapsed time for k-mer : " << eTIME << "\n";
}




// Serial (Revised)

void KmerCluster::runContig_rev1(std::vector<Contig> ctgList){
	ccList.resize(ctgList.size());
	for(int i=0; i<ctgList.size(); i++){
		ContigCluster cc;
		cc.noClones = ctgList.at(i).noClones;
		cc.seq = ctgList.at(i).seq;
		cc.contigList.push_back(ctgList.at(i));
		cc.ContigIdxList = ctgList.at(i).idxList;
		ccList.at(i) = cc;
	}

	std::cout << "Create ccList -> Done\n";

	std::vector<int> hitLine;
	for(int i=0; i<ccList.size()-1; i++){
		if(ccList.at(i).noClones == 0) continue;
		
		KmerSimple ks(i,ccList);
		hitLine = ks.runkmer(5);
		
		if(hitLine.size() == 0)	continue;
		
		int lastHitIdx = getValHitIdx(hitLine, ccList, i);
		
		if(lastHitIdx == 0 && hitLine.at(0) == i){
			ccList.at(i).solo = true;
			continue;
		}
		
		int idx_hit;
		for(int j=0; j<=lastHitIdx; j++){
			idx_hit = hitLine.at(j);
			if(idx_hit < i){
				ccList.at(idx_hit).contigList.insert(ccList.at(idx_hit).contigList.end(), ccList.at(i).contigList.begin(), ccList.at(i).contigList.end());
				ccList.at(idx_hit).noClones += ccList.at(i).noClones;
				// ccList.at(idx_hit).contigIdxList.. not done!!
				ccList.at(i).noClones = 0;
			}else if(idx_hit > i){
				ccList.at(i).contigList.insert(ccList.at(i).contigList.end(), ccList.at(idx_hit).contigList.begin(), ccList.at(idx_hit).contigList.end());
				ccList.at(i).noClones += ccList.at(idx_hit).noClones;
				//ccList.at(i).contigIdxList.. not done!!
				ccList.at(idx_hit).noClones = 0;
			}
		}
	}

	std::cout << "\nMerge -> Done\n";

	totalClones = 0;
	totalOTUs = 0;

	for(int i=0; i<ccList.size(); i++){
		totalClones += ccList.at(i).noClones;
		if(ccList.at(i).noClones > 0) totalOTUs++;
	}
	
	std::cout << "\nOTUS = " + int2str(totalOTUs);
	std::cout << "\nClones = " + int2str(totalClones) + "\n";
}




// serial
int KmerCluster::getValHitIdx(std::vector<int> hitLine, std::vector<ContigCluster> ccList, int qIdx){
	int hitLine_size = hitLine.size();
	std::vector<int> hitLineFlags(hitLine_size,-1);

	int current_idx = (int)(hitLine_size/2);
	int left_idx = 0;
	int right_idx = hitLine_size-1;
	double sim = 0.0;

	int t_len = ccList.at(0).seq.size()*2;
	

	// SSW
		
	std::string aligned1,aligned2;
	aligned1.reserve(t_len);
	aligned2.reserve(t_len);
	
	StripedSmithWaterman::Aligner sw(10,9,15,6);
	StripedSmithWaterman::Alignment al;
	StripedSmithWaterman::Filter filter;

	while(1){
		
		sw.Align(ccList.at(hitLine.at(current_idx)).seq.c_str(), ccList.at(qIdx).seq.c_str(), ccList.at(qIdx).seq.size(), filter, &al, aligned1, aligned2);
		sim = Tools::calculate(aligned1, aligned2);

		// flags, left, right update
		if(sim >= sim_cutoff){
			hitLineFlags.at(current_idx) = 1;
			left_idx = current_idx;
		}else{
			hitLineFlags.at(current_idx) = 0;
			right_idx = current_idx;
		}

		// current idx update
		if(left_idx+1 == right_idx){
			if(hitLineFlags.at(left_idx) == -1) current_idx = left_idx;
			else if(hitLineFlags.at(right_idx) == -1) current_idx = right_idx;
		}else{
			current_idx = (int)((left_idx+right_idx)/2);
		}

		// terminate condition
		if(hitLineFlags.at(0) == 0){
			current_idx = -1;
			break;
		}
		if(hitLineFlags.at(hitLine_size-1) == 1){
			current_idx = hitLine_size-1;
			break;
		}
		if(hitLineFlags.at(current_idx) == 1 && hitLineFlags.at(current_idx + 1) == 0) break;
		if(hitLineFlags.at(current_idx) == 0 && hitLineFlags.at(current_idx - 1) == 1){
			current_idx--;
			break;
		}
	}
	

	// Needle
	/*	
	char* aligned1 = new char[t_len];
	char* aligned2 = new char[t_len];
	Needle nwdist(10,0.1);

	while(1){
		nwdist.Align(ccList.at(qIdx).seq.c_str(), ccList.at(hitLine.at(current_idx)).seq.c_str(), aligned1, aligned2);
		sim = Tools::calculate(aligned1, aligned2);

		// flags, left, right update
		if(sim >= sim_cutoff){
			hitLineFlags.at(current_idx) = 1;
			left_idx = current_idx;
		}else{
			hitLineFlags.at(current_idx) = 0;
			right_idx = current_idx;
		}

		// current idx update
		if(left_idx+1 == right_idx){
			if(hitLineFlags.at(left_idx) == -1) current_idx = left_idx;
			else if(hitLineFlags.at(right_idx) == -1) current_idx = right_idx;
		}else{
			current_idx = (int)((left_idx+right_idx)/2);
		}

		// terminate condition
		if(hitLineFlags.at(0) == 0){
			current_idx = -1;
			break;
		}
		if(hitLineFlags.at(hitLine_size-1) == 1){
			current_idx = hitLine_size-1;
			break;
		}
		if(hitLineFlags.at(current_idx) == 1 && hitLineFlags.at(current_idx + 1) == 0) break;
		if(hitLineFlags.at(current_idx) == 0 && hitLineFlags.at(current_idx - 1) == 1){
			current_idx--;
			break;
		}
	}

	delete [] aligned1;
	delete [] aligned2;
	*/
	

	return current_idx;
}


std::vector<int> KmerCluster::getValHitIdxAll(std::vector<int> hitLine, std::vector<ContigCluster> ccList, int qIdx){
	int hitLine_size = hitLine.size();

	std::vector<int> valHitLine;
	double sim[hitLine_size];

	int t_len = ccList.at(0).seq.size()*2;

	// SSW
	
	#pragma omp parallel for schedule(static)
	for(int i=0; i<hitLine.size(); i++){
		
		std::string aligned1,aligned2;
		aligned1.reserve(t_len);
		aligned2.reserve(t_len);

		StripedSmithWaterman::Aligner sw(10,9,15,6);
		StripedSmithWaterman::Alignment al;
		StripedSmithWaterman::Filter filter;

		sw.Align(ccList.at(hitLine.at(i)).seq.c_str(), ccList.at(qIdx).seq.c_str(), ccList.at(qIdx).seq.size(), filter, &al, aligned1, aligned2);
		sim[i] = Tools::calculate(aligned1, aligned2);
	}
	

	/*
	#pragma omp parallel for schedule(auto)
	for(int i=0; i<hitLine.size(); i++){
			char* aligned1 = new char[t_len];
			char* aligned2 = new char[t_len];
			Needle nwdist(15,6.6);

			nwdist.Align(ccList.at(hitLine.at(i)).seq.c_str(), ccList.at(qIdx).seq.c_str(), aligned1, aligned2);
			sim[i] = Tools::calculate(aligned1,aligned2);

			delete [] aligned1;
			delete [] aligned2;
	}
	*/


	for(int i=0; i<hitLine.size(); i++){
		if(sim[i] >= sim_cutoff)	valHitLine.push_back(i);
		else continue;
	}

	return valHitLine;
}


// parallel
int KmerCluster::getValHitIdx_p1(std::vector<int> hitLine, std::vector<ContigCluster> ccList, int qIdx){
	const int hitLine_size = hitLine.size();
	//const int numThreads = omp_get_num_threads();
	const int numThreads = 1;


	int* hitLineFlags = new int[hitLine_size];
	memset(hitLineFlags,0,hitLine_size*sizeof(int));

	
	int current_idx = 0;
	int t_len = ccList.at(0).seq.size()*2;

	// last element test

	// SSW
	StripedSmithWaterman::Aligner sw(10,9,15,6);
	StripedSmithWaterman::Alignment al;
	StripedSmithWaterman::Filter filter;

	std::cout << ".";

	std::string aligned1, aligned2;
	aligned1.reserve(t_len);
	aligned2.reserve(t_len);

	sw.Align(ccList[hitLine[hitLine_size-1]].seq.c_str(), ccList[qIdx].seq.c_str(), ccList[qIdx].seq.size(), filter, &al, aligned1, aligned2);
	int sim = Tools::calculate(aligned1, aligned2);

	if(sim >= sim_cutoff) hitLineFlags[hitLine_size-1] = 1;
	else hitLineFlags[hitLine_size-1] = -1;
	


	// Needle
	/*
	char* aligned1 = new char[t_len];
	char* aligned2 = new char[t_len];
	
	Needle nwdist(10,0.1);
	nwdist.Align(ccList[qIdx].seq.c_str(), ccList[hitLine[hitLine_size-1]].seq.c_str(), aligned1, aligned2);
	int sim = Tools::calculate(aligned1, aligned2);

	if(sim >= sim_cutoff) hitLineFlags[hitLine_size-1] = 1;
	else hitLineFlags[hitLine_size-1] = -1;
		
	delete [] aligned1;
	delete [] aligned2;
	*/


	if(hitLineFlags[hitLine_size-1] == 1)	return current_idx = hitLine_size-1;

	// remain elements test
	for(int i=hitLine_size-2; i>=0;){
		int termFlags = false;
		int j;
		for(j=i; j>=0 && j-i<numThreads; j--){
			
			// SSW
			
			StripedSmithWaterman::Aligner sw(10,9,15,6);
			StripedSmithWaterman::Alignment al;
			StripedSmithWaterman::Filter filter;

			std::string aligned1, aligned2;
			aligned1.reserve(t_len);
			aligned2.reserve(t_len);

			sw.Align(ccList.at(hitLine[j]).seq.c_str(), ccList.at(qIdx).seq.c_str(), ccList.at(qIdx).seq.size(), filter, &al, aligned1, aligned2);
			int sim = Tools::calculate(aligned1,aligned2);

			if(sim >= sim_cutoff) hitLineFlags[j] = 1;
			else hitLineFlags[j] = -1;
			
			
			
			// Needle
			/*
			char* aligned1 = new char[t_len];
			char* aligned2 = new char[t_len];

			Needle nwdist(10,0.1);
			nwdist.Align(ccList.at(qIdx).seq.c_str(), ccList.at(hitLine[j]).seq.c_str(), aligned1, aligned2);
			int sim = Tools::calculate(aligned1, aligned2);

			if(sim >= sim_cutoff) hitLineFlags[j] = 1;
			else hitLineFlags[j] = -1;
		
			delete [] aligned1;
			delete [] aligned2;
			*/
		}
		
		for(int k=i; k>=j; k--,i--){
			if(hitLineFlags[k] == 1 && hitLineFlags[k+1] == -1){
				current_idx = k;
				termFlags = true;
				break;
			}
		}

		if(termFlags) break;
	}

	delete [] hitLineFlags;

	return current_idx;
}





std::vector<int> KmerCluster::getValHitIdxAllUpdate(std::vector<int> hitLine, std::vector<ContigCluster> ccList, int qIdx){
	int hitLine_size = hitLine.size();

	std::vector<int> valHitLine;
	double sim = 0.0;

	int t_len = ccList.at(0).seq.size()*2;

	// SSW
	std::string aligned1,aligned2;
	aligned1.reserve(t_len);
	aligned2.reserve(t_len);

	StripedSmithWaterman::Aligner sw(10,9,15,6);
	StripedSmithWaterman::Alignment al;
	StripedSmithWaterman::Filter filter;
	std::string sumseq = ccList.at(qIdx).seq;

	for(int i=0; i<hitLine.size(); i++){
		
		sw.Align(ccList.at(hitLine.at(i)).seq.c_str(), sumseq.c_str(), sumseq.size(), filter, &al, aligned1, aligned2);
		sim = Tools::calculate(aligned1, aligned2);

		// flags, left, right update
		if(sim >= sim_cutoff){
				valHitLine.push_back(i);
				sumseq = AssembleGlobal::SumTwoSeq(aligned2,aligned1);
		}
		else continue;
	}

	return valHitLine;
}




