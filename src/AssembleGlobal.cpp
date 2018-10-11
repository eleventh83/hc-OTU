#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>


#include "Nucleotide.h"
#include "Seqtools.h"
#include "FastaSeq.h"
#include "Contig.h"
#include "Tools.h"
#include "needle.h"
#include "ssw.h"
#include "ssw_cpp.h"
#include "AssembleGlobal.h"


bool compareSeqLen(std::string a, std::string b){	return a.size()>b.size();}


AssembleGlobal::AssembleGlobal(int initialCapacity){}

AssembleGlobal::AssembleGlobal(std::vector<std::string> inSeq){
	seq = inSeq;
}

AssembleGlobal::~AssembleGlobal(){}

std::string AssembleGlobal::getContig(){
	return contig;
}

std::string AssembleGlobal::getContigWithoutGap(){
	return Seqtools::unAlign(contig);
}

void AssembleGlobal::addSeq(std::string oneSeq){
	seq.push_back(oneSeq);
}


std::string AssembleGlobal::SumTwoSeq(std::string seq1, std::string seq2){
	std::string res(seq1.size(),0);
//#pragma omp parallel num_threads(8)
//#pragma omp parallel for
	for(int i=0; i<(signed)seq1.size(); i++){
		if(seq1.at(i) == seq2.at(i)){
			res.at(i) = seq1.at(i);
			continue;
		}
		if(seq1.at(i)=='-'){
			res.at(i) = seq2.at(i);
			continue;
		}
		if(seq2.at(i)=='-'){
			res.at(i) = seq1.at(i);
			continue;
		}
		res.at(i) = seq1.at(i);
	}
	return res;
}

void AssembleGlobal::doAssemble(){
		

	StripedSmithWaterman::Aligner sw(10,9,15,6);
	StripedSmithWaterman::Alignment al;
	StripedSmithWaterman::Filter filter;

	std::string sumSeq = seq.at(0);
	int t_len = sumSeq.size() * 2;

	std::string aligned1, aligned2;
	aligned1.reserve(t_len);
	aligned2.reserve(t_len);

	// sorting
	std::sort(seq.begin(), seq.end(),compareSeqLen);

	for(int i=1; i<(signed)seq.size(); i++){
		sw.Align(seq.at(i).c_str(), sumSeq.c_str(), sumSeq.size(), filter, &al, aligned1, aligned2);
		//std::cout << i << "-th " << Tools::calculate(aligned1,aligned2) << std::endl;
		if( Tools::calculate(aligned1, aligned2) >= 0.997 ){
				sumSeq = SumTwoSeq(aligned1,aligned2);
		}else{
				Outlier.push_back(i);
		}
	}
	
	contig = sumSeq;


	/*
	for(int i=0; i<(signed)seq.size(); i++){
		sw.Align(seq.at(i).c_str(), contig.c_str(), contig.size(), filter, &al, aligned1, aligned2);
		seq[i] = aligned2;
		//std::cout << "b" ;
	}

	doConsensus();	

	*/

	// Needleman
	/*
	Needle nw(10,1);
	
	char* aligned1 = new char[t_len];
	char* aligned2 = new char[t_len];

	for(int i=1; i<(signed)seq.size(); i++){

		if(seq.size() > 1800 && i == 1001){
			std::cout << "sumSeq	" << sumSeq << "\n";
			std::cout << "seqi	" << seq.at(i) << "\n";
		}

		nw.Align(sumSeq.c_str(), seq.at(i).c_str(), aligned1, aligned2);
		sumSeq = SumTwoSeq(aligned1,aligned2);

		if(seq.size() > 1800 && i==1001){
			std::cout << "align1	" << aligned1 << "\n";
			std::cout << "align2	" << aligned2 << "\n";
			std::cout << "sumSeq2	" << sumSeq << "\n";
		}


	}

	delete [] aligned1;
	delete [] aligned2;

	contig = sumSeq;


	for(int i=0; i<(signed)seq.size(); i++){
		char* naligned1 = new char[t_len];
		char* naligned2 = new char[t_len];

		nw.Align(contig.c_str(), seq.at(i).c_str(),naligned1, naligned2);
		seq[i] = naligned2;

		delete [] naligned1;
		delete [] naligned2;
	}

	doConsensus();
	*/

}


void AssembleGlobal::doConsensus(){
	for(int i=0; i<(signed)seq.size(); i++)	seq.at(i) = Seqtools::padEndGaps(seq.at(i),'~');
	nt.resize(contig.size());
	//for(int n=0; n<(signed)contig.size(); n++)	nt.push_back(Nucleotide());

	for(int i=0; i<(signed)seq.size(); i++){
		for(int n=0; n<(signed)contig.size(); n++){
			nt[n].total++;
			switch(seq.at(i)[n]){
				case 'A' : nt[n].A++; break;
				case 'C' : nt[n].C++; break;
				case 'G' : nt[n].G++; break;
				case 'T' : nt[n].T++; break;
				case '-' : nt[n].Gap++; break;
				case '~' : nt[n].endGap++; break;  
				default: nt[n].N++; break;
			}
		}
	}

	for(int i=0; i<contig.size(); i++)
		contig[i] = nt[i].getMajority_pyro();
}

void AssembleGlobal::print(){
	std::cout << ">contig\n" << contig << "\n";
	for(int i=0; i<(signed)seq.size(); i++)	std::cout << ">" << i << "\n" << seq.at(i) << "\n";
	for(int i=0; i<(signed)contig.size(); i++)
		std::cout << "pos " << i << " A:" << nt.at(i).A << " C:" << nt.at(i).C << " G:" << nt.at(i).G << " T:" << nt.at(i).T << " N:" << nt.at(i).N << " -:" << nt.at(i).Gap;
}
