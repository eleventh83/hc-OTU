#include <iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <climits>
#include <math.h>
#include <sys/time.h>
#include <omp.h>


#include "generalFunctions.h"
#include "Fasta.h"
#include "Fastq.h"
#include "CLSQSequence.h"
#include "CLSQFastaReader.h"
#include "AlignNoGap.h"
#include "AlignNoGapWrap.h" // for parallel 
#include "Nucleotide.h"
#include "AssembleGlobal.h"
#include "CLHeader.h"
#include "FastaSeq.h"
#include "Contig.h"
#include "Tools.h"
#include "ContigCluster.h"
#include "KmerDist.h"
#include "KmerHitLine.h"
#include "needle.h"
#include "KmerCluster.h"
#include "TKmerCluster.h"
#include "TBCclass.h"




bool compareCLSQ(CLSQSequence a, CLSQSequence b){	return a.clonesFor.size()>b.clonesFor.size();	}
bool compareIdx(int a, int b) { return a<b; }

std::vector<double> elapsed_time;
double BEGIN, END;



static inline double myclock(){
	struct timeval tval;
	gettimeofday(&tval, NULL);
	return (tval.tv_sec *1000 + tval.tv_usec / 1000);
}


TBCclass::TBCclass(){
	nSingleton = 0;
	nOtus = 0;
	inmc = "";
	nIter = 1000;
}

TBCclass::~TBCclass(){}

void TBCclass::assembleWithRead(std::string inFile_, double similarity_, int numThreads, double kmerthreshold ){

	// threads set
	omp_set_num_threads(numThreads);

	std::string inFile = inFile_;
	double sim_thold = similarity_;
	std::string outfile = "";
	int noMismatchAllowed = 0; 
	int min_overlap = 200;
	double min_overlap_rate = 0.8; //temporal
	int noSeqRead = 0;

	// set file name
	if(!inFile.empty()){
		inmc = replace(".fasta","",inFile);
		//std::cout << inmc.c_str() << "\n";
		//std::cout << inFile.c_str() << "\n";
		
		std::stringstream t_noMismatchAllowed;
		t_noMismatchAllowed << noMismatchAllowed;
		outfile = inmc + "." + t_noMismatchAllowed.str() + ".contig";
		std::cout << "outfile : " << outfile.c_str() << "\n";

	}

	std::vector<CLSQSequence> seqList;
	std::vector<CLSQSequence> oriSeqList;

	
	BEGIN = myclock();
	// Step 1 - Start	-------------------------------------------------------------------
	CLSQFastaReader sqFastaReader = CLSQFastaReader(inFile);
	sqFastaReader.read();
	seqList = sqFastaReader.getList();
	oriSeqList = seqList;
	noSeqRead = seqList.size();
	for(int i=0; i<(signed)seqList.size(); i++)	seqList.at(i).clonesFor.push_back(i);
	//clusterFoStream = new FileOutputStream(outfile); // contig file out !!! - not Done
	// Step 1 - End	-----------------------------------------------------------------------
	printf("Step 1 End\n");
	END = myclock();
	elapsed_time.push_back((END-BEGIN)/1000);


	// calculate the mean length of condensed sequences
	int totallen = 0;
	for(int i=0; i<(signed)seqList.size(); i++) totallen += seqList.at(i).getCondensedSeq().size();
	min_overlap = (int)(totallen/seqList.size()*min_overlap_rate);


	BEGIN = myclock();
	// Step 2 - Start	-------------------------------------------------------------------
	int tmp_index = 0;
		
	
	// original
	
	/*
	int num_op = 0;

	AlignNoGap ang;
	ang.setMinimumOverlap(min_overlap);

	for(int i=0; i<(signed)seqList.size(); i++){
		//if(tmp_index++ == 100) {
		//	tmp_index = 0;
		//	std::cout << ".";
		//}
		if(seqList[i].getCloneNum() == 0)	continue;

		for(int j=i+1; j<(signed)seqList.size(); j++){
			if(seqList[j].getCloneNum() == 0)	continue;
			ang.setSequences(seqList[i].getCondensedSeq(), seqList[j].getCondensedSeq());
			num_op++;


			if(ang.doAlign_mistch_allow(0) != -1){
				seqList[i].addClone();
				seqList[i].addCloneFor(j);
				seqList[i].setCondensedSeq(ang.getSumSeq());
				seqList[j].setCloneNum(0);
			}
		}
	}
	printf("# operation is %d\n",num_op);
	*/



	/*
	// parallel - new 1
	
	int chunk_size = 80;
	std::vector<int> angFlags(chunk_size,-1); // -2:skip -1:not aligned o/w:aligned 
	std::vector<AlignNoGap> angList(chunk_size);
	for(int i=0; i<chunk_size; i++) angList[i].setMinimumOverlap(min_overlap);

	
	for(int i=0; i<(signed)seqList.size(); i++){
		if(seqList[i].getCloneNum() == 0) continue;
		
		int jobSize = chunk_size;
		
		for(int j=i+1; j<(signed)seqList.size();){

			//std::cout << "i " << i << " j " << j << "\n";
			
			if(j+jobSize-1 >= (signed)seqList.size()) jobSize = seqList.size()-j; 
			
			for(int idx =0; idx <jobSize; idx++) angFlags[idx] = -1; // reset flags
			
			#pragma omp parallel for schedule(auto)
			for(int k=0; k<jobSize ; k++){
				if(seqList[k+j].getCloneNum() == 0) continue;
				
				angList[k].setSequences(seqList[i].getCondensedSeq(), seqList[k+j].getCondensedSeq());
				//angList[k].setSequences(seqList[i].getSequence() , seqList[k+j].getSequence());

				angFlags[k] = angList[k].doAlign_mistch_allow(0);
			}

			for(int l=0; l<jobSize; l++,j++){
				if(angFlags[l] > -1){
					seqList[i].addClone();
					seqList[i].addCloneFor(j);
					seqList[i].setCondensedSeq(angList[l].getSumSeqMod2());
					seqList[j].setCloneNum(0);

					break;
				}
			}

		}
	}
	*/

/*
	// parallel - another (faster)
	for(int i=0; i<(signed)seqList.size(); i++){
			if(seqList[i].getCloneNum() == 0) continue;

			std::vector<int> toAddIdx;

			#pragma omp parallel for schedule(dynamic)
			for(int j=i+1; j<(signed)seqList.size(); j++){
				if(seqList[j].getCloneNum() == 0)	continue;

				AlignNoGap ang(seqList[i].getCondensedSeq(), seqList[j].getCondensedSeq(), min_overlap);
				if(ang.doAlign_mistch_allow(0) > -1){
					#pragma omp critical (CTG_UPDATE)
						{
							//seqList[i].addClone();
							//seqList[i].addCloneFor(j);
							toAddIdx.push_back(j);
							seqList[j].setCloneNum(0);
						}
				}
			}

			// update (forward)
			std::sort(toAddIdx.begin(), toAddIdx.end(), compareIdx);
			for(int k=0; k<toAddIdx.size(); k++){
					seqList[i].addClone();
					seqList[i].addCloneFor(toAddIdx[k]);
			}
	}

*/
	// parallel - kmer
	for(int i=0; i<(signed)seqList.size(); i++){
			if(seqList[i].getCloneNum() == 0) continue;

			std::vector<int> toAddIdx;

//			KmerDist tmp_kmer(5);
//			tmp_kmer.setFreqA(seqList[i].getCondensedSeq());

			#pragma omp parallel for schedule(dynamic)
			for(int j=i+1; j<(signed)seqList.size();j++){
					if(seqList[j].getCloneNum() == 0) continue;
					
					KmerDist tmp_kmer(5);
					tmp_kmer.setFreqA(seqList[i].getCondensedSeq());
					tmp_kmer.setFreqB(seqList[j].getCondensedSeq());
					if(tmp_kmer.getsim() < kmerthreshold){
						#pragma omp critical (CTG_UPDATE)
							{
									toAddIdx.push_back(j);
									seqList[j].setCloneNum(0);
							}
					}
			}

			std::sort(toAddIdx.begin(), toAddIdx.end(), compareIdx);
			for(int k=0; k<toAddIdx.size(); k++){
					seqList[i].addClone();
					seqList[i].addCloneFor(toAddIdx[k]);
			}
	}

	




	// Step 2 - End	-----------------------------------------------------------------------
	printf("Step 2 End\n");
	END = myclock();
	elapsed_time.push_back((END-BEGIN)/1000);





	BEGIN = myclock();
	// Step 3 - Start	-------------------------------------------------------------------
	/*
	tmp_index = 0;
	#pragma omp parallel for schedule(auto)
	for(int i=0; i<(signed)seqList.size(); i++){
		//if(tmp_index++ == 9){
		//	tmp_index = 0;
		//	std::cout << "\n" << (i+1) ;
		//}

		if( seqList[i].getCloneNum() == 0 )	continue;
		if( seqList[i].getCloneNum() == 1 ){
			seqList[i].setContig(seqList[i].getSequence());
			//seqList[i].setKmercount(KmerDist::getFreqFromNT(seqList[i].contig,5));

			continue;
		}
		
		AssembleGlobal ass(seqList.at(i).getCloneNum());
		for(int j=0; j<(signed)seqList.at(i).clonesFor.size(); j++){

			int pos = seqList.at(i).clonesFor.at(j);
			ass.addSeq(seqList.at(pos).getSequence());
		}
		ass.doAssemble();

		seqList.at(i).setContig(ass.getContigWithoutGap());
		//seqList.at(i).setKmercount(KmerDist::getFreqFromNT(seqList.at(i).contig,5));

	}
	*/
	// Step 3 - End	-----------------------------------------------------------------------




	// Step 3 (revised) - Start	-------------------------------------------------------------------
	tmp_index = 0;
	std::vector<int> noneSingleton;
	for(int i=0; i<(signed)seqList.size(); i++){
		//if(tmp_index++ == 9){
		//	tmp_index = 0;
		//	std::cout << "\n" << (i+1) ;
		//}

		if( seqList[i].getCloneNum() == 0 )	continue;
		if( seqList[i].getCloneNum() == 1 ){
				seqList[i].setContig(seqList[i].getSequence());
				continue;
		}
		noneSingleton.push_back(i);
	}
		
	#pragma omp parallel for schedule(dynamic)
	for(int i=0; i<(signed)noneSingleton.size(); i++){
		AssembleGlobal ass(seqList.at(noneSingleton[i]).getCloneNum());
		for(int j=0; j<(signed)seqList.at(noneSingleton[i]).clonesFor.size(); j++){
			int pos = seqList.at(noneSingleton[i]).clonesFor.at(j);
			ass.addSeq(seqList.at(pos).getSequence());
		}
		ass.doAssemble();

		seqList.at(noneSingleton[i]).setContig(ass.getContigWithoutGap());
		
		//std::cout << ass.Outlier.size() << std::endl;
		/*
		for(int j=ass.Outlier.size(); j>0; j--){
				int pos = seqList.at(noneSingleton[i]).clonesFor.at(ass.Outlier.at(j-1));
				seqList.at(noneSingleton[i]).clonesFor.erase(seqList.at(noneSingleton[i]).clonesFor.begin()+ass.Outlier.at(j-1));
				seqList.at(pos).setCloneNum(1);
				seqList.at(pos).setContig(seqList[pos].getSequence());
		}
		*/
		


	}
	// Step 3 (revised) - End	-----------------------------------------------------------------------
	
	printf("Step 3 End\n");
	END = myclock();
	elapsed_time.push_back((END-BEGIN)/1000);



	BEGIN = myclock();
	// Step 4 - Start	-------------------------------------------------------------------
	std::sort(seqList.begin(), seqList.end(),compareCLSQ);
	// Step 4 - End	-----------------------------------------------------------------------
	printf("Step 4 End\n");
	END = myclock();
	elapsed_time.push_back((END-BEGIN)/1000);



	// contig file out - Start 	-----------------------------------------------------------
	std::ofstream contigfile;
	contigfile.open(outfile.c_str());

	std::ostringstream ClusterOutStream;

	int noClones = 0;
	for(int i=0; i<(signed)seqList.size(); i++){
		if(seqList.at(i).getCloneNum() == 0) continue;
		
		ClusterOutStream << "%" + seqList.at(i).getTitle() + "|ME|" + int2str(seqList.at(i).getCloneNum()) + "\n";
		ClusterOutStream << seqList.at(i).getContig() + "\t";
		ClusterOutStream << int2str(seqList.at(i).clonesFor.at(0)) + "\n";
		ClusterOutStream << "//" << "\n";

		noClones++;

		if(seqList.at(i).getCloneNum() > 1){
			ClusterOutStream << "Start" << "\n";
			for(int c=0; c<(signed)seqList.at(i).clonesFor.size(); c++){
				ClusterOutStream << oriSeqList.at(seqList.at(i).clonesFor.at(c)).getTitle() + "\t";
				ClusterOutStream << oriSeqList.at(seqList.at(i).clonesFor.at(c)).getSequence() + "\t";
				ClusterOutStream << int2str(seqList.at(i).clonesFor.at(c)) + "\n";
			}
			ClusterOutStream << "//\n";
		}
	}


	//std::ostringstream HeaderStream;
	CLHeader newHeader;
	newHeader.setReadCount(noSeqRead);
	newHeader.setContigCount(noClones);
	int nTmpMax = INT_MIN;
	int nTmpMin = INT_MAX;

	long tmpSum = 0;

	for(int i=0; i< noSeqRead; i++){
		nTmpMax = std::max(nTmpMax, (int)seqList.at(i).getSequence().size());
		nTmpMin = std::min(nTmpMin, (int)seqList.at(i).getSequence().size());
		tmpSum += seqList.at(i).getSequence().size();
	}

	newHeader.setMaxLen(nTmpMax);
	newHeader.setMinLen(nTmpMin);
	newHeader.setAvgLen(tmpSum / (double) (seqList.size() * 100));
	newHeader.setInfo("/!type", "contig");
	
	contigfile << newHeader.writeToByteArray();
	contigfile << ClusterOutStream.str();
	contigfile << "/!end" << std::endl;

	contigfile.close();

	// contif file out - End 	-----------------------------------------------------------

	// KmerCount file out - Start ------------------------------------------------
//	std::ofstream Kmercountfile;
//	Kmercountfile.open(kmerFreqFile.c_str());
//	std::ostringstream KmerFreqOutStream;
//
//	std::ostringstream KmerFreqOutStream;
//
//	Kmercountfile.open(kmerFreqFile.c_str());
//	std::ostringstream KmerFreqOutStream;
//
//			KmerFreqOutStream << "%" + seqList.at(i).getTitle() + "|ME|" + int2str(seqList.at(i).getCloneNum()) + "\n";
//			KmerFreqOutStream << seqList.at(i).kmercount[0];
//			for(int j=1; j<(int)pow((double)4,5); j++){
//				KmerFreqOutStream << " " << + seqList.at(i).kmercount[j];
//			}
//			KmerFreqOutStream << "\n";
//	}
//
//	Kmercountfile << KmerFreqOutStream.str();
//	Kmercountfile << "/!end" << std::endl;
//	Kmercountfile.close();

	// KmerCount file out - End --------------------------------------------------

	// Profile file out - Start --------------------------------------------
	
	// KmerCount file out - End -------------------------------------------------



	//omp_set_num_threads(1);

	BEGIN = myclock();
	// Step 5 - Start	-------------------------------------------------------------------
	TKmerCluster tbc(inmc,sim_thold);
	// Step 5 - End	-----------------------------------------------------------------------
	printf("Step 5 End\n");
	END = myclock();
	elapsed_time.push_back((END-BEGIN)/1000);




	

	// time display
	double tottime = 0;
	std::cout << "\n";
	for(int i=0; i<(signed)elapsed_time.size(); i++){
		std::cout << "Step " << i+1 << " : " << elapsed_time.at(i) << "\n";
		tottime += elapsed_time.at(i);
	}
	std::cout << "Total : " << tottime << "\n";


	return;



}


