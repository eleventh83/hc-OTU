#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <omp.h>


#include "generalFunctions.h"
#include "FastaSeq.h"
#include "Contig.h"
#include "Tools.h"
#include "ContigCluster.h"
#include "KmerCluster.h"
#include "CLHeader.h"
#include "TKmerCluster.h"


TKmerCluster::TKmerCluster(){
	infile = "";
	outfile = "";
	nSingleton = 0;
	nOtus = 0;
	
	tmpHeader = "";
	tmpTbc = "";
	myoutfile = "";
	otu_rabund = "";
	idxfile = "";
}


TKmerCluster::TKmerCluster(double cut):KmerCluster(cut){
	infile = "";
	outfile = "";
	nSingleton = 0;
	nOtus = 0;
	
	tmpHeader = "";
	tmpTbc = "";
	myoutfile = "";
	otu_rabund = "";
	idxfile = "";
}


// this one!
TKmerCluster::TKmerCluster(std::string cl_id, double sim_cutoff):KmerCluster(sim_cutoff){
	infile = cl_id + ".0.contig";
	outfile = cl_id + ".tbc";
	tmpHeader = cl_id + ".theader.tmp";
	tmpTbc = cl_id + ".ttmp";
	myoutfile = cl_id + ".tbc";
	otu_rabund = cl_id + ".trabund";
	idxfile = cl_id + ".idx.txt";
	//kmercntfile = cl_id + ".0.kmercount";

	std::ofstream writer;

	std::cout << "\nProcessing " + infile;
	std::vector<Contig> ctgList = Tools::ContigFile2ContigList(infile);


	std::cout << "\nctgList : " + int2str(ctgList.size()) + "\n";

	std::cout << "Run Contig Process\n";

	// merging contigs
	runContig_all(ctgList);

	// assign singletons
	

	std::cout << "End Contig Process\n";

	
	// idxfile write
	writer.open(idxfile.c_str());
	int* clusteridx = new int[totalClones];
	int noCluster = 1;
	for(int i=0; i<ccList.size(); i++){
		if(ccList[i].noClones == 0) continue;
		for(int j=0; j<ccList[i].noClones; j++){
			clusteridx[ccList[i].ContigIdxList[j]] = noCluster;
		}
		noCluster++;
	}
	for(int i=0; i<totalClones ; i++) writer << int2str(clusteridx[i]) + "\n";
	writer.flush();
	writer.close();
	delete [] clusteridx;

	
	// trabund file write
	writer.open(otu_rabund.c_str());
	writer << dbl2str(sim_cutoff) + "	" + int2str(totalOTUs) + "	";
	for(int i=0; i<ccList.size(); i++){
		if(ccList[i].noClones > 0){
			writer << int2str(ccList[i].noClones) + "	";
			if(ccList[i].noClones == 1) nSingleton++;
		}
		else break;
	}
	writer.flush();
	writer.close();

	// tbc file write
	writer.open(outfile.c_str());
	
	writer << "/!reads|" + int2str(totalClones) + "\n";
	writer << "/!maxlen|\n";
	writer << "/!minlen|\n";
	writer << "/!avglen|\n";
	writer << "/!clusters|" + int2str(ctgList.size()) + "\n";
	writer << "/!type|cluster\n";
	writer << "/!method|kmer\n";
	writer << "/!cutoff|" + dbl2str(sim_cutoff) + "\n";
	writer << "/!comm_ver|0.0\n";	


	for(int i=0; i<ccList.size(); i++){
		if(ccList[i].noClones < 1) continue;
		
		int totalClone = 0;
		std::string clusterRep;
		bool repSeleted = false;
		std::string byte0s;

		for(int j=0; j<ccList[i].contigList.size(); j++){
			if(!repSeleted){
				clusterRep = ccList[i].contigList[j].title;
				repSeleted = true;
			}

			totalClone += ccList[i].contigList[j].seqList.size();

			for(int k=0; k<ccList[i].contigList[j].seqList.size(); k++)
				byte0s += ccList[i].contigList[j].seqList[k].title + "\t" + ccList[i].contigList[j].seqList[k].sequence + "\n";
		}

		writer << "%" + clusterRep + "|ME|" + int2str(totalClone) + "\n";
		writer << ccList[i].seq + "\n";

		if(totalClone > 1){
			writer << "//\nBegin\n";
			writer << byte0s;
		}
		writer << "//\n";
	}
	writer << "/!end\n";
	writer.flush();
	writer.close();

}

TKmerCluster::~TKmerCluster(){}



int TKmerCluster::getSingleton(){
	return nSingleton;
}

int TKmerCluster::getOtus(){
	return nOtus;
}

//OTUList TKmerCluster::getOTUList(){};	// not done!!

void TKmerCluster::clustering(std::vector<Contig> ctgList, std::string cl_id){

	otu_rabund = cl_id + ".rabund";
	outfile = cl_id + ".otu";

	std::ofstream writer(otu_rabund.c_str());

	runContig(ctgList);
	std::cout << "\nOTUS = " + int2str(totalOTUs) + "\nClones = " + int2str(totalClones);

	nOtus = totalOTUs;

	writer << dbl2str(sim_cutoff) + "	" + int2str(totalOTUs) + "	";
	for(int i=0; i<ccList.size(); i++){
		if(ccList.at(i).noClones > 0){
			writer << int2str(ccList.at(i).noClones) + "	";
			if(ccList.at(i).noClones == 1) nSingleton++;
		}else
			break;
	}

	writer.close();

	writer.open(outfile.c_str());
	for(int i=0; i<ccList.size(); i++){
		ContigCluster cc = ccList.at(i);
		if(cc.noClones < 1)		continue;
		Contig ctg = cc.contigList.at(0);
		writer << ctg.title;
		for(int j=1; j<cc.contigList.size(); j++){
			ctg = cc.contigList.at(j);
			writer << "\t" + ctg.title;
		}
		writer << "\n";
		writer.flush();
	}
	writer.close();
}

void TKmerCluster::clustering(std::vector<Contig> ctgList, std::string cl_id,std::vector<int> otuCounts){
	otu_rabund = cl_id + ".rabund";
	outfile = cl_id + ".otu";

	std::vector<int> merged;

	std::ofstream writer(otu_rabund.c_str());
	runContig(ctgList);
	for(int i=0; i<ccList.size(); i++){
		if(ccList.at(i).noClones > 0)	merged.push_back(ccList.at(i).noClones);
		else	break;
	}

	merged.insert(merged.end(), otuCounts.begin(), otuCounts.end());
	std::sort(merged.begin(), merged.end());
	
	writer << dbl2str(sim_cutoff) + "	" + int2str(totalOTUs + otuCounts.size()) + "	";
	for(int i=merged.size()-1; i>=0; i--){
		writer << int2str(merged.at(i)) + "	";
		if(merged.at(i) == 1)	nSingleton++;
	}
	writer.flush();
	writer.close();
}









