#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdlib>



#include "generalFunctions.h"
#include "FastaSeq.h"
#include "Contig.h"
#include "Tools.h"

double Tools::calculate(std::string s1, std::string s2){
	int gap = 0;
	int match = 0;
	int mismatch = 0;
	int endgap = 0;
	int insertion = 0;
	int deletion = 0;
	int sum_match_mismatch = 0;


	if (s1.size() != s2.size()) return -1;
	if ( s1[0]=='-' || s1[s1.size()-1] == '-') s1 = Tools::padEndGaps(s1, '.');
	if ( s2[0]=='-' || s2[s2.size()-1] == '-') s2 = Tools::padEndGaps(s2, '.');

	gap = match = mismatch = endgap = 0;

	for (int i=0; i<s1.size(); i++){
		if (s1[i]=='-' && s2[i]=='-') continue;
		if (s1[i]=='.' && s2[i]=='.') continue;
		if (s1[i]=='.' && s2[i]=='-') continue;
		if (s1[i]=='-' && s2[i]=='.') continue;
		if (s1[i]=='.') { endgap++; continue; }
		if (s2[i]=='.') { endgap++; continue; }
		if (s1[i]=='-') { gap++; insertion++; continue; }
		if (s2[i]=='-') { gap++;  deletion++; continue; }
		if (s1[i]=='A'){
			if (s2[i]=='A') { match++; continue; }
			if (s2[i]=='C') { mismatch++; continue; }
			if (s2[i]=='G') { mismatch++; continue; }
			if (s2[i]=='T') { mismatch++; continue; }
		}
		if (s1[i]=='C'){
			if (s2[i]=='A') { mismatch++; continue; }
			if (s2[i]=='C') { match++; continue; }
			if (s2[i]=='G') { mismatch++; continue; }
			if (s2[i]=='T') { mismatch++; continue; }
		}
		if (s1[i]=='G'){
			if (s2[i]=='A') { mismatch++; continue; }
			if (s2[i]=='C') { mismatch++; continue; }
			if (s2[i]=='G') { match++; continue; }
			if (s2[i]=='T') { mismatch++; continue; }
		}
		if (s1[i]=='T'){
			if (s2[i]=='A') { mismatch++; continue; }
			if (s2[i]=='C') { mismatch++; continue; }
			if (s2[i]=='G') { mismatch++; continue; }
			if (s2[i]=='T') { match++; continue; }
		}
	}
	sum_match_mismatch = match+mismatch;
	return ( (double)match / ( (double) mismatch + (double) match));
	//return ( (double)match / (double)s1.size() );
}

std::string Tools::padEndGaps(std::string seq, char fill){
	for(int i=0; i<seq.size(); i++){
		if (seq[i]!='-') break;
		seq[i]=fill;
	}
	for (int i=seq.size()-1; i>-1; i--){
		if (seq[i]!='-') break;
		seq[i]=fill;  
	}
	return seq;
}

std::string Tools::padString(std::string s, int n, char c, bool paddingLeft){
	std::string pad(n,c);
	//std::string pad;
	//for(int i=0; i<n; i++) pad[i] = c;
	
	if(paddingLeft)	s = pad + s;
	else	s = s + pad;
	
	return s;
}

double Tools::getTimeStampNow(){
	return 0.0;
}	//not done!!


std::vector<Contig> Tools::ContigFile2ContigList(std::string contig_fn){
	std::vector<Contig> sList;

	std::ifstream inFile(contig_fn.c_str());
	std::string line;

	while(!inFile.eof()){
		
		std::getline(inFile,line);
		if(line.substr(0,5) == "/!end")	break;
		if(line.substr(0,2) == "/!")	continue;
		if(line.at(0) == '%'){
			Contig ctg;
			std::vector<std::string> s = StringSplit(line,"\\|",4);
			if(s.size() == 3){
				ctg.title = s.at(0).substr(1);
				ctg.noClones = atoi(s.at(2).c_str());
			}else{
				ctg.title = line.substr(1);
				ctg.noClones = 1;
			}


			std::string tline;
			std::getline(inFile,tline);
			s = StringSplit(tline,"\t",3);
			ctg.seq = s.at(0);
			if(ctg.noClones == 1) ctg.idxList.push_back(atoi(s.at(1).c_str()));

			std::getline(inFile,tline);
			if(ctg.noClones > 1){
				tline.clear();
				std::getline(inFile,tline);

				while(1){
					line.clear();
					std::getline(inFile,line);
					if(line == "") break;

					if(line.substr(0,2) == "//" )	break;
					s = StringSplit(line,"\t",4);
					FastaSeq fa;
					fa.title = s.at(0);
					fa.sequence = s.at(1);
					ctg.seqList.push_back(fa);
					ctg.idxList.push_back(atoi(s.at(2).c_str()));
				}
				if(ctg.noClones != ctg.seqList.size()){
					std::cout << "Illegal Contig file. Invalid clone numbers\n";
				}
			}
			else if(ctg.noClones == 1){
				FastaSeq fa;
				fa.title = ctg.title;
				fa.sequence = ctg.seq;	// hmmmmmmmmmmmmmmmmmmmm
				ctg.seqList.push_back(fa);
			}
			sList.push_back(ctg);
		}
	}

	inFile.close();

	return sList;
}


