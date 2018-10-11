#include <string>
#include <vector>

#include "Seqtools.h"
#include "generalFunctions.h"

Seqtools::Seqtools(){};
Seqtools::~Seqtools(){};


double Seqtools::calSimilarity(std::string seq1, std::string seq2){
	double total = 0;
	double match = 0;
	int len = seq1.size();
	if(seq1.size() > seq2.size()) len = seq2.size();

	for(int i=0; i<len; i++){
		if(seq1.at(i) == '.') continue;
		if(seq2.at(i) == '.') continue;
		if(seq1.at(i) == '-') continue;
		if(seq2.at(i) == '-') continue;
		if(isMatchedDegenerate(seq1.at(i), seq2.at(i))>0) match++;
		total++;
	}
	return match/total;
}

int Seqtools::isMatchedDegenerate(char base1, char base2){
	if (base1=='-') return 0;
	if (base2=='-') return 0;
	switch (base1) {
	case 'A': 	
		if (base2=='A') return 2;
		if (base2=='R') return 1;
		if (base2=='M') return 1;
		if (base2=='W') return 1;
		if (base2=='D') return 1;
		if (base2=='H') return 1;
		if (base2=='V') return 1;
		if (base2=='N') return 1;
		return 0;
	case 'C': 
		if (base2=='C') return 2;
		if (base2=='Y') return 1;
		if (base2=='M') return 1;
		if (base2=='S') return 1;
		if (base2=='H') return 1;
		if (base2=='B') return 1;
		if (base2=='V') return 1;
		if (base2=='N') return 1;
		return 0;
	case 'G':
		if (base2=='G') return 2;
		if (base2=='R') return 1;
		if (base2=='S') return 1;
		if (base2=='K') return 1;
		if (base2=='D') return 1;
		if (base2=='B') return 1;
		if (base2=='V') return 1;
		if (base2=='N') return 1;
		return 0;
	case 'T': 	
		if (base2=='T') return 2;
		if (base2=='Y') return 1;
		if (base2=='W') return 1;
		if (base2=='K') return 1;
		if (base2=='D') return 1;
		if (base2=='H') return 1;
		if (base2=='B') return 1;
		if (base2=='N') return 1;
		return 0;
	default :
		break;
	}
	if (base1=='R'){
		if (base2=='Y') return 0;
		if (base2=='C') return 0;
		if (base2=='T') return 0;
		return 1;
	}
	if (base1=='Y'){
		if (base2=='R') return 0;
		if (base2=='A') return 0;
		if (base2=='G') return 0;
		return 1;
	}
	if (base1=='M'){
		if (base2=='K') return 0;
		if (base2=='T') return 0;
		if (base2=='G') return 0;				
		return 1;
	}
	if (base1=='K'){
		if (base2=='A') return 0;
		if (base2=='C') return 0;
		if (base2=='M') return 0;
		return 1;
	}
	if (base1=='S'){
		if (base2=='A') return 0;
		if (base2=='T') return 0;
		if (base2=='W') return 0;
		return 1;
	}
	if (base1=='W'){
		if (base2=='S') return 0;
		if (base2=='C') return 0;
		if (base2=='G') return 0;
		return 1;
	}
	if (base1=='V'){
		if (base2=='T') return 0;
		return 1;
	}
	if (base1=='B'){
		if (base2=='A') return 0;
		return 1;
	}
	if (base1=='H'){
		if (base2=='G') return 0;
		return 1;
	}
	if (base1=='D'){
		if (base2=='C') return 0;
		return 1;
	}
	if (base1=='N')	return 1;
	
	return 0;
}



std::string Seqtools::condense(std::string seq){

	char* s = new char[seq.size()];
	int index = 0;
	char char1 = seq[0];
	s[index++] = char1;
	
	for(int i=1; i<(signed)seq.size(); i++){
		if(char1 == seq[i])	continue;
		char1 = seq[i];
		s[index++] = char1;
	}


	std::string r_str = s;
	delete [] s;
	return r_str.substr(0,index);
}

std::string Seqtools::unAlign(std::string seq){
	std::string str;
	seq = toUpperCase(seq);

	for(int i=0; i<(signed)seq.size(); i++){
		switch (seq[i])
		{
		case 'A' : str.append(1,seq[i]); break;
		case 'C' : str.append(1,seq[i]); break;
		case 'G' : str.append(1,seq[i]); break;
		case 'T' : str.append(1,seq[i]); break;
		case 'B' : str.append(1,seq[i]); break;
		case 'V' : str.append(1,seq[i]); break;
		case 'D' : str.append(1,seq[i]); break;
		case 'H' : str.append(1,seq[i]); break;
		case 'K' : str.append(1,seq[i]); break;
		case 'M' : str.append(1,seq[i]); break;
		case 'S' : str.append(1,seq[i]); break;
		case 'W' : str.append(1,seq[i]); break;
		case 'Y' : str.append(1,seq[i]); break;
		case 'R' : str.append(1,seq[i]); break;
		case 'N' : str.append(1,seq[i]); break;
		default: break;
		} 			
	}
	return str;
}

std::string Seqtools::padEndGaps(std::string seq, char fill){
	for(int i=0; i<seq.size();i++){
		if(seq[i]!='-' && seq[i]!='.')	break;
		seq[i] = fill;
	}
	for(int i=seq.size()-1; i>-1; i--){
		if(seq[i]!='-' && seq[i]!='.')	break;
		seq[i]=fill;
	}
	return seq;
}




