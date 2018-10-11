#include "Nucleotide.h"

Nucleotide::Nucleotide(){
	A=0;
	C=0;
	G=0;
	T=0;
	Gap=0;
	endGap=0;
	N=0;
	total=0;
}

Nucleotide::~Nucleotide(){}

char Nucleotide::getMajority_pyro(){
	int MinOccurence = 3;

	double sum=A+C+G+T+Gap;
	double AA=A/sum;
	double CC=C/sum;
	double GG=G/sum;
	double TT=T/sum;
	double vGap=Gap/sum;


	if(sum<10) MinOccurence=1;

	if(AA>=0.5) if  (A>=MinOccurence) return 'A'; else return '-';
	if(CC>=0.5) if  (C>=MinOccurence) return 'C'; else return '-';
	if(GG>=0.5) if  (G>=MinOccurence) return 'G'; else return '-';
	if(TT>=0.5) if  (T>=MinOccurence) return 'T'; else return '-';
	if(vGap>0.5 && Gap>=MinOccurence) return '-';

	return 'N';
}

char Nucleotide::getMajority_normal(){
	double sum=A+C+G+T+Gap;
	double sum2=sum+endGap;

	double AA=A/sum;
	double CC=C/sum;
	double GG=G/sum;
	double TT=T/sum;
	double vGap=Gap/sum;

	if (AA>0.5) if (A==1) return 'a'; else return 'A';
	if (CC>0.5) if (C==1) return 'c'; else return 'C';
	if (GG>0.5) if (G==1) return 'g'; else return 'G';
	if (TT>0.5) if (T==1) return 't'; else return 'T';
	if (vGap>0.5) return '-';

	return 'N';
}