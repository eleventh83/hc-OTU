#include "KmerHitLine.h"

KmerHitLine::KmerHitLine(){
	contigidx = -1;
	kmerscore = 1;
}

KmerHitLine::KmerHitLine(int idx, double score){
	contigidx = idx;
	kmerscore = score;
}

KmerHitLine::~KmerHitLine(){}

void KmerHitLine::setIdx(int idx){
	contigidx = idx;
}

void KmerHitLine::setScore(double score){
	kmerscore = score;
}


//int kmerHitLine::compareTo(kmerHitLine o){};	// not done !!
