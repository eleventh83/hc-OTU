/***************************************************************************************
ESPRIT: Estimating Species Richness Using Large Collections of 16S rRNA Shotgun Sequences  
Yijun Sun, Yunpeng Cai, Li Liu, Fahong Yu, Michael L. Farrell, William McKendree, William Farmerie
Nucleic Acids Research 2009

Programed by Yunpeng Cai
Copyright Hold by University of Florida, Gainesville, FL32610

****************************************************************************************/
/********************************General Routines************************/


#ifndef UTIL_H
#define UTIL_H

//#ifdef __cplusplus
//extern "C" {
//#endif

#define Sqr(x) ((x)*(x))

void* Malloc(size_t ssize);

inline int Min(int a,int b)
{
	return (a<b)?a:b;
}

inline int Max(int a,int b)
{
	return (a>b)?a:b;
}

inline double Min(double a,double b)
{
	return (a<b)?a:b;	
}

inline double Max(double a,double b)
{
	return (a>b)?a:b;	
}

//#ifdef __cplusplus
//}
//#endif

void AddSuff(char *output,char *input, const char *suff);
void RepExt(char *output, const char *ext);
void GetPath(char *output,char *input);

#define Max_Buf 65535
#define MAXPATH 2000
#define MAX_SEQS 1000000

#endif
