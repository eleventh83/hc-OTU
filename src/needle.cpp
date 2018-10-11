/***************************************************************************************
ESPRIT: Estimating Species Richness Using Large Collections of 16S rRNA Shotgun Sequences  
Yijun Sun, Yunpeng Cai, Li Liu, Fahong Yu, Michael L. Farrell, William McKendree, William Farmerie
Nucleic Acids Research 2009

Programed by Yunpeng Cai
Copyright Hold by University of Florida, Gainesville, FL32610

****************************************************************************************/
/************************Needleman Wunsch Alignment***************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "needle.h"
#include "util.h"

#define MinFloat -1e100

DIRS VTracks::GetMaxInc(double ia, double ib, double ic, double &newval)
{
	DIRS dirs=DIR_UP;
	newval=vals[DIR_UP]+ic;
	if (vals[DIR_LEFT] + ib > newval)
		{
			newval= vals[DIR_LEFT] + ib;
			dirs=DIR_LEFT;
		}	
	if (vals[DIR_DIAG] + ia > newval)
		{
			newval= vals[DIR_DIAG] + ia;
			dirs=DIR_DIAG;
		}	
	return dirs;
}

Needle::Needle(double go, double gl)
{
		gap_open=-1*go;
		gap_len=-1*gl;
}

double Needle::Score(int row, int col){
    if (row==4 && col==4) return 0;
    if (row==col) return 10;
    return -9;
}

int* Needle::StrToCode(const char* str)
{
	int len=strlen(str);
	int *l=(int *)Malloc(strlen(str)*sizeof(int));
	
	for (int i = 0; i < len; i++)
	{
		switch(toupper(str[i]))
		{
		case 'A':
				l[i]=0;
				break;
		case 'G':
				l[i]=1;
				break;
		case 'C':
				l[i]=2;
				break;
		case 'T':
				l[i]=3;
				break;
		case 'N':
			  	l[i]=4;		
				break;
		default:
				l[i]=-1;
				break;
		}
	}
	return l;
}

double Needle::CalculateMatrix(int* source, int* dest, int slen, int dlen,DIRTracks *pr)
{
  int i,j,t;
  double inc,kdist;
 // FILE *fp;
  
	VTracks *ar = new VTracks[dlen+1];
  	VTracks lefttrack,thistrack;
	DIRS dirs;
 
 	//fp=fopen("nul","w");
 	
 	ar[0].SetVal(DIR_DIAG,0);
 	ar[0].SetVal(DIR_LEFT,MinFloat);
 	ar[0].SetVal(DIR_UP,MinFloat);
  
  for (i = 1; i <= dlen; i++)
  {
  	ar[i].SetVal(DIR_DIAG,MinFloat);
  	ar[i].SetVal(DIR_UP,MinFloat);
    ar[i].SetVal(DIR_LEFT,gap_open + gap_len* (i -1));
  }
        
        
  for (i = 1; i <= slen; i++)
  {
  	lefttrack.SetVal(DIR_DIAG,MinFloat);
  	lefttrack.SetVal(DIR_UP,gap_open + gap_len* (i -1));
  	lefttrack.SetVal(DIR_LEFT,MinFloat);
  	
    for (j = 1; j <= dlen; j++)
      {  
  			 dirs=ar[j-1].GetMaxInc(0,0,0,kdist);
         kdist+=Score(source[i-1],dest[j-1]);                 
				 pr[i*(dlen+1)+j].SetVal(DIR_DIAG,dirs);
				 thistrack.SetVal(DIR_DIAG,kdist);
  			 dirs=lefttrack.GetMaxInc(gap_open,gap_len,gap_open,kdist);
				 pr[i*(dlen+1)+j].SetVal(DIR_LEFT,dirs);
				 thistrack.SetVal(DIR_LEFT,kdist);
  			 dirs=ar[j].GetMaxInc(gap_open,gap_open,gap_len,kdist);
				 pr[i*(dlen+1)+j].SetVal(DIR_UP,dirs);
				 thistrack.SetVal(DIR_UP,kdist);
				 ar[j-1]=lefttrack;
				 lefttrack=thistrack;
      }
    ar[dlen]=thistrack;  

//	  for (t=1;t<=dlen;t++)
 //	 {
 // 		fprintf(fp,"(%5.1f%5.1f%5.1f)",ar[t].GetVal(DIR_DIAG),ar[t].GetVal(DIR_LEFT),ar[t].GetVal(DIR_UP)); 
 // 	}      
//	  for (t=1;t<=dlen;t++)
//	 	{
//  		fprintf(fp,"(%5d%5d%5d)",pr[i][t].GetVal(DIR_DIAG),pr[i][t].GetVal(DIR_LEFT),pr[i][t].GetVal(DIR_UP)); 
 // 	}      
//  	fprintf(fp,"\n");
  }
  
  dirs=ar[dlen].GetMaxInc(0,0,0,kdist);
  pr[slen*(dlen+1)+dlen].SetEntry(dirs);
//	fprintf(fp,"%d %f\n",dirs,kdist);
//  fclose(fp);

  delete [] ar;
  return kdist;
}
  
void Needle::GetAlignments(DIRTracks *pr, int dlen, const char *sA, const char *sB, char *alA,char *alB)
{
  char *buf,*ptrA,*ptrB,*pbuf;
	int i,j;
	DIRS endpt,newpt;
	
  ptrA=alA;
  ptrB=alB;       
         
  i=strlen(sA);
  j=strlen(sB);          

	buf=(char *)Malloc((i+j+10)*sizeof(char));

  while (i > 0 && j > 0)
  {
    endpt=pr[i*(dlen+1)+j].GetEntry();
    newpt=pr[i*(dlen+1)+j].GetEntryVal();
    if (endpt==DIR_DIAG)
       {
         *(ptrA++) = sA[i-1];
         *(ptrB++) = sB[j-1];
         i--;j--;                
       }
    else if (endpt==DIR_UP)
      {
        *(ptrA++) = sA[i-1];
        *(ptrB++) = '-';
        i--;
      }
    else if (endpt==DIR_LEFT)
      {
        *(ptrA++) = '-';
        *(ptrB++) = sB[j-1];
        j--;
      }
    else
    {	
    	fprintf(stderr,"Error No Match %d %d\n",i,j);
    	exit(0);
    }  
    pr[i*(dlen+1)+j].SetEntry(newpt);
  }
  while(i > 0)
  {
    *(ptrA++) = sA[i-1];
    *(ptrB++) = '-';
    i--;            
  }
  while(j > 0)
  {
    *(ptrA++) = '-';
    *(ptrB++) = sB[j - 1];
    j--;            
  }
	ptrA--;ptrB--;

	pbuf=buf;
	do
	{
		*(pbuf++)=*(ptrA--);
	}while (ptrA >=alA);
	*(pbuf)=0;
	strcpy(alA,buf);
	  
	pbuf=buf;
	do
	{
		*(pbuf++)=*(ptrB--);
	}while (ptrB >=alB);
	*(pbuf)=0;
	strcpy(alB,buf);
	free (buf);
}

void Needle::Align(const char *seq1,const char *seq2,char *al1,char *al2)
{
	DIRTracks *pr;
	double score;
	int len1,len2;
	int *code1, *code2;
	int i,j;
	
	len1=strlen(seq1);
	len2=strlen(seq2);
	code1=StrToCode(seq1);
	code2=StrToCode(seq2);
	
	pr=(DIRTracks *)Malloc((len1+1)*(len2+1)*sizeof(DIRTracks));
	
  //for (i=0;i<=len1;i++)
  //{
 	//	pr[i]=new DIRTracks[len2+1];
	//} 		

	score=CalculateMatrix(code1, code2, len1, len2, pr);
	GetAlignments(pr, len2, seq1, seq2, al1, al2);

/*	for (i=0;i<=len1;i++)
	{
		delete [] pr[i];
	}*/
	
	free(pr);
	free(code1);
	free(code2);
}
