/***************************************************************************************
ESPRIT: Estimating Species Richness Using Large Collections of 16S rRNA Shotgun Sequences  
Yijun Sun, Yunpeng Cai, Li Liu, Fahong Yu, Michael L. Farrell, William McKendree, William Farmerie
Nucleic Acids Research 2009

Programed by Yunpeng Cai
Copyright Hold by University of Florida, Gainesville, FL32610

****************************************************************************************/
/************************Auxiliary Functions***************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void* Malloc(size_t size)
{
	void *ptr;
	if ((ptr=malloc(size))==NULL)
		{
			fprintf(stderr,"Out of Memory!\n");
			exit(0);
		}
  return ptr;
}

void AddSuff(char *output,char *input, const char *suff)
{
	char *pch,*qch;
	int i;
	qch=input;
	while (*qch=='.' || *qch=='/' || *qch=='\\')
	    qch++;
	qch=strchr(qch,'.');

	qch=strchr(input,'.');
	if (qch ==NULL)
		{
			sprintf(output,"%s%s",input,suff);
			return;
		}
	do	
	{
		pch=qch;
		qch=strchr(pch+1,'.');
	}while(qch !=NULL);
	
	for (i=0;input+i < pch;i++)
		output[i]=input[i];
	output[i]=0;
	strcat(output,suff);
	strcat(output,pch);
}

void RepExt(char *output, const char *ext)
{
	char *pch,*qch;
	int i;
	qch=output;
	while (*qch=='.' || *qch=='/' || *qch=='\\')
	    qch++;
	qch=strchr(qch,'.');
	if (qch ==NULL)
		{
			strcat(output,".");
			strcat(output,ext);
			return;
		}
	do	
	{
		pch=qch;
		qch=strchr(pch+1,'.');
	}while(qch !=NULL);
	
	pch++;
	strcpy(pch,ext);
}

void GetPath(char *output,char *input)
{
	char *pch,*qch;
	if (strchr(input,'/')==NULL && strchr(input,'\\')==NULL)
	{
		*output=0;
		return;
	}
	pch=input+strlen(input);
	while (pch >input && *pch !='/' && *pch !='\\') 
	{pch--;}
	for (qch=input;qch <pch;qch++)
		*(output++)=*qch;
	*output=0;		
}
