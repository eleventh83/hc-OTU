/***************************************************************************************
ESPRIT: Estimating Species Richness Using Large Collections of 16S rRNA Shotgun Sequences  
Yijun Sun, Yunpeng Cai, Li Liu, Fahong Yu, Michael L. Farrell, William McKendree, William Farmerie
Nucleic Acids Research 2009

Programed by Yunpeng Cai
Copyright Hold by University of Florida, Gainesville, FL32610

****************************************************************************************/
/*******************************Needleman Wunsch Algorithm******************************/


#ifndef Needle_H
#define Needle_H

enum DIRS{DIR_DIAG=0,DIR_LEFT,DIR_UP};

class DIRTracks{
 protected:
 	DIRS defVal;
 	DIRS dirs[3];
 
 public:
 	DIRTracks(){defVal=DIR_DIAG;};
 	~DIRTracks(){};
 	DIRS GetVal(DIRS id) { return dirs[id];};
	void SetVal(DIRS id,DIRS val) { dirs[id]=val;};
	DIRS GetEntry(){return defVal;};
	DIRS GetEntryVal(){return dirs[defVal];};
	void SetEntry(DIRS id){defVal=id;};
};

class VTracks{
 protected:
 	double vals[3];
 
 public:
 	VTracks(){};
 	~VTracks(){};
 	double GetVal(DIRS id) { return vals[id];};
	void SetVal(DIRS id,double val) { vals[id]=val;};
	
	VTracks & operator = (const VTracks  & v)
	{
		for (int i=0;i<3;i++) 
		  vals[i]=v.vals[i];
		return *this;
	}
	DIRS GetMaxInc(double ia, double ib, double ic, double &newval);
};


class Needle{
	private:
	double gap_open;
	double gap_len; 
	double score;
	
	void GetAlignments(DIRTracks *pr, int dlen,const char *sA, const char *sB, char *alA,char *alB);
	double CalculateMatrix(int* source, int* dest, int slen, int dlen, DIRTracks *pr);
	
	protected: //overridable components
	virtual int* StrToCode(const char* str); 
	// Index characters into a code, the index table must be zero-based
	virtual double Score(int row, int col); 
	// Provide the award/penalty for one pair of indexed characters,
	// A positive score favorates the match while a negative score discourage it                                   

	public:
		Needle(double go, double gl); //set the penalty for starting and extending a gap
		~Needle(){};
		void Align(const char *seq1, const char *seq2,char *al1, char *al2); 
		//seq1, seq2: input sequences; al1, al2: Output aligned sequances
};

#endif
