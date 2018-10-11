#include <iostream>
#include <stdlib.h>
#include <string>

#include "Defs.h"
#include "TBCclass.h"

using namespace std;

int main(int argc, char* argv[]){

	// usage
	if(argc < 4){
		cout << printUsage;
		return 0;
	}

	// fasta file
	string infile = argv[1];

	// threshold
	double sim_thold = 0.97;
	sim_thold = atof(argv[2]);

	// threads
	int num_threads = 1;
	num_threads = atoi(argv[3]);

	// clustering
	TBCclass tbc;
	if(argc == 4)	tbc.assembleWithRead(infile,sim_thold,num_threads,0.06);
	else if(argc == 5){
			tbc.assembleWithRead(infile,sim_thold,num_threads,atof(argv[4]));
	}

	return 1;
}

