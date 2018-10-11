#include <string>
#include <vector>

#include "FastaSeq.h"
#include "Contig.h"
#include "ContigCluster.h"


ContigCluster::ContigCluster(){
	seq = "";
	noClones = 0;
	solo = false;
}


ContigCluster::~ContigCluster(){};

//int ContigCluster::compareTo(ContigCluster cc){};
