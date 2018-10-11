#ifndef CONTIGCLUSTER_H_
#define CONTIGCLUSTER_H_

class ContigCluster{
public:
	ContigCluster();
	~ContigCluster();

	std::string seq;
	int noClones;
	bool solo;

	std::vector<Contig> contigList;
	std::vector<int> ContigIdxList;
	std::vector<int> kmercount;

	//int compareTo(ContigCluster cc)

};


#endif /* CONTIGCLUSTER_H_ */
