#ifndef TKMERCLUSTER_H_
#define TKMERCLUSTER_H_

class TKmerCluster:KmerCluster{
private:
	std::string infile;
	std::string outfile;
	int nSingleton;
	int nOtus;

public:
	TKmerCluster();
	TKmerCluster(double cut);
	TKmerCluster(std::string cl_id, double sim_cutoff);
	virtual ~TKmerCluster();

	std::string tmpHeader;
	std::string tmpTbc;
	std::string myoutfile;
	std::string otu_rabund;
	std::string idxfile;
	//std::string kmercntfile;

	//PrintWriter writer; // not done!!

	int getSingleton();
	int getOtus();

	//OTUList getOTUList();	// not done!!
	void clustering(std::vector<Contig> ctgList, std::string cl_id, std::vector<int> otuCounts);
	void clustering(std::vector<Contig> ctgList, std::string cl_id);

};

#endif /* TKMERCLUSTER_H_ */
