#ifndef KMERCLUSTER_H_
#define KMERCLUSTER_H_

class KmerCluster{
public:
	KmerCluster();
	KmerCluster(double sim_cutoff);
	virtual ~KmerCluster();


	//double stemp;	// not done!!
	//std::string blastFileName;	// not done!!
	std::string temp_dir;
	double sim_cutoff;
	double sim_cutoff_margin;
	double solo_cutoff;
	std::vector<ContigCluster> ccList;
	int totalClones;
	int totalOTUs;
	int blast_type;	// not done!!
	std::string blast_path;	// not done!!


	//int reset_blastdb(int pos);	// not done!!
	void runContig(std::vector<Contig> ctgList, double sim_cutoff, double sim_cutoff_margin, std::string temp_dir, int blast_type, std::string blast_path); // not done!!
	void runContig(std::vector<Contig> ctgList);
	void runContig_rev1(std::vector<Contig> ctgList);
	void runContig_p1(std::vector<Contig> ctgList);
	void runContig_all(std::vector<Contig> ctgList);
	//OTUList getOTUList();	// not done!!
	int getValHitIdx(std::vector<int> hitLine, std::vector<ContigCluster> ccList, int qIdx);
	std::vector<int> getValHitIdxAll(std::vector<int> hitLine, std::vector<ContigCluster> ccList, int qIdx);
	std::vector<int> getValHitIdxAllUpdate(std::vector<int> hitLine, std::vector<ContigCluster> ccList, int qId);
	int getValHitIdx_p1(std::vector<int> hitLine, std::vector<ContigCluster> ccList, int qIdx);

};

#endif /* KMERCLUSTER_H_ */
