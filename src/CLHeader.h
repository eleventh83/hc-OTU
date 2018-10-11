#ifndef CLHEADER_H_
#define CLHEADER_H_

class CLHeader{
private:
	int nTotalRawRead;
	int nUnsortedRawRead;
	int nSortedRawRead;
	int nDroppedRawRead;
	int nContigs;
	int nClusters;
	std::string strFormat;
	std::string strMethod;
	std::string strCutoff;
	int user_uid;  

public:
	CLHeader();
	~CLHeader();


	std::string meta;
	std::string db_info;
	int reads;

	int max_len;
	int min_len;
	double avg_len;

	std::string barcode;
	std::string primer_for;
	std::string linker_for;
	std::string primer_rev;
	std::string linker_rev;
	std::string machine;


	std::string length_stat;

	int cdhit_otus;
	int cdhit_singleton;
	std::string cdhit_rf;
	double cdhit_chao;
	double cdhit_chao_lci;
	double cdhit_chao_hci;
	double cdhit_ace;
	double cdhit_ace_lci;
	double cdhit_ace_hci;
	double cdhit_jackknife;
	double cdhit_jackknife_lci;
	double cdhit_jackknife_hci;
	double cdhit_shannon;
	double cdhit_shannon_lci;
	double cdhit_shannon_hci;
	double cdhit_npshannon;
	double cdhit_simpson;
	double cdhit_simpson_lci;
	double cdhit_simpson_hci;

	int reads_target;
	int cdhit_otus_target;
	int cdhit_singleton_target;
	std::string cdhit_rf_target;
	double cdhit_chao_target;
	double cdhit_chao_lci_target;
	double cdhit_chao_hci_target;
	double cdhit_ace_target;
	double cdhit_ace_lci_target;
	double cdhit_ace_hci_target;
	double cdhit_jackknife_target;
	double cdhit_jackknife_lci_target;
	double cdhit_jackknife_hci_target;
	double cdhit_shannon_target;
	double cdhit_shannon_lci_target;
	double cdhit_shannon_hci_target;
	double cdhit_npshannon_target;
	double cdhit_simpson_target;
	double cdhit_simpson_lci_target;
	double cdhit_simpson_hci_target;

	double cdhit_good_library;
	double cdhit_good_library_target;

	int blast_otus;
	std::string blast_rf;
	int blast_singleton;
	double blast_chao;
	double blast_chao_lci;
	double blast_chao_hci;
	double blast_ace;
	double blast_ace_lci;
	double blast_ace_hci;
	double blast_jackknife;
	double blast_jackknife_lci;
	double blast_jackknife_hci;
	double blast_shannon;
	double blast_shannon_lci;
	double blast_shannon_hci;
	double blast_npshannon;
	double blast_simpson;
	double blast_simpson_lci;
	double blast_simpson_hci;

	int blast_otus_target;
	std::string blast_rf_target;
	int blast_singleton_target;
	double blast_chao_target;
	double blast_chao_lci_target;
	double blast_chao_hci_target;
	double blast_ace_target;
	double blast_ace_lci_target;
	double blast_ace_hci_target;
	double blast_jackknife_target;
	double blast_jackknife_lci_target;
	double blast_jackknife_hci_target;
	double blast_shannon_target;
	double blast_shannon_lci_target;
	double blast_shannon_hci_target;
	double blast_npshannon_target;
	double blast_simpson_target;
	double blast_simpson_lci_target;
	double blast_simpson_hci_target;

	double blast_good_library;
	double blast_good_library_target;

	int user_group_uid;
	std::string user_sample_name;
	std::string sample_description;
	std::string sampling_date;
	int data_version;
	int n_chimera;
	std::string coord_ns;
	std::string coord_ew;
	double comm_ver;
	std::string phylum_comp;

	std::string source_keyword;

	std::string chunlab_uid;

	/////
	void setInfo(std::string line);
	void setInfo(std::string tag, std::string value);
	std::string getInfo(std::string tag);	
	//void writeToWriter(std::ostream writer);
	//void writeToByteArray(ByteArrayOutputStream byteArray);
	//	ByteArrayOutputStream writeToByteArray();
	std::string writeToByteArray();

	void setReadCount( int nTmp );
	void setMaxLen( int nTmp ); 
	void setMinLen( int nTmp );	
	void setAvgLen( double dTmp );
	void setContigCount( int nTmp ); 

	int getDroppedRawReads(); 
	int getUnsortedRawReads(); 
	int getSortedRawReads(); 
	int getTotalRawReads(); 

	//public ArrayList<Field> getPublicFields()
};

#endif /* CLHEADER_H_ */
