#ifndef CLTAG_H_
#define CLTAG_H_

namespace CLTag{
	int FILE_TYPE_FASTA = 0;
	int FILE_TYPE_FASTQ = 1;
	int FILE_TYPE_CLF = 2; 

	/* tag */
	std::string COMMENT		= "//";
	std::string RETRACTOR		= "\\";
	std::string DELIMITER		= "|";
	std::string DBINFO			= "/!gene";
	std::string COUNT			= "/!reads";

	std::string TOTAL_RAW_READS_COUNT     = "/!total_raw_reads";
	std::string UNSORTED_RAW_READS_COUNT     = "/!unsorted_raw_reads";
	std::string SORTED_RAW_READS_COUNT     = "/!sorted_raw_reads";
	std::string DROPPED_RAW_READS_COUNT     = "/!dropped_raw_reads";

	std::string MAX_LENGTH		= "/!maxlen";
	std::string MIN_LENGTH		= "/!minlen";
	std::string AVG_LENGTH		= "/!avglen";
	std::string CLUSTER		= "/!clusters";
	std::string START			= "/!begin";
	std::string END			= "/!end";
	std::string BARCODE		= "/!barcode";
	std::string FLINKER		= "/!linker_for";
	std::string RLINKER		= "/!linker_rev";
	std::string FPRIMER		= "/!primer_for";
	std::string RPRIMER		= "/!primer_rev";
	std::string TYPE			= "/!type";
	std::string METHOD			= "/!method";
	std::string CONTIGS		= "/!contigs";
	std::string MACHINE		= "/!machine";
	std::string CUTOFF			= "/!cutoff";

	std::string FPRIMER_SEQ = "GAGTTTGATCMTGGCTCAG";
	std::string RPRIMER_SEQ = "WTTACCGCGGCTGCTGG";

	/* tag value */
	std::string TYPE_SEQ		= "name_seq";
	std::string TYPE_SEQ_QUAL	= "name_seq_qual";
	std::string TYPE_CLF		= "clf";
	std::string TYPE_CM		= "cm";
	std::string TYPE_CONTIG	= "contig";
	std::string TYPE_FASTA	= "fasta";
	std::string TYPE_CLUSTER	= "cluster";
	std::string TITANIUM		= "454-TITANIUM";
	std::string FLX			= "454-FLX";
	std::string GS20			= "GS20";
	std::string GA2X    = "GA2X";
	std::string HISEQ    = "HISEQ";
	std::string MISEQ    = "MISEQ";


	std::string FASTQ_INIT		= "@";
	std::string FASTA_INIT		= ">";
	std::string CONTIG_INIT	= "%";
	std::string COMMUNITY_INIT	= "CN ";
	std::string COMMUNITY_CNT	= "ME ";
	std::string COMMUNITY_SEQ	= "SQ ";
	std::string COMMUNITY_ID	= "ID ";

	std::string NEW_LINE = "\n";
	std::string TAB = "\t";

	std::string EXTENSION			= ".clf";
	std::string CM_EXTENSION		= ".cm";
	std::string CONTIG_EXTENSION	= ".contig";
	std::string SORT_SUFFIX		= ".sort";
	std::string REPOSITORY			= "../repository/";

	std::string LENGTH_STAT		="/!lengthstat";

	/* cdhit values */
	std::string CDHIT_SINGLETON	= "/!cdhit_singleton";
	std::string CDHIT_OTUS			= "/!cdhit_otus";
	std::string CDHIT_RF			= "/!cdhit_rf";

	std::string CDHIT_CHAO			= "/!cdhit_chao";
	std::string CDHIT_CHAO_LCI		= "/!cdhit_chao_lci";
	std::string CDHIT_CHAO_HCI		= "/!cdhit_chao_hci";

	std::string CDHIT_ACE			= "/!cdhit_ace";
	std::string CDHIT_ACE_LCI		= "/!cdhit_ace_lci";
	std::string CDHIT_ACE_HCI		= "/!cdhit_ace_hci";

	std::string CDHIT_JACKKNIFE			= "/!cdhit_jackknife";
	std::string CDHIT_JACKKNIFE_LCI		= "/!cdhit_jackknife_lci";
	std::string CDHIT_JACKKNIFE_HCI		= "/!cdhit_jackknife_hci";

	std::string CDHIT_SHANNON				= "/!cdhit_shannon";
	std::string CDHIT_SHANNON_LCI			= "/!cdhit_shannon_lci";
	std::string CDHIT_SHANNON_HCI			= "/!cdhit_shannon_hci";
	std::string CDHIT_NPSHANNON			= "/!cdhit_npshannon";


	std::string CDHIT_SIMPSON				= "/!cdhit_simpson";
	std::string CDHIT_SIMPSON_LCI			= "/!cdhit_simpson_lci";
	std::string CDHIT_SIMPSON_HCI			= "/!cdhit_simpson_hci";


	/* cdhit value for target taxon */

	std::string READS_TARGET = "/!reads_target";
	std::string CDHIT_SINGLETON_TARGET = "/!cdhit_singleton_target";
	std::string CDHIT_OTUS_TARGET     = "/!cdhit_otus_target";
	std::string CDHIT_RF_TARGET     = "/!cdhit_rf_target";

	std::string CDHIT_CHAO_TARGET     = "/!cdhit_chao_target";
	std::string CDHIT_CHAO_LCI_TARGET   = "/!cdhit_chao_lci_target";
	std::string CDHIT_CHAO_HCI_TARGET   = "/!cdhit_chao_hci_target";

	std::string CDHIT_ACE_TARGET      = "/!cdhit_ace_target";
	std::string CDHIT_ACE_LCI_TARGET    = "/!cdhit_ace_lci_target";
	std::string CDHIT_ACE_HCI_TARGET    = "/!cdhit_ace_hci_target";

	std::string CDHIT_JACKKNIFE_TARGET      = "/!cdhit_jackknife_target";
	std::string CDHIT_JACKKNIFE_LCI_TARGET    = "/!cdhit_jackknife_lci_target";
	std::string CDHIT_JACKKNIFE_HCI_TARGET    = "/!cdhit_jackknife_hci_target";

	std::string CDHIT_SHANNON_TARGET        = "/!cdhit_shannon_target";
	std::string CDHIT_SHANNON_LCI_TARGET      = "/!cdhit_shannon_lci_target";
	std::string CDHIT_SHANNON_HCI_TARGET      = "/!cdhit_shannon_hci_target";
	std::string CDHIT_NPSHANNON_TARGET      = "/!cdhit_npshannon_target";


	std::string CDHIT_SIMPSON_TARGET        = "/!cdhit_simpson_target";
	std::string CDHIT_SIMPSON_LCI_TARGET      = "/!cdhit_simpson_lci_target";
	std::string CDHIT_SIMPSON_HCI_TARGET      = "/!cdhit_simpson_hci_target";

	std::string CDHIT_GOOD_LIBRARY ="/!cdhit_good_library";
	std::string CDHIT_GOOD_LIBRARY_TARGET ="/!cdhit_good_library_target";

	/* end of cdhit */

	/* blast values */
	std::string BLAST_OTUS					= "/!blast_otus";
	std::string BLAST_RF					= "/!blast_rf";
	std::string BLAST_SINGLETON = "/!blast_singleton";

	std::string BLAST_CHAO					= "/!blast_chao";
	std::string BLAST_CHAO_LCI				= "/!blast_chao_lci";
	std::string BLAST_CHAO_HCI				= "/!blast_chao_hci";

	std::string BLAST_ACE					= "/!blast_ace";
	std::string BLAST_ACE_LCI				= "/!blast_ace_lci";
	std::string BLAST_ACE_HCI				= "/!blast_ace_hci";

	std::string BLAST_JACKKNIFE			= "/!blast_jackknife";
	std::string BLAST_JACKKNIFE_LCI		= "/!blast_jackknife_lci";
	std::string BLAST_JACKKNIFE_HCI		= "/!blast_jackknife_hci";

	std::string BLAST_SHANNON				= "/!blast_shannon";
	std::string BLAST_SHANNON_LCI			= "/!blast_shannon_lci";
	std::string BLAST_SHANNON_HCI			= "/!blast_shannon_hci";
	std::string BLAST_NPSHANNON			= "/!blast_npshannon";

	std::string BLAST_SIMPSON				= "/!blast_simpson";
	std::string BLAST_SIMPSON_LCI			= "/!blast_simpson_lci";
	std::string BLAST_SIMPSON_HCI			= "/!blast_simpson_hci";

	std::string BLAST_GOOD_LIBRARY ="/!cdhit_good_library";


	std::string BLAST_OTUS_TARGET          = "/!blast_otus_target";
	std::string BLAST_RF_TARGET         = "/!blast_rf_target";
	std::string BLAST_SINGLETON_TARGET = "/!blast_singleton_target";

	std::string BLAST_CHAO_TARGET         = "/!blast_chao_target";
	std::string BLAST_CHAO_LCI_TARGET       = "/!blast_chao_lci_target";
	std::string BLAST_CHAO_HCI_TARGET       = "/!blast_chao_hci_target";

	std::string BLAST_ACE_TARGET          = "/!blast_ace_target";
	std::string BLAST_ACE_LCI_TARGET        = "/!blast_ace_lci_target";
	std::string BLAST_ACE_HCI_TARGET        = "/!blast_ace_hci_target";

	std::string BLAST_JACKKNIFE_TARGET      = "/!blast_jackknife_target";
	std::string BLAST_JACKKNIFE_LCI_TARGET    = "/!blast_jackknife_lci_target";
	std::string BLAST_JACKKNIFE_HCI_TARGET    = "/!blast_jackknife_hci_target";

	std::string BLAST_SHANNON_TARGET        = "/!blast_shannon_target";
	std::string BLAST_SHANNON_LCI_TARGET      = "/!blast_shannon_lci_target";
	std::string BLAST_SHANNON_HCI_TARGET      = "/!blast_shannon_hci_target";
	std::string BLAST_NPSHANNON_TARGET      = "/!blast_npshannon_target";

	std::string BLAST_SIMPSON_TARGET        = "/!blast_simpson_target";
	std::string BLAST_SIMPSON_LCI_TARGET      = "/!blast_simpson_lci_target";
	std::string BLAST_SIMPSON_HCI_TARGET      = "/!blast_simpson_hci_target";

	std::string BLAST_GOOD_LIBRARY_TARGET ="/!cdhit_good_library_target";

	/* end of blast */

	std::string DATA_VERSION  = "/!data_version"; // extaxon db version


	std::string USER_UID = "/!user_uid";
	std::string USER_GROUP_UID = "/!user_group_uid";
	std::string SAMPLING_DATE = "/!sampling_date";
	std::string USER_SAMPLE_NAME = "/!user_sample_name";
	std::string SAMPLE_DESCRIPTION = "/!sample_description";
	std::string COUNTRY = "/!country";
	std::string COORDS_NS = "/!coord_ns";
	std::string COORDS_EW = "/!coord_ew";
	std::string N_CHIMERA= "/!n_chimera"; // number of chimera
	std::string COMM_VER = "/!comm_ver"; // data version( table structure )
	std::string PHYLUM_COMP = "/!phylum_comp"; // phylum composition
	std::string SOURCE_KEYWORD = "/!source_keyword"; // keywor for sample
	std::string CHUNLAB_UID = "/!chunlab_uid"; // chunlab uid begin with C_


	//std::string SEPERATOR = System.getProperty("file.separator");

	//std::string os = System.getProperty("os.name");

	std::string DATA_DIR = "/jc/cldb/run/";
	std::string OUTPUT_DIR = "/jc/cldb/comm/";
	std::string EUK_DBNAME  = "/jc/db/lsu_d1d2.fasta";
	std::string DBNAME	= "/jc/db/ref16s_prok.fasta";
	std::string SMALLDBNAME = "/jc/db/filter_16s_prok.fasta";
	std::string EZDBNAME	= "/jc/db/ref16s_prok.fasta";
	std::string CHIMERA_FREE_DB = "/jc/db/ref16s_non_chimera.fasta";

	std::string TMP_DIR = "/tmp";
	std::string TMP_DIR_WIN = "C:\\tmp";

}


#endif /* CLTAG_H_ */
