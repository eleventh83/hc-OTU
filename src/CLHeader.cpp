#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <climits>

#include "CLTag.h"
#include "CLHeader.h"
#include "generalFunctions.h"


using namespace CLTag;

CLHeader::CLHeader(){
	// private
	nTotalRawRead = -1;
	nUnsortedRawRead = -1;
	nSortedRawRead = -1;
	nDroppedRawRead = -1;
	nContigs = -1;
	nClusters = -1;
	strFormat = "";
	strMethod = "";
	strCutoff = "";
	user_uid = -1;  



	// public
	meta = "";
	db_info = "";
	reads = -1;

	max_len	= INT_MAX;
	min_len	= INT_MIN;
	avg_len	= -1;

	barcode	= "";
	primer_for	= "";
	linker_for	= "";
	primer_rev	= "";
	linker_rev	= "";
	machine	= "";


	length_stat	= "";

	cdhit_otus = -1;
	cdhit_singleton = -1;
	cdhit_rf = "";
	cdhit_chao = -1;
	cdhit_chao_lci =-1;
	cdhit_chao_hci =-1;
	cdhit_ace = -1;
	cdhit_ace_lci = -1;
	cdhit_ace_hci   = -1;
	cdhit_jackknife = -1;
	cdhit_jackknife_lci   = -1;
	cdhit_jackknife_hci   = -1;
	cdhit_shannon   = -1;
	cdhit_shannon_lci     = -1;
	cdhit_shannon_hci     = -1;
	cdhit_npshannon = -1;
	cdhit_simpson   = -1;
	cdhit_simpson_lci     = -1;
	cdhit_simpson_hci = -1;

	reads_target = -1;
	cdhit_otus_target = -1;
	cdhit_singleton_target = -1;
	cdhit_rf_target = "";
	cdhit_chao_target = -1;
	cdhit_chao_lci_target =-1;
	cdhit_chao_hci_target =-1;
	cdhit_ace_target = -1;
	cdhit_ace_lci_target = -1;
	cdhit_ace_hci_target   = -1;
	cdhit_jackknife_target = -1;
	cdhit_jackknife_lci_target   = -1;
	cdhit_jackknife_hci_target   = -1;
	cdhit_shannon_target   = -1;
	cdhit_shannon_lci_target     = -1;
	cdhit_shannon_hci_target     = -1;
	cdhit_npshannon_target = -1;
	cdhit_simpson_target   = -1;
	cdhit_simpson_lci_target     = -1;
	cdhit_simpson_hci_target = -1;

	cdhit_good_library = -1;
	cdhit_good_library_target = -1;

	blast_otus = -1;
	blast_rf = "";
	blast_singleton= -1;
	blast_chao = -1;
	blast_chao_lci =-1;
	blast_chao_hci =-1;
	blast_ace = -1;
	blast_ace_lci = -1;
	blast_ace_hci   = -1;
	blast_jackknife = -1;
	blast_jackknife_lci   = -1;
	blast_jackknife_hci   = -1;
	blast_shannon   = -1;
	blast_shannon_lci     = -1;
	blast_shannon_hci     = -1;
	blast_npshannon = -1;
	blast_simpson   = -1;
	blast_simpson_lci     = -1;
	blast_simpson_hci = -1;

	blast_otus_target = -1;
	blast_rf_target = "";
	blast_singleton_target= -1;
	blast_chao_target = -1;
	blast_chao_lci_target =-1;
	blast_chao_hci_target =-1;
	blast_ace_target = -1;
	blast_ace_lci_target = -1;
	blast_ace_hci_target = -1;
	blast_jackknife_target = -1;
	blast_jackknife_lci_target = -1;
	blast_jackknife_hci_target = -1;
	blast_shannon_target = -1;
	blast_shannon_lci_target = -1;
	blast_shannon_hci_target = -1;
	blast_npshannon_target = -1;
	blast_simpson_target = -1;
	blast_simpson_lci_target = -1;
	blast_simpson_hci_target = -1;;

	blast_good_library = -1;
	blast_good_library_target = -1;

	user_group_uid = -1;
	user_sample_name = "";
	sample_description = "";
	sampling_date = "";
	data_version = -1;
	n_chimera = -1;
	coord_ns = "";
	coord_ew = "";
	comm_ver = 0.0;
	phylum_comp = "";

	source_keyword = "";

	chunlab_uid = "";
}

CLHeader::~CLHeader(){};

void CLHeader::setInfo(std::string line){
	std::vector<std::string> tagInfo;


	if( (line.find(CDHIT_RF) != -1 ) || (line.find(DBINFO) != -1) || (line.find(LENGTH_STAT) != -1) || (line.find(PHYLUM_COMP) != -1) || (line.find(SOURCE_KEYWORD) != -1) ) tagInfo = StringSplit(line, RETRACTOR + DELIMITER, 2);
	else tagInfo = StringSplit(line,RETRACTOR + DELIMITER, 0);

	if(line.substr(0,2) == "//" && tagInfo.size() < 2 ) return;

	if( tagInfo[0] == DBINFO ) db_info = bothtrim(tagInfo[1]);
	else if( tagInfo[0]== COUNT ) reads = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== MAX_LENGTH ) max_len = atoi(bothtrim(tagInfo[1]).c_str());

	else if( tagInfo[0]== MIN_LENGTH ) min_len = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== AVG_LENGTH ) avg_len = atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CONTIGS ) nContigs = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CLUSTER ) nClusters = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BARCODE ) barcode = bothtrim(tagInfo[1]);
	else if( tagInfo[0]== FLINKER ) linker_for =  bothtrim(tagInfo[1]);
	else if( tagInfo[0]== FPRIMER ) primer_for =  bothtrim(tagInfo[1]);
	else if( tagInfo[0]== RLINKER ) linker_rev =  bothtrim(tagInfo[1]);
	else if( tagInfo[0]== RPRIMER ) primer_rev =  bothtrim(tagInfo[1]);
	else if( tagInfo[0]== MACHINE ) machine =  bothtrim(tagInfo[1]);
	else if( tagInfo[0]== TYPE ) 	strFormat =  bothtrim(tagInfo[1]);
	else if( tagInfo[0]== METHOD )	strMethod =  bothtrim(tagInfo[1]);
	else if( tagInfo[0]== CUTOFF )	strCutoff =  bothtrim(tagInfo[1]);

	else if( tagInfo[0]== TOTAL_RAW_READS_COUNT )  nTotalRawRead =  atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== UNSORTED_RAW_READS_COUNT )  nUnsortedRawRead =  atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== SORTED_RAW_READS_COUNT )  nSortedRawRead =  atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== DROPPED_RAW_READS_COUNT )  nDroppedRawRead =  atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SINGLETON )	 	cdhit_singleton = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_OTUS ) 			cdhit_otus = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_RF )				cdhit_rf = bothtrim(tagInfo[1]) ;
	else if( tagInfo[0]== CDHIT_CHAO )			cdhit_chao = atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_CHAO_LCI )		cdhit_chao_lci = atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_CHAO_HCI )		cdhit_chao_hci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_ACE )			cdhit_ace =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_ACE_LCI )		cdhit_ace_lci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_ACE_HCI )		cdhit_ace_hci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_JACKKNIFE )		cdhit_jackknife =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_JACKKNIFE_LCI )	cdhit_jackknife_lci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_JACKKNIFE_HCI )	cdhit_jackknife_hci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SHANNON ) 		cdhit_shannon =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SHANNON_LCI )	cdhit_shannon_lci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SHANNON_HCI ) 	cdhit_shannon_hci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_NPSHANNON ) 		cdhit_npshannon =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SIMPSON ) 		cdhit_simpson =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SIMPSON_LCI ) 	cdhit_simpson_lci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SIMPSON_HCI )	cdhit_simpson_hci =  atof(bothtrim(tagInfo[1]).c_str());

	else if( tagInfo[0]== READS_TARGET )   reads_target = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SINGLETON_TARGET )   cdhit_singleton_target = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_OTUS_TARGET )      cdhit_otus_target = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_RF_TARGET )        cdhit_rf_target = bothtrim(tagInfo[1]) ;
	else if( tagInfo[0]== CDHIT_CHAO_TARGET )      cdhit_chao_target = atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_CHAO_LCI_TARGET )    cdhit_chao_lci_target = atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_CHAO_HCI_TARGET )    cdhit_chao_hci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_ACE_TARGET )     cdhit_ace_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_ACE_LCI_TARGET )   cdhit_ace_lci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_ACE_HCI_TARGET )   cdhit_ace_hci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_JACKKNIFE_TARGET )   cdhit_jackknife_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_JACKKNIFE_LCI_TARGET ) cdhit_jackknife_lci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_JACKKNIFE_HCI_TARGET ) cdhit_jackknife_hci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SHANNON_TARGET )     cdhit_shannon_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SHANNON_LCI_TARGET ) cdhit_shannon_lci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SHANNON_HCI_TARGET )   cdhit_shannon_hci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_NPSHANNON_TARGET )     cdhit_npshannon_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SIMPSON_TARGET )     cdhit_simpson_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SIMPSON_LCI_TARGET )   cdhit_simpson_lci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_SIMPSON_HCI_TARGET ) cdhit_simpson_hci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_GOOD_LIBRARY)   cdhit_good_library=  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== CDHIT_GOOD_LIBRARY_TARGET ) cdhit_good_library_target =  atof(bothtrim(tagInfo[1]).c_str());

	else if( tagInfo[0]== BLAST_OTUS ) 			blast_otus = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_RF )				blast_rf = bothtrim(tagInfo[1]) ;
	else if( tagInfo[0]== BLAST_SINGLETON )        blast_singleton = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_CHAO )			blast_chao = atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_CHAO_LCI )		blast_chao_lci = atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_CHAO_HCI )		blast_chao_hci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_ACE )			blast_ace =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_ACE_LCI )		blast_ace_lci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_ACE_HCI )		blast_ace_hci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_JACKKNIFE )		blast_jackknife =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_JACKKNIFE_LCI )	blast_jackknife_lci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_JACKKNIFE_HCI )	blast_jackknife_hci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SHANNON ) 		blast_shannon =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SHANNON_LCI )	blast_shannon_lci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SHANNON_HCI ) 	blast_shannon_hci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_NPSHANNON ) 	blast_npshannon =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SIMPSON ) 		blast_simpson =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SIMPSON_LCI ) 	blast_simpson_lci =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SIMPSON_HCI )	blast_simpson_hci =  atof(bothtrim(tagInfo[1]).c_str());

	else if( tagInfo[0]== BLAST_OTUS_TARGET )      blast_otus_target = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_RF_TARGET )        blast_rf_target = bothtrim(tagInfo[1]) ;
	else if( tagInfo[0]== BLAST_SINGLETON_TARGET )        blast_singleton_target = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_CHAO_TARGET )      blast_chao_target = atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_CHAO_LCI_TARGET )    blast_chao_lci_target = atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_CHAO_HCI_TARGET )    blast_chao_hci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_ACE_TARGET )     blast_ace_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_ACE_LCI_TARGET )   blast_ace_lci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_ACE_HCI_TARGET )   blast_ace_hci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_JACKKNIFE_TARGET )   blast_jackknife_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_JACKKNIFE_LCI_TARGET ) blast_jackknife_lci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_JACKKNIFE_HCI_TARGET ) blast_jackknife_hci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SHANNON_TARGET )     blast_shannon_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SHANNON_LCI_TARGET ) blast_shannon_lci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SHANNON_HCI_TARGET )   blast_shannon_hci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_NPSHANNON_TARGET )   blast_npshannon_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SIMPSON_TARGET )     blast_simpson_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SIMPSON_LCI_TARGET )   blast_simpson_lci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_SIMPSON_HCI_TARGET ) blast_simpson_hci_target =  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_GOOD_LIBRARY)   blast_good_library=  atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== BLAST_GOOD_LIBRARY_TARGET ) blast_good_library_target =  atof(bothtrim(tagInfo[1]).c_str());

	else if( tagInfo[0]== LENGTH_STAT )			length_stat	 =  bothtrim(tagInfo[1]);
	else if( tagInfo[0]== USER_GROUP_UID ) user_group_uid = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== USER_UID ) user_uid = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== USER_SAMPLE_NAME ) user_sample_name = bothtrim(tagInfo[1]);
	else if( tagInfo[0]== SAMPLE_DESCRIPTION ) sample_description = bothtrim(tagInfo[1]);
	else if( tagInfo[0]== USER_SAMPLE_NAME ) sampling_date = bothtrim(tagInfo[1]);
	else if( tagInfo[0]== COORDS_NS ) coord_ns = bothtrim(tagInfo[1]);
	else if( tagInfo[0]== COORDS_EW ) coord_ew = bothtrim(tagInfo[1]);
	else if( tagInfo[0]== DATA_VERSION ) data_version = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== N_CHIMERA ) n_chimera = atoi(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== COMM_VER ) comm_ver = atof(bothtrim(tagInfo[1]).c_str());
	else if( tagInfo[0]== PHYLUM_COMP ) phylum_comp = bothtrim(tagInfo[1]);
	else if( tagInfo[0]== SOURCE_KEYWORD ) source_keyword = bothtrim(tagInfo[1]);
	else if( tagInfo[0]== CHUNLAB_UID ) chunlab_uid = bothtrim(tagInfo[1]);

	else 
	{
		if( line.substr(0,2) == "//" ) meta += (line + NEW_LINE);
	}
}


void CLHeader::setInfo(std::string tag, std::string value){
	if( value == "" ) return;

	if( tag == DBINFO ) db_info = bothtrim(value);
	else if( tag == COUNT ) reads = atoi(bothtrim(value).c_str());
	else if( tag == MAX_LENGTH ) max_len = atoi(bothtrim(value).c_str());
	else if( tag == MIN_LENGTH ) min_len = atoi(bothtrim(value).c_str());
	else if( tag == AVG_LENGTH ) avg_len = atof(bothtrim(value).c_str());
	else if( tag == CONTIGS )	nContigs = atoi(bothtrim(value).c_str());
	else if( tag == CLUSTER ) nClusters = atoi(bothtrim(value).c_str());
	else if( tag == BARCODE )	barcode = bothtrim(value);
	else if( tag == FLINKER )	linker_for =  bothtrim(value);
	else if( tag == FPRIMER )	primer_for =  bothtrim(value);
	else if( tag == RLINKER )	linker_rev =  bothtrim(value);
	else if( tag == RPRIMER )	primer_rev =  bothtrim(value);
	else if( tag == MACHINE )	machine =  bothtrim(value);
	else if( tag == TYPE )		strFormat =  bothtrim(value);
	else if( tag == METHOD )		strMethod =  bothtrim(value);
	else if( tag == CUTOFF )		strCutoff =  bothtrim(value);
	else if( tag == LENGTH_STAT )	 	length_stat = bothtrim(value);

	else if( tag == CDHIT_SINGLETON )	 	cdhit_singleton = atoi(bothtrim(value).c_str());
	else if( tag == CDHIT_OTUS ) 			cdhit_otus = atoi(bothtrim(value).c_str());
	else if( tag == CDHIT_RF )				cdhit_rf = bothtrim(value) ;
	else if( tag == CDHIT_CHAO )			cdhit_chao = atof(bothtrim(value).c_str());
	else if( tag == CDHIT_CHAO_LCI )		cdhit_chao_lci = atof(bothtrim(value).c_str());
	else if( tag == CDHIT_CHAO_HCI )		cdhit_chao_hci =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_ACE )			cdhit_ace =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_ACE_LCI )		cdhit_ace_lci =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_ACE_HCI )		cdhit_ace_hci =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_JACKKNIFE )		cdhit_jackknife =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_JACKKNIFE_LCI )	cdhit_jackknife_lci =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_JACKKNIFE_HCI ) 	cdhit_jackknife_hci =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SHANNON ) 		cdhit_shannon =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SHANNON_LCI )	cdhit_shannon_lci =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SHANNON_HCI ) 	cdhit_shannon_hci =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_NPSHANNON ) 	cdhit_npshannon =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SIMPSON ) 		cdhit_simpson =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SIMPSON_LCI ) 	cdhit_simpson_lci =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SIMPSON_HCI )	cdhit_simpson_hci =  atof(bothtrim(value).c_str());

	else if( tag == READS_TARGET )    reads_target = atoi(bothtrim(value).c_str());
	else if( tag == CDHIT_SINGLETON_TARGET )    cdhit_singleton_target = atoi(bothtrim(value).c_str());
	else if( tag == CDHIT_OTUS_TARGET )       cdhit_otus_target = atoi(bothtrim(value).c_str());
	else if( tag == CDHIT_RF_TARGET )       cdhit_rf_target = bothtrim(value) ;
	else if( tag == CDHIT_CHAO_TARGET )     cdhit_chao_target = atof(bothtrim(value).c_str());
	else if( tag == CDHIT_CHAO_LCI_TARGET )   cdhit_chao_lci_target = atof(bothtrim(value).c_str());
	else if( tag == CDHIT_CHAO_HCI_TARGET )   cdhit_chao_hci_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_ACE_TARGET )      cdhit_ace_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_ACE_LCI_TARGET )    cdhit_ace_lci_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_ACE_HCI_TARGET )    cdhit_ace_hci_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_JACKKNIFE_TARGET )    cdhit_jackknife_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_JACKKNIFE_LCI_TARGET )  cdhit_jackknife_lci_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_JACKKNIFE_HCI_TARGET )  cdhit_jackknife_hci_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SHANNON_TARGET )    cdhit_shannon_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SHANNON_LCI_TARGET )  cdhit_shannon_lci_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SHANNON_HCI_TARGET )  cdhit_shannon_hci_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_NPSHANNON_TARGET )  cdhit_npshannon_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SIMPSON_TARGET )    cdhit_simpson_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SIMPSON_LCI_TARGET )  cdhit_simpson_lci_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_SIMPSON_HCI_TARGET )  cdhit_simpson_hci_target =  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_GOOD_LIBRARY)   cdhit_good_library=  atof(bothtrim(value).c_str());
	else if( tag == CDHIT_GOOD_LIBRARY_TARGET ) cdhit_good_library_target =  atof(bothtrim(value).c_str());

	else if( tag == BLAST_OTUS ) 		blast_otus = atoi(bothtrim(value).c_str());
	else if( tag == BLAST_RF )			blast_rf = bothtrim(value) ;
	else if( tag == BLAST_SINGLETON)     blast_singleton = atoi(bothtrim(value).c_str()) ;
	else if( tag == BLAST_CHAO )			blast_chao = atof(bothtrim(value).c_str());
	else if( tag == BLAST_CHAO_LCI )		blast_chao_lci = atof(bothtrim(value).c_str());
	else if( tag == BLAST_CHAO_HCI )		blast_chao_hci =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_ACE )			blast_ace =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_ACE_LCI )		blast_ace_lci =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_ACE_HCI )		blast_ace_hci =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_JACKKNIFE )		blast_jackknife =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_JACKKNIFE_LCI )	blast_jackknife_lci =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_JACKKNIFE_HCI )	blast_jackknife_hci =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SHANNON ) 		blast_shannon =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SHANNON_LCI )	blast_shannon_lci =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SHANNON_HCI ) 	blast_shannon_hci =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_NPSHANNON ) 	blast_npshannon =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SIMPSON ) 		blast_simpson =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SIMPSON_LCI ) 	blast_simpson_lci =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SIMPSON_HCI )	blast_simpson_hci =  atof(bothtrim(value).c_str());

	else if( tag == BLAST_OTUS_TARGET )      blast_otus_target = atoi(bothtrim(value).c_str());
	else if( tag == BLAST_RF_TARGET )        blast_rf_target = bothtrim(value) ;
	else if( tag == BLAST_SINGLETON_TARGET)     blast_singleton_target = atoi(bothtrim(value).c_str()) ;
	else if( tag == BLAST_CHAO_TARGET )      blast_chao_target = atof(bothtrim(value).c_str());
	else if( tag == BLAST_CHAO_LCI_TARGET )    blast_chao_lci_target = atof(bothtrim(value).c_str());
	else if( tag == BLAST_CHAO_HCI_TARGET )    blast_chao_hci_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_ACE_TARGET )     blast_ace_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_ACE_LCI_TARGET )   blast_ace_lci_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_ACE_HCI_TARGET )   blast_ace_hci_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_JACKKNIFE_TARGET )   blast_jackknife_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_JACKKNIFE_LCI_TARGET ) blast_jackknife_lci_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_JACKKNIFE_HCI_TARGET ) blast_jackknife_hci_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SHANNON_TARGET )     blast_shannon_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SHANNON_LCI_TARGET ) blast_shannon_lci_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SHANNON_HCI_TARGET )   blast_shannon_hci_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_NPSHANNON_TARGET )   blast_npshannon_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SIMPSON_TARGET )     blast_simpson_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SIMPSON_LCI_TARGET )   blast_simpson_lci_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_SIMPSON_HCI_TARGET ) blast_simpson_hci_target =  atof(bothtrim(value).c_str());
	else if( tag == BLAST_GOOD_LIBRARY)   blast_good_library=  atof(bothtrim(value).c_str());
	else if( tag == BLAST_GOOD_LIBRARY_TARGET ) blast_good_library_target =  atof(bothtrim(value).c_str());

	else if( tag == USER_GROUP_UID ) user_group_uid = atoi(bothtrim(value).c_str());
	else if( tag == USER_UID ) user_uid = atoi(bothtrim(value).c_str());
	else if( tag == USER_SAMPLE_NAME ) user_sample_name = bothtrim(value);
	else if( tag == SAMPLE_DESCRIPTION) sample_description = bothtrim(value);
	else if( tag == USER_SAMPLE_NAME) sampling_date = bothtrim(value);
	else if( tag == COORDS_NS ) coord_ns = bothtrim(value);
	else if( tag == COORDS_EW ) coord_ew = bothtrim(value);
	else if( tag == DATA_VERSION) data_version = atoi(bothtrim(value).c_str());
	else if( tag == N_CHIMERA) n_chimera = atoi(bothtrim(value).c_str());
	else if( tag == COMM_VER) comm_ver = atof(bothtrim(value).c_str());
	else if( tag == PHYLUM_COMP ) phylum_comp = bothtrim(value);
	else if( tag == SOURCE_KEYWORD ) source_keyword= bothtrim(value);
	else if( tag == CHUNLAB_UID ) chunlab_uid = bothtrim(value);

}

std::string CLHeader::getInfo(std::string tag){
	if( tag == DBINFO ) return db_info;
	if( tag == COUNT ) return int2str(reads);
	if( tag == MAX_LENGTH ) return int2str(max_len);
	if( tag == MIN_LENGTH ) return int2str(min_len);
	if( tag == AVG_LENGTH ) return dbl2str(avg_len);
	if( tag == CONTIGS ) return int2str(nContigs);
	if( tag == CLUSTER ) return int2str(nClusters);
	if( tag == BARCODE ) return barcode;
	if( tag == FLINKER ) return linker_for;
	if( tag == FPRIMER ) return primer_for;
	if( tag == RLINKER ) return linker_rev;
	if( tag == RPRIMER ) return primer_rev;
	if( tag == MACHINE ) return machine;
	if( tag == TYPE ) return strFormat;
	if( tag == METHOD ) return strMethod;
	if( tag == CUTOFF ) return strCutoff;

	if( tag.find_first_of( LENGTH_STAT ) )	 	return length_stat;

	if( tag.find_first_of( READS_TARGET ) != -1)   return int2str( reads_target );
	if( tag.find_first_of( CDHIT_SINGLETON ) != -1 )	 	return int2str( cdhit_singleton );
	if( tag.find_first_of( CDHIT_OTUS ) != -1 ) 			return int2str( cdhit_otus );
	if( tag.find_first_of( CDHIT_RF ) != -1 )				return cdhit_rf;
	if( tag.find_first_of( CDHIT_CHAO ) != -1 )			return dbl2str( cdhit_chao );
	if( tag.find_first_of( CDHIT_CHAO_LCI ) != -1 )		return dbl2str( cdhit_chao_lci );
	if( tag.find_first_of( CDHIT_CHAO_HCI ) != -1 )		return dbl2str( cdhit_chao_hci );
	if( tag.find_first_of( CDHIT_ACE ) != -1 )				return dbl2str( cdhit_ace );
	if( tag.find_first_of( CDHIT_ACE_LCI ) != -1 )			return dbl2str( cdhit_ace_lci );
	if( tag.find_first_of( CDHIT_ACE_HCI ) != -1 )			return dbl2str( cdhit_ace_hci );
	if( tag.find_first_of( CDHIT_JACKKNIFE ) != -1 )		return dbl2str( cdhit_jackknife );
	if( tag.find_first_of( CDHIT_JACKKNIFE_LCI ) != -1 )	return dbl2str( cdhit_jackknife_lci );
	if( tag.find_first_of( CDHIT_JACKKNIFE_HCI ) != -1 )	return dbl2str( cdhit_jackknife_hci );
	if( tag.find_first_of( CDHIT_SHANNON ) != -1 ) 		return dbl2str( cdhit_shannon );
	if( tag.find_first_of( CDHIT_SHANNON_LCI ) != -1 )		return dbl2str( cdhit_shannon_lci );
	if( tag.find_first_of( CDHIT_SHANNON_HCI ) != -1 ) 	return dbl2str( cdhit_shannon_hci );
	if( tag.find_first_of( CDHIT_NPSHANNON ) != -1 ) 		return dbl2str( cdhit_npshannon );
	if( tag.find_first_of( CDHIT_SIMPSON ) != -1 ) 		return dbl2str( cdhit_simpson );
	if( tag.find_first_of( CDHIT_SIMPSON_LCI ) != -1 ) 	return dbl2str( cdhit_simpson_lci );
	if( tag.find_first_of( CDHIT_SIMPSON_HCI ) != -1 )		return dbl2str( cdhit_simpson_hci );
	if( tag.find_first_of( CDHIT_GOOD_LIBRARY) != -1 )   return dbl2str( cdhit_good_library );
	if( tag.find_first_of( CDHIT_GOOD_LIBRARY_TARGET ) != -1 ) return dbl2str( cdhit_good_library_target );

	if( tag.find_first_of( BLAST_OTUS ) != -1 ) 			return int2str( blast_otus );
	if( tag.find_first_of( BLAST_RF ) != -1 )				return blast_rf;
	if( tag.find_first_of( BLAST_SINGLETON ) != -1 )        return blast_singleton +"";
	if( tag.find_first_of( BLAST_CHAO ) != -1 )			return dbl2str( blast_chao );
	if( tag.find_first_of( BLAST_CHAO_LCI ) != -1 )		return dbl2str( blast_chao_lci );
	if( tag.find_first_of( BLAST_CHAO_HCI ) != -1 )		return dbl2str( blast_chao_hci );
	if( tag.find_first_of( BLAST_ACE ) != -1 )				return dbl2str( blast_ace );
	if( tag.find_first_of( BLAST_ACE_LCI ) != -1 )			return dbl2str( blast_ace_lci );
	if( tag.find_first_of( BLAST_ACE_HCI ) != -1 )			return dbl2str( blast_ace_hci );
	if( tag.find_first_of( BLAST_JACKKNIFE ) != -1 )		return dbl2str( blast_jackknife );
	if( tag.find_first_of( BLAST_JACKKNIFE_LCI ) != -1 )	return dbl2str( blast_jackknife_lci );
	if( tag.find_first_of( BLAST_JACKKNIFE_HCI ) != -1 )	return dbl2str( blast_jackknife_hci );
	if( tag.find_first_of( BLAST_SHANNON ) != -1 ) 		return dbl2str( blast_shannon );
	if( tag.find_first_of( BLAST_SHANNON_LCI ) != -1 )		return dbl2str( blast_shannon_lci );
	if( tag.find_first_of( BLAST_SHANNON_HCI ) != -1 ) 	return dbl2str( blast_shannon_hci );
	if( tag.find_first_of( BLAST_NPSHANNON ) != -1 ) 		return dbl2str( blast_npshannon );
	if( tag.find_first_of( BLAST_SIMPSON ) != -1 ) 		return dbl2str( blast_simpson );
	if( tag.find_first_of( BLAST_SIMPSON_LCI ) != -1 ) 	return dbl2str( blast_simpson_lci );
	if( tag.find_first_of( BLAST_SIMPSON_HCI ) != -1 )		return dbl2str( blast_simpson_hci );

	if( tag.find_first_of( BLAST_OTUS_TARGET ) != -1 )      return int2str(blast_otus_target) +"";
	if( tag.find_first_of( BLAST_RF_TARGET ) != -1 )        return blast_rf_target +"";
	if( tag.find_first_of( BLAST_SINGLETON_TARGET ) != -1 )        return blast_singleton_target +"";
	if( tag.find_first_of( BLAST_CHAO_TARGET ) != -1 )      return dbl2str(blast_chao_target) +"";
	if( tag.find_first_of( BLAST_CHAO_LCI_TARGET ) != -1 )    return dbl2str(blast_chao_lci_target)+"";
	if( tag.find_first_of( BLAST_CHAO_HCI_TARGET ) != -1 )    return dbl2str(blast_chao_hci_target)+"";
	if( tag.find_first_of( BLAST_ACE_TARGET ) != -1 )     return dbl2str(blast_ace_target)+"";
	if( tag.find_first_of( BLAST_ACE_LCI_TARGET ) != -1 )   return dbl2str(blast_ace_lci_target)+"";
	if( tag.find_first_of( BLAST_ACE_HCI_TARGET )  != -1)   return dbl2str(blast_ace_hci_target)+"";
	if( tag.find_first_of( BLAST_JACKKNIFE_TARGET ) != -1 )   return dbl2str(blast_jackknife_target)+"";
	if( tag.find_first_of( BLAST_JACKKNIFE_LCI_TARGET ) != -1 ) return dbl2str(blast_jackknife_lci_target)+"";
	if( tag.find_first_of( BLAST_JACKKNIFE_HCI_TARGET ) != -1 ) return dbl2str(blast_jackknife_hci_target)+"";
	if( tag.find_first_of( BLAST_SHANNON_TARGET ) != -1 )     return dbl2str(blast_shannon_target)+"";
	if( tag.find_first_of( BLAST_SHANNON_LCI_TARGET ) != -1 ) return dbl2str(blast_shannon_lci_target)+"";
	if( tag.find_first_of( BLAST_SHANNON_HCI_TARGET ) != -1 )   return dbl2str(blast_shannon_hci_target)+"";
	if( tag.find_first_of( BLAST_NPSHANNON_TARGET ) != -1 )   return dbl2str(blast_npshannon_target)+"";
	if( tag.find_first_of( BLAST_SIMPSON_TARGET ) != -1 )     return dbl2str(blast_simpson_target)+"";
	if( tag.find_first_of( BLAST_SIMPSON_LCI_TARGET ) != -1 )   return dbl2str(blast_simpson_lci_target)+"";
	if( tag.find_first_of( BLAST_SIMPSON_HCI_TARGET ) != -1 ) return dbl2str(blast_simpson_hci_target)+"";
	if( tag.find_first_of( BLAST_GOOD_LIBRARY) != -1 )   return dbl2str(blast_good_library)+"";
	if( tag.find_first_of( BLAST_GOOD_LIBRARY_TARGET ) != -1 ) return dbl2str(blast_good_library_target)+"";

	if( tag.find_first_of( USER_GROUP_UID ) != -1 ) return int2str(user_group_uid);
	if( tag.find_first_of( USER_UID ) != -1 ) return int2str(user_uid);
	if( tag.find_first_of( USER_SAMPLE_NAME ) != -1 ) return user_sample_name;
	if( tag.find_first_of( SAMPLE_DESCRIPTION) != -1 ) return sample_description;
	if( tag.find_first_of( USER_SAMPLE_NAME) != -1 ) return sampling_date;
	if( tag.find_first_of( COORDS_NS ) != -1 ) return coord_ns ;
	if( tag.find_first_of( COORDS_EW ) != -1 ) return coord_ew ;
	if( tag.find_first_of( DATA_VERSION) != -1 ) return int2str( data_version );
	if( tag.find_first_of( N_CHIMERA ) != -1 ) return int2str( n_chimera );
	if( tag.find_first_of( COMM_VER ) != -1 ) return dbl2str( comm_ver );
	if( tag.find_first_of( PHYLUM_COMP ) != -1 ) return phylum_comp;
	if( tag.find_first_of( SOURCE_KEYWORD ) != -1 ) return source_keyword;
	if( tag.find_first_of( CHUNLAB_UID ) != -1 ) return chunlab_uid;



	return "";
}


//void CLHeader::writeToWriter(std::ostream writer){
//
//	if( meta != "" && meta.size() > 0 ) writer << meta;
//
//	if( db_info != "" ) writer << ( DBINFO + DELIMITER + db_info + NEW_LINE );
//	if( reads > -1 ) writer << ( COUNT + DELIMITER + int2str(reads) + NEW_LINE ); 
//	if( max_len > -1 && max_len != INT_MIN ) writer << ( MAX_LENGTH + DELIMITER + int2str(max_len) + NEW_LINE );
//	if( min_len > -1 && min_len != INT_MAX ) writer << ( MIN_LENGTH + DELIMITER + int2str(min_len) + NEW_LINE );
//	if( avg_len > -1 ) writer << ( AVG_LENGTH + DELIMITER + dbl2str(avg_len) + NEW_LINE );
//	if( nContigs > -1 ) writer << ( CONTIGS + DELIMITER + dbl2str(nContigs) + NEW_LINE );
//	if( nClusters > -1 ) writer << ( CLUSTER + DELIMITER + dbl2str(nClusters) + NEW_LINE );
//
//	if( barcode != "" ) writer << ( BARCODE + DELIMITER + barcode + NEW_LINE  );
//	if( linker_for != "" ) writer << ( FLINKER + DELIMITER + linker_for + NEW_LINE );
//	if( primer_for != "" ) writer << ( FPRIMER + DELIMITER + primer_for + NEW_LINE );
//	if( linker_rev != "" ) writer << ( RLINKER + DELIMITER + linker_rev + NEW_LINE );
//	if( primer_rev != "" ) writer << ( RPRIMER + DELIMITER + primer_rev + NEW_LINE );
//	if( machine != "" ) writer << ( MACHINE + DELIMITER + machine + NEW_LINE  );
//	if( strFormat != "" ) writer << ( TYPE + DELIMITER + strFormat + NEW_LINE );
//	if( strMethod != "" ) writer << ( METHOD + DELIMITER + strMethod + NEW_LINE );
//	if( strCutoff != "" ) writer << ( CUTOFF + DELIMITER + strCutoff + NEW_LINE );
//
//	if( length_stat != "" )      writer << ( LENGTH_STAT      + DELIMITER + length_stat + NEW_LINE );
//	if( cdhit_singleton  > -1 )      writer << ( CDHIT_SINGLETON    + DELIMITER + int2str( cdhit_singleton ) + NEW_LINE );      
//	if( cdhit_otus  > -1 )         writer << ( CDHIT_OTUS       + DELIMITER + int2str( cdhit_otus )    + NEW_LINE );            
//	if( cdhit_rf != "" )         writer << ( CDHIT_RF       + DELIMITER + cdhit_rf                + NEW_LINE );                       
//	if( cdhit_chao  > -1 )         writer << ( CDHIT_CHAO       + DELIMITER + dbl2str( cdhit_chao )     + NEW_LINE );           
//	if( cdhit_chao_lci  > -1 )       writer << ( CDHIT_CHAO_LCI     + DELIMITER + dbl2str( cdhit_chao_lci )   + NEW_LINE );       
//	if( cdhit_chao_hci  > -1 )       writer << ( CDHIT_CHAO_HCI     + DELIMITER + dbl2str( cdhit_chao_hci )   + NEW_LINE );       
//	if( cdhit_ace  > -1 )        writer << ( CDHIT_ACE      + DELIMITER + dbl2str( cdhit_ace )      + NEW_LINE );    
//	if( cdhit_ace_lci  > -1 )      writer << ( CDHIT_ACE_LCI    + DELIMITER + dbl2str( cdhit_ace_lci )    + NEW_LINE );
//	if( cdhit_ace_hci  > -1 )      writer << ( CDHIT_ACE_HCI    + DELIMITER + dbl2str( cdhit_ace_hci )    + NEW_LINE );
//	if( cdhit_jackknife  > -1 )      writer << ( CDHIT_JACKKNIFE    + DELIMITER + dbl2str( cdhit_jackknife )  + NEW_LINE );      
//	if( cdhit_jackknife_lci  > -1 )    writer << ( CDHIT_JACKKNIFE_LCI  + DELIMITER + dbl2str( cdhit_jackknife_lci )+ NEW_LINE );          
//	if( cdhit_jackknife_hci  > -1 )    writer << ( CDHIT_JACKKNIFE_HCI  + DELIMITER + dbl2str( cdhit_jackknife_hci )+ NEW_LINE );          
//	if( cdhit_shannon  > -1 )        writer << ( CDHIT_SHANNON    + DELIMITER + dbl2str( cdhit_shannon )    + NEW_LINE );        
//	if( cdhit_shannon_lci  > -1 )    writer << ( CDHIT_SHANNON_LCI  + DELIMITER + dbl2str( cdhit_shannon_lci )  + NEW_LINE );    
//	if( cdhit_shannon_hci  > -1 )      writer << ( CDHIT_SHANNON_HCI  + DELIMITER + dbl2str( cdhit_shannon_hci )  + NEW_LINE );            
//	if( cdhit_npshannon  > -1 )      writer << ( CDHIT_NPSHANNON    + DELIMITER + dbl2str( cdhit_npshannon )  + NEW_LINE );      
//	if( cdhit_simpson  > -1 )        writer << ( CDHIT_SIMPSON    + DELIMITER + dbl2str( cdhit_simpson )    + NEW_LINE );        
//	if( cdhit_simpson_lci  > -1 )      writer << ( CDHIT_SIMPSON_LCI  + DELIMITER + dbl2str( cdhit_simpson_lci )  + NEW_LINE );            
//	if( cdhit_simpson_hci  > -1 )    writer << ( CDHIT_SIMPSON_HCI  + DELIMITER + dbl2str( cdhit_simpson_hci )  + NEW_LINE );
//
//	if( reads_target  > -1 )       writer << ( READS_TARGET + DELIMITER + int2str( reads_target ) + NEW_LINE );
//	if( cdhit_singleton_target  > -1 )       writer << ( CDHIT_SINGLETON_TARGET    + DELIMITER + int2str( cdhit_singleton_target ) + NEW_LINE );      
//	if( cdhit_otus_target  > -1 )         writer << ( CDHIT_OTUS_TARGET       + DELIMITER + int2str( cdhit_otus_target )    + NEW_LINE);            
//	if( cdhit_rf_target != "" )         writer << ( CDHIT_RF_TARGET       + DELIMITER + cdhit_rf_target                + NEW_LINE );                       
//	if( cdhit_chao_target  > -1 )         writer << ( CDHIT_CHAO_TARGET       + DELIMITER + dbl2str( cdhit_chao_target )     + NEW_LINE );           
//	if( cdhit_chao_lci_target  > -1 )       writer << ( CDHIT_CHAO_LCI_TARGET     + DELIMITER + dbl2str( cdhit_chao_lci_target )   + NEW_LINE );       
//	if( cdhit_chao_hci_target  > -1 )       writer << ( CDHIT_CHAO_HCI_TARGET     + DELIMITER + dbl2str( cdhit_chao_hci_target )   + NEW_LINE );       
//	if( cdhit_ace_target  > -1 )        writer << ( CDHIT_ACE_TARGET      + DELIMITER + dbl2str( cdhit_ace_target )      + NEW_LINE );    
//	if( cdhit_ace_lci_target  > -1 )      writer << ( CDHIT_ACE_LCI_TARGET    + DELIMITER + dbl2str( cdhit_ace_lci_target )    + NEW_LINE );
//	if( cdhit_ace_hci_target  > -1 )      writer << ( CDHIT_ACE_HCI_TARGET    + DELIMITER + dbl2str( cdhit_ace_hci_target )    + NEW_LINE );
//	if( cdhit_jackknife_target  > -1 )      writer << ( CDHIT_JACKKNIFE_TARGET    + DELIMITER + dbl2str( cdhit_jackknife_target )  + NEW_LINE );      
//	if( cdhit_jackknife_lci_target  > -1 )    writer << ( CDHIT_JACKKNIFE_LCI_TARGET  + DELIMITER + dbl2str( cdhit_jackknife_lci_target )+ NEW_LINE );          
//	if( cdhit_jackknife_hci_target  > -1 )    writer << ( CDHIT_JACKKNIFE_HCI_TARGET  + DELIMITER + dbl2str( cdhit_jackknife_hci_target )+ NEW_LINE );          
//	if( cdhit_shannon_target  > -1 )        writer << ( CDHIT_SHANNON_TARGET    + DELIMITER + dbl2str( cdhit_shannon_target )    + NEW_LINE );        
//	if( cdhit_shannon_lci_target  > -1 )    writer << ( CDHIT_SHANNON_LCI_TARGET  + DELIMITER + dbl2str( cdhit_shannon_lci_target )  + NEW_LINE );    
//	if( cdhit_shannon_hci_target  > -1 )      writer << ( CDHIT_SHANNON_HCI_TARGET  + DELIMITER + dbl2str( cdhit_shannon_hci_target )  + NEW_LINE );            
//	if( cdhit_npshannon_target  > -1 )      writer << ( CDHIT_NPSHANNON_TARGET    + DELIMITER + dbl2str( cdhit_npshannon_target )  + NEW_LINE );      
//	if( cdhit_simpson_target  > -1 )        writer << ( CDHIT_SIMPSON_TARGET    + DELIMITER + dbl2str( cdhit_simpson_target )    + NEW_LINE );        
//	if( cdhit_simpson_lci_target  > -1 )      writer << ( CDHIT_SIMPSON_LCI_TARGET  + DELIMITER + dbl2str( cdhit_simpson_lci_target )  + NEW_LINE );            
//	if( cdhit_simpson_hci_target  > -1 )    writer << ( CDHIT_SIMPSON_HCI_TARGET  + DELIMITER + dbl2str( cdhit_simpson_hci_target )  + NEW_LINE );
//
//	if( cdhit_good_library > -1 )    writer << ( CDHIT_GOOD_LIBRARY + DELIMITER + dbl2str( cdhit_good_library)  + NEW_LINE );
//	if( cdhit_good_library_target  > -1 )      writer << ( CDHIT_GOOD_LIBRARY_TARGET  + DELIMITER + dbl2str( cdhit_good_library_target)  + NEW_LINE );            
//
//
//	if( blast_otus  > -1 )         writer << ( BLAST_OTUS       + DELIMITER + int2str( blast_otus )    + NEW_LINE );            
//	if( blast_rf != "" )         writer << ( BLAST_RF       + DELIMITER + blast_rf                + NEW_LINE );                       
//	if( blast_singleton > -1 )         writer << ( BLAST_SINGLETON       + DELIMITER + int2str(blast_singleton)                + NEW_LINE );
//	if( blast_chao  > -1 )         writer << ( BLAST_CHAO       + DELIMITER + dbl2str( blast_chao )     + NEW_LINE );           
//	if( blast_chao_lci  > -1 )       writer << ( BLAST_CHAO_LCI     + DELIMITER + dbl2str( blast_chao_lci )   + NEW_LINE );       
//	if( blast_chao_hci  > -1 )       writer << ( BLAST_CHAO_HCI     + DELIMITER + dbl2str( blast_chao_hci )   + NEW_LINE );       
//	if( blast_ace  > -1 )        writer << ( BLAST_ACE      + DELIMITER + dbl2str( blast_ace )      + NEW_LINE );    
//	if( blast_ace_lci  > -1 )      writer << ( BLAST_ACE_LCI    + DELIMITER + dbl2str( blast_ace_lci )    + NEW_LINE );
//	if( blast_ace_hci  > -1 )      writer << ( BLAST_ACE_HCI    + DELIMITER + dbl2str( blast_ace_hci )    + NEW_LINE );
//	if( blast_jackknife  > -1 )      writer << ( BLAST_JACKKNIFE    + DELIMITER + dbl2str( blast_jackknife )  + NEW_LINE );      
//	if( blast_jackknife_lci  > -1 )    writer << ( BLAST_JACKKNIFE_LCI  + DELIMITER + dbl2str( blast_jackknife_lci )+ NEW_LINE );          
//	if( blast_jackknife_hci  > -1 )    writer << ( BLAST_JACKKNIFE_HCI  + DELIMITER + dbl2str( blast_jackknife_hci )+ NEW_LINE );          
//	if( blast_shannon  > -1 )        writer << ( BLAST_SHANNON    + DELIMITER + dbl2str( blast_shannon )    + NEW_LINE );        
//	if( blast_shannon_lci  > -1 )    writer << ( BLAST_SHANNON_LCI  + DELIMITER + dbl2str( blast_shannon_lci )  + NEW_LINE );    
//	if( blast_shannon_hci  > -1 )      writer << ( BLAST_SHANNON_HCI  + DELIMITER + dbl2str( blast_shannon_hci )  + NEW_LINE );            
//	if( blast_npshannon  > -1 )      writer << ( BLAST_NPSHANNON    + DELIMITER + dbl2str( blast_npshannon )  + NEW_LINE );      
//	if( blast_simpson  > -1 )        writer << ( BLAST_SIMPSON    + DELIMITER + dbl2str( blast_simpson )    + NEW_LINE );        
//	if( blast_simpson_lci  > -1 )      writer << ( BLAST_SIMPSON_LCI  + DELIMITER + dbl2str( blast_simpson_lci )  + NEW_LINE );            
//	if( blast_simpson_hci  > -1 )    writer << ( BLAST_SIMPSON_HCI  + DELIMITER + dbl2str( blast_simpson_hci )  + NEW_LINE );
//
//	if( blast_otus_target  > -1 )         writer << ( BLAST_OTUS_TARGET       + DELIMITER + int2str( blast_otus_target )    + NEW_LINE );            
//	if( blast_rf_target != "" )         writer << ( BLAST_RF_TARGET       + DELIMITER + blast_rf_target                + NEW_LINE );
//	if( blast_singleton_target > -1 )         writer << ( BLAST_SINGLETON_TARGET       + DELIMITER + int2str(blast_singleton_target) + NEW_LINE );
//	if( blast_chao_target  > -1 )         writer << ( BLAST_CHAO_TARGET       + DELIMITER + dbl2str( blast_chao_target )     + NEW_LINE );           
//	if( blast_chao_lci_target  > -1 )       writer << ( BLAST_CHAO_LCI_TARGET     + DELIMITER + dbl2str( blast_chao_lci_target )   + NEW_LINE );       
//	if( blast_chao_hci_target  > -1 )       writer << ( BLAST_CHAO_HCI_TARGET     + DELIMITER + dbl2str( blast_chao_hci_target )   + NEW_LINE );       
//	if( blast_ace_target  > -1 )        writer << ( BLAST_ACE_TARGET      + DELIMITER + dbl2str( blast_ace_target )      + NEW_LINE );    
//	if( blast_ace_lci_target  > -1 )      writer << ( BLAST_ACE_LCI_TARGET    + DELIMITER + dbl2str( blast_ace_lci_target )    + NEW_LINE );
//	if( blast_ace_hci_target  > -1 )      writer << ( BLAST_ACE_HCI_TARGET    + DELIMITER + dbl2str( blast_ace_hci_target )    + NEW_LINE );
//	if( blast_jackknife_target  > -1 )      writer << ( BLAST_JACKKNIFE_TARGET    + DELIMITER + dbl2str( blast_jackknife_target )  + NEW_LINE );      
//	if( blast_jackknife_lci_target  > -1 )    writer << ( BLAST_JACKKNIFE_LCI_TARGET  + DELIMITER + dbl2str( blast_jackknife_lci_target )+ NEW_LINE );          
//	if( blast_jackknife_hci_target  > -1 )    writer << ( BLAST_JACKKNIFE_HCI_TARGET  + DELIMITER + dbl2str( blast_jackknife_hci_target )+ NEW_LINE );          
//	if( blast_shannon_target  > -1 )        writer << ( BLAST_SHANNON_TARGET    + DELIMITER + dbl2str( blast_shannon_target )    + NEW_LINE );        
//	if( blast_shannon_lci_target  > -1 )    writer << ( BLAST_SHANNON_LCI_TARGET  + DELIMITER + dbl2str( blast_shannon_lci_target )  + NEW_LINE );    
//	if( blast_shannon_hci_target  > -1 )      writer << ( BLAST_SHANNON_HCI_TARGET  + DELIMITER + dbl2str( blast_shannon_hci_target )  + NEW_LINE );            
//	if( blast_npshannon_target  > -1 )      writer << ( BLAST_NPSHANNON_TARGET    + DELIMITER + dbl2str( blast_npshannon_target )  + NEW_LINE );      
//	if( blast_simpson_target  > -1 )        writer << ( BLAST_SIMPSON_TARGET    + DELIMITER + dbl2str( blast_simpson_target )    + NEW_LINE );        
//	if( blast_simpson_lci_target  > -1 )      writer << ( BLAST_SIMPSON_LCI_TARGET  + DELIMITER + dbl2str( blast_simpson_lci_target )  + NEW_LINE );            
//	if( blast_simpson_hci_target  > -1 )    writer << ( BLAST_SIMPSON_HCI_TARGET  + DELIMITER + dbl2str( blast_simpson_hci_target )  + NEW_LINE );
//
//	if( blast_good_library  > -1 )    writer << ( BLAST_GOOD_LIBRARY + DELIMITER + dbl2str( blast_good_library )  + NEW_LINE );
//	if( blast_good_library_target  > -1 )      writer << ( BLAST_GOOD_LIBRARY_TARGET  + DELIMITER + dbl2str( blast_good_library_target )  + NEW_LINE );
//
//	if( user_group_uid > -1 ) writer << ( USER_GROUP_UID + DELIMITER + int2str(user_group_uid) + NEW_LINE );
//	if( user_uid > -1 ) writer << ( USER_UID + DELIMITER + int2str(user_uid) + NEW_LINE );
//	if( user_sample_name != "" ) writer << ( USER_SAMPLE_NAME+ DELIMITER + user_sample_name + NEW_LINE );
//	if( sample_description != "" ) writer << ( SAMPLE_DESCRIPTION + DELIMITER + sample_description + NEW_LINE );
//	if( sampling_date != "" ) writer << ( USER_SAMPLE_NAME+ DELIMITER + sampling_date + NEW_LINE );
//	if( coord_ns != "" ) writer << ( COORDS_NS + DELIMITER + coord_ns + NEW_LINE );
//	if( coord_ew != "" ) writer << ( COORDS_EW+ DELIMITER + coord_ew + NEW_LINE );
//	if( data_version > -1 ) writer << ( DATA_VERSION + DELIMITER + int2str(data_version) + NEW_LINE );
//	if( n_chimera > -1 ) writer << ( N_CHIMERA + DELIMITER + int2str(n_chimera ) + NEW_LINE );
//	if( comm_ver > -1 ) writer << ( COMM_VER + DELIMITER + dbl2str(comm_ver) + NEW_LINE );
//	if( phylum_comp != "" ) writer << ( PHYLUM_COMP + DELIMITER + phylum_comp + NEW_LINE );
//	if( source_keyword != "" ) writer << ( SOURCE_KEYWORD + DELIMITER + source_keyword + NEW_LINE );
//	if( chunlab_uid != "" ) writer << ( CHUNLAB_UID + DELIMITER + chunlab_uid + NEW_LINE );
//}


std::string CLHeader::writeToByteArray(){
	std::ostringstream writer;
	//writer.setf(std::ios::showpoint);

	if( meta != "" && meta.size() > 0 ) writer << meta;

	if( db_info != "" ) writer << ( DBINFO + DELIMITER + db_info + NEW_LINE );
	if( reads > -1 ) writer << ( COUNT + DELIMITER + int2str(reads) + NEW_LINE ); 
	if( max_len > -1 && max_len != INT_MIN ) writer << ( MAX_LENGTH + DELIMITER + int2str(max_len) + NEW_LINE );
	if( min_len > -1 && min_len != INT_MAX ) writer << ( MIN_LENGTH + DELIMITER + int2str(min_len) + NEW_LINE );
	if( avg_len > -1 ) writer << ( AVG_LENGTH + DELIMITER + dbl2str(avg_len) + NEW_LINE );
	if( nContigs > -1 ) writer << ( CONTIGS + DELIMITER + dbl2str(nContigs) + NEW_LINE );
	if( nClusters > -1 ) writer << ( CLUSTER + DELIMITER + dbl2str(nClusters) + NEW_LINE );

	if( barcode != "" ) writer << ( BARCODE + DELIMITER + barcode + NEW_LINE  );
	if( linker_for != "" ) writer << ( FLINKER + DELIMITER + linker_for + NEW_LINE );
	if( primer_for != "" ) writer << ( FPRIMER + DELIMITER + primer_for + NEW_LINE );
	if( linker_rev != "" ) writer << ( RLINKER + DELIMITER + linker_rev + NEW_LINE );
	if( primer_rev != "" ) writer << ( RPRIMER + DELIMITER + primer_rev + NEW_LINE );
	if( machine != "" ) writer << ( MACHINE + DELIMITER + machine + NEW_LINE  );
	if( strFormat != "" ) writer << ( TYPE + DELIMITER + strFormat + NEW_LINE );
	if( strMethod != "" ) writer << ( METHOD + DELIMITER + strMethod + NEW_LINE );
	if( strCutoff != "" ) writer << ( CUTOFF + DELIMITER + strCutoff + NEW_LINE );

	if( length_stat != "" )      writer << ( LENGTH_STAT      + DELIMITER + length_stat + NEW_LINE );
	if( cdhit_singleton  > -1 )      writer << ( CDHIT_SINGLETON    + DELIMITER + int2str( cdhit_singleton ) + NEW_LINE );      
	if( cdhit_otus  > -1 )         writer << ( CDHIT_OTUS       + DELIMITER + int2str( cdhit_otus )    + NEW_LINE );            
	if( cdhit_rf != "" )         writer << ( CDHIT_RF       + DELIMITER + cdhit_rf                + NEW_LINE );                       
	if( cdhit_chao  > -1 )         writer << ( CDHIT_CHAO       + DELIMITER + dbl2str( cdhit_chao )     + NEW_LINE );           
	if( cdhit_chao_lci  > -1 )       writer << ( CDHIT_CHAO_LCI     + DELIMITER + dbl2str( cdhit_chao_lci )   + NEW_LINE );       
	if( cdhit_chao_hci  > -1 )       writer << ( CDHIT_CHAO_HCI     + DELIMITER + dbl2str( cdhit_chao_hci )   + NEW_LINE );       
	if( cdhit_ace  > -1 )        writer << ( CDHIT_ACE      + DELIMITER + dbl2str( cdhit_ace )      + NEW_LINE );    
	if( cdhit_ace_lci  > -1 )      writer << ( CDHIT_ACE_LCI    + DELIMITER + dbl2str( cdhit_ace_lci )    + NEW_LINE );
	if( cdhit_ace_hci  > -1 )      writer << ( CDHIT_ACE_HCI    + DELIMITER + dbl2str( cdhit_ace_hci )    + NEW_LINE );
	if( cdhit_jackknife  > -1 )      writer << ( CDHIT_JACKKNIFE    + DELIMITER + dbl2str( cdhit_jackknife )  + NEW_LINE );      
	if( cdhit_jackknife_lci  > -1 )    writer << ( CDHIT_JACKKNIFE_LCI  + DELIMITER + dbl2str( cdhit_jackknife_lci )+ NEW_LINE );          
	if( cdhit_jackknife_hci  > -1 )    writer << ( CDHIT_JACKKNIFE_HCI  + DELIMITER + dbl2str( cdhit_jackknife_hci )+ NEW_LINE );          
	if( cdhit_shannon  > -1 )        writer << ( CDHIT_SHANNON    + DELIMITER + dbl2str( cdhit_shannon )    + NEW_LINE );        
	if( cdhit_shannon_lci  > -1 )    writer << ( CDHIT_SHANNON_LCI  + DELIMITER + dbl2str( cdhit_shannon_lci )  + NEW_LINE );    
	if( cdhit_shannon_hci  > -1 )      writer << ( CDHIT_SHANNON_HCI  + DELIMITER + dbl2str( cdhit_shannon_hci )  + NEW_LINE );            
	if( cdhit_npshannon  > -1 )      writer << ( CDHIT_NPSHANNON    + DELIMITER + dbl2str( cdhit_npshannon )  + NEW_LINE );      
	if( cdhit_simpson  > -1 )        writer << ( CDHIT_SIMPSON    + DELIMITER + dbl2str( cdhit_simpson )    + NEW_LINE );        
	if( cdhit_simpson_lci  > -1 )      writer << ( CDHIT_SIMPSON_LCI  + DELIMITER + dbl2str( cdhit_simpson_lci )  + NEW_LINE );            
	if( cdhit_simpson_hci  > -1 )    writer << ( CDHIT_SIMPSON_HCI  + DELIMITER + dbl2str( cdhit_simpson_hci )  + NEW_LINE );

	if( reads_target  > -1 )       writer << ( READS_TARGET + DELIMITER + int2str( reads_target ) + NEW_LINE );
	if( cdhit_singleton_target  > -1 )       writer << ( CDHIT_SINGLETON_TARGET    + DELIMITER + int2str( cdhit_singleton_target ) + NEW_LINE );      
	if( cdhit_otus_target  > -1 )         writer << ( CDHIT_OTUS_TARGET       + DELIMITER + int2str( cdhit_otus_target )    + NEW_LINE);            
	if( cdhit_rf_target != "" )         writer << ( CDHIT_RF_TARGET       + DELIMITER + cdhit_rf_target                + NEW_LINE );                       
	if( cdhit_chao_target  > -1 )         writer << ( CDHIT_CHAO_TARGET       + DELIMITER + dbl2str( cdhit_chao_target )     + NEW_LINE );           
	if( cdhit_chao_lci_target  > -1 )       writer << ( CDHIT_CHAO_LCI_TARGET     + DELIMITER + dbl2str( cdhit_chao_lci_target )   + NEW_LINE );       
	if( cdhit_chao_hci_target  > -1 )       writer << ( CDHIT_CHAO_HCI_TARGET     + DELIMITER + dbl2str( cdhit_chao_hci_target )   + NEW_LINE );       
	if( cdhit_ace_target  > -1 )        writer << ( CDHIT_ACE_TARGET      + DELIMITER + dbl2str( cdhit_ace_target )      + NEW_LINE );    
	if( cdhit_ace_lci_target  > -1 )      writer << ( CDHIT_ACE_LCI_TARGET    + DELIMITER + dbl2str( cdhit_ace_lci_target )    + NEW_LINE );
	if( cdhit_ace_hci_target  > -1 )      writer << ( CDHIT_ACE_HCI_TARGET    + DELIMITER + dbl2str( cdhit_ace_hci_target )    + NEW_LINE );
	if( cdhit_jackknife_target  > -1 )      writer << ( CDHIT_JACKKNIFE_TARGET    + DELIMITER + dbl2str( cdhit_jackknife_target )  + NEW_LINE );      
	if( cdhit_jackknife_lci_target  > -1 )    writer << ( CDHIT_JACKKNIFE_LCI_TARGET  + DELIMITER + dbl2str( cdhit_jackknife_lci_target )+ NEW_LINE );          
	if( cdhit_jackknife_hci_target  > -1 )    writer << ( CDHIT_JACKKNIFE_HCI_TARGET  + DELIMITER + dbl2str( cdhit_jackknife_hci_target )+ NEW_LINE );          
	if( cdhit_shannon_target  > -1 )        writer << ( CDHIT_SHANNON_TARGET    + DELIMITER + dbl2str( cdhit_shannon_target )    + NEW_LINE );        
	if( cdhit_shannon_lci_target  > -1 )    writer << ( CDHIT_SHANNON_LCI_TARGET  + DELIMITER + dbl2str( cdhit_shannon_lci_target )  + NEW_LINE );    
	if( cdhit_shannon_hci_target  > -1 )      writer << ( CDHIT_SHANNON_HCI_TARGET  + DELIMITER + dbl2str( cdhit_shannon_hci_target )  + NEW_LINE );            
	if( cdhit_npshannon_target  > -1 )      writer << ( CDHIT_NPSHANNON_TARGET    + DELIMITER + dbl2str( cdhit_npshannon_target )  + NEW_LINE );      
	if( cdhit_simpson_target  > -1 )        writer << ( CDHIT_SIMPSON_TARGET    + DELIMITER + dbl2str( cdhit_simpson_target )    + NEW_LINE );        
	if( cdhit_simpson_lci_target  > -1 )      writer << ( CDHIT_SIMPSON_LCI_TARGET  + DELIMITER + dbl2str( cdhit_simpson_lci_target )  + NEW_LINE );            
	if( cdhit_simpson_hci_target  > -1 )    writer << ( CDHIT_SIMPSON_HCI_TARGET  + DELIMITER + dbl2str( cdhit_simpson_hci_target )  + NEW_LINE );

	if( cdhit_good_library > -1 )    writer << ( CDHIT_GOOD_LIBRARY + DELIMITER + dbl2str( cdhit_good_library)  + NEW_LINE );
	if( cdhit_good_library_target  > -1 )      writer << ( CDHIT_GOOD_LIBRARY_TARGET  + DELIMITER + dbl2str( cdhit_good_library_target)  + NEW_LINE );            


	if( blast_otus  > -1 )         writer << ( BLAST_OTUS       + DELIMITER + int2str( blast_otus )    + NEW_LINE );            
	if( blast_rf != "" )         writer << ( BLAST_RF       + DELIMITER + blast_rf                + NEW_LINE );                       
	if( blast_singleton > -1 )         writer << ( BLAST_SINGLETON       + DELIMITER + int2str(blast_singleton)                + NEW_LINE );
	if( blast_chao  > -1 )         writer << ( BLAST_CHAO       + DELIMITER + dbl2str( blast_chao )     + NEW_LINE );           
	if( blast_chao_lci  > -1 )       writer << ( BLAST_CHAO_LCI     + DELIMITER + dbl2str( blast_chao_lci )   + NEW_LINE );       
	if( blast_chao_hci  > -1 )       writer << ( BLAST_CHAO_HCI     + DELIMITER + dbl2str( blast_chao_hci )   + NEW_LINE );       
	if( blast_ace  > -1 )        writer << ( BLAST_ACE      + DELIMITER + dbl2str( blast_ace )      + NEW_LINE );    
	if( blast_ace_lci  > -1 )      writer << ( BLAST_ACE_LCI    + DELIMITER + dbl2str( blast_ace_lci )    + NEW_LINE );
	if( blast_ace_hci  > -1 )      writer << ( BLAST_ACE_HCI    + DELIMITER + dbl2str( blast_ace_hci )    + NEW_LINE );
	if( blast_jackknife  > -1 )      writer << ( BLAST_JACKKNIFE    + DELIMITER + dbl2str( blast_jackknife )  + NEW_LINE );      
	if( blast_jackknife_lci  > -1 )    writer << ( BLAST_JACKKNIFE_LCI  + DELIMITER + dbl2str( blast_jackknife_lci )+ NEW_LINE );          
	if( blast_jackknife_hci  > -1 )    writer << ( BLAST_JACKKNIFE_HCI  + DELIMITER + dbl2str( blast_jackknife_hci )+ NEW_LINE );          
	if( blast_shannon  > -1 )        writer << ( BLAST_SHANNON    + DELIMITER + dbl2str( blast_shannon )    + NEW_LINE );        
	if( blast_shannon_lci  > -1 )    writer << ( BLAST_SHANNON_LCI  + DELIMITER + dbl2str( blast_shannon_lci )  + NEW_LINE );    
	if( blast_shannon_hci  > -1 )      writer << ( BLAST_SHANNON_HCI  + DELIMITER + dbl2str( blast_shannon_hci )  + NEW_LINE );            
	if( blast_npshannon  > -1 )      writer << ( BLAST_NPSHANNON    + DELIMITER + dbl2str( blast_npshannon )  + NEW_LINE );      
	if( blast_simpson  > -1 )        writer << ( BLAST_SIMPSON    + DELIMITER + dbl2str( blast_simpson )    + NEW_LINE );        
	if( blast_simpson_lci  > -1 )      writer << ( BLAST_SIMPSON_LCI  + DELIMITER + dbl2str( blast_simpson_lci )  + NEW_LINE );            
	if( blast_simpson_hci  > -1 )    writer << ( BLAST_SIMPSON_HCI  + DELIMITER + dbl2str( blast_simpson_hci )  + NEW_LINE );

	if( blast_otus_target  > -1 )         writer << ( BLAST_OTUS_TARGET       + DELIMITER + int2str( blast_otus_target )    + NEW_LINE );            
	if( blast_rf_target != "" )         writer << ( BLAST_RF_TARGET       + DELIMITER + blast_rf_target                + NEW_LINE );
	if( blast_singleton_target > -1 )         writer << ( BLAST_SINGLETON_TARGET       + DELIMITER + int2str(blast_singleton_target) + NEW_LINE );
	if( blast_chao_target  > -1 )         writer << ( BLAST_CHAO_TARGET       + DELIMITER + dbl2str( blast_chao_target )     + NEW_LINE );           
	if( blast_chao_lci_target  > -1 )       writer << ( BLAST_CHAO_LCI_TARGET     + DELIMITER + dbl2str( blast_chao_lci_target )   + NEW_LINE );       
	if( blast_chao_hci_target  > -1 )       writer << ( BLAST_CHAO_HCI_TARGET     + DELIMITER + dbl2str( blast_chao_hci_target )   + NEW_LINE );       
	if( blast_ace_target  > -1 )        writer << ( BLAST_ACE_TARGET      + DELIMITER + dbl2str( blast_ace_target )      + NEW_LINE );    
	if( blast_ace_lci_target  > -1 )      writer << ( BLAST_ACE_LCI_TARGET    + DELIMITER + dbl2str( blast_ace_lci_target )    + NEW_LINE );
	if( blast_ace_hci_target  > -1 )      writer << ( BLAST_ACE_HCI_TARGET    + DELIMITER + dbl2str( blast_ace_hci_target )    + NEW_LINE );
	if( blast_jackknife_target  > -1 )      writer << ( BLAST_JACKKNIFE_TARGET    + DELIMITER + dbl2str( blast_jackknife_target )  + NEW_LINE );      
	if( blast_jackknife_lci_target  > -1 )    writer << ( BLAST_JACKKNIFE_LCI_TARGET  + DELIMITER + dbl2str( blast_jackknife_lci_target )+ NEW_LINE );          
	if( blast_jackknife_hci_target  > -1 )    writer << ( BLAST_JACKKNIFE_HCI_TARGET  + DELIMITER + dbl2str( blast_jackknife_hci_target )+ NEW_LINE );          
	if( blast_shannon_target  > -1 )        writer << ( BLAST_SHANNON_TARGET    + DELIMITER + dbl2str( blast_shannon_target )    + NEW_LINE );        
	if( blast_shannon_lci_target  > -1 )    writer << ( BLAST_SHANNON_LCI_TARGET  + DELIMITER + dbl2str( blast_shannon_lci_target )  + NEW_LINE );    
	if( blast_shannon_hci_target  > -1 )      writer << ( BLAST_SHANNON_HCI_TARGET  + DELIMITER + dbl2str( blast_shannon_hci_target )  + NEW_LINE );            
	if( blast_npshannon_target  > -1 )      writer << ( BLAST_NPSHANNON_TARGET    + DELIMITER + dbl2str( blast_npshannon_target )  + NEW_LINE );      
	if( blast_simpson_target  > -1 )        writer << ( BLAST_SIMPSON_TARGET    + DELIMITER + dbl2str( blast_simpson_target )    + NEW_LINE );        
	if( blast_simpson_lci_target  > -1 )      writer << ( BLAST_SIMPSON_LCI_TARGET  + DELIMITER + dbl2str( blast_simpson_lci_target )  + NEW_LINE );            
	if( blast_simpson_hci_target  > -1 )    writer << ( BLAST_SIMPSON_HCI_TARGET  + DELIMITER + dbl2str( blast_simpson_hci_target )  + NEW_LINE );

	if( blast_good_library  > -1 )    writer << ( BLAST_GOOD_LIBRARY + DELIMITER + dbl2str( blast_good_library )  + NEW_LINE );
	if( blast_good_library_target  > -1 )      writer << ( BLAST_GOOD_LIBRARY_TARGET  + DELIMITER + dbl2str( blast_good_library_target )  + NEW_LINE );

	if( user_group_uid > -1 ) writer << ( USER_GROUP_UID + DELIMITER + int2str(user_group_uid) + NEW_LINE );
	if( user_uid > -1 ) writer << ( USER_UID + DELIMITER + int2str(user_uid) + NEW_LINE );
	if( user_sample_name != "" ) writer << ( USER_SAMPLE_NAME+ DELIMITER + user_sample_name + NEW_LINE );
	if( sample_description != "" ) writer << ( SAMPLE_DESCRIPTION + DELIMITER + sample_description + NEW_LINE );
	if( sampling_date != "" ) writer << ( USER_SAMPLE_NAME+ DELIMITER + sampling_date + NEW_LINE );
	if( coord_ns != "" ) writer << ( COORDS_NS + DELIMITER + coord_ns + NEW_LINE );
	if( coord_ew != "" ) writer << ( COORDS_EW+ DELIMITER + coord_ew + NEW_LINE );
	if( data_version > -1 ) writer << ( DATA_VERSION + DELIMITER + int2str(data_version) + NEW_LINE );
	if( n_chimera > -1 ) writer << ( N_CHIMERA + DELIMITER + int2str(n_chimera ) + NEW_LINE );
	if( comm_ver > -1 ) writer << COMM_VER + DELIMITER  << dbl2str(comm_ver) <<  NEW_LINE ;
	if( phylum_comp != "" ) writer << ( PHYLUM_COMP + DELIMITER + phylum_comp + NEW_LINE );
	if( source_keyword != "" ) writer << ( SOURCE_KEYWORD + DELIMITER + source_keyword + NEW_LINE );
	if( chunlab_uid != "" ) writer << ( CHUNLAB_UID + DELIMITER + chunlab_uid + NEW_LINE );

	return writer.str();
}

void CLHeader::setReadCount(int nTmp ){
	reads = nTmp;
}

void CLHeader::setMaxLen(int nTmp){
	max_len = nTmp;
}

void CLHeader::setMinLen(int nTmp){
	min_len = nTmp;
}

void CLHeader::setAvgLen(double dTmp){
	avg_len = dTmp;
}

void CLHeader::setContigCount(int nTmp){
	nContigs = nTmp;
}

int CLHeader::getDroppedRawReads(){
	return nDroppedRawRead;
}

int CLHeader::getUnsortedRawReads(){
	return nUnsortedRawRead;
}

int CLHeader::getSortedRawReads(){
	return nSortedRawRead;
}

int CLHeader::getTotalRawReads(){
	return nTotalRawRead;
}

