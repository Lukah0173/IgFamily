// * * IgFamily.h * *
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *


#include <cstdlib>
#include <string> // provides - std::string

#ifndef IgFamily
#define	IgFamily

namespace IgFamily {

	using std::string;

	const string version = "v0.7.9";

	const bool FILESYSTEM_UPDATE_ALL = 1;
	const bool NOVOR_DENOVO = 0;
	const bool MAP_FOUT_BY_SCORE = 1;
	const bool MAP_FOUT_BY_DISTINCT = 1;
	const bool MAP_FOUT_PEPTIDE_SUMMARY_BY_SPECTRALCOUNT = 1;
	const bool DEBUG_MODE = 0;
	const bool OUTPUT_FASTA = 1;

	const string IGFAMILY_ROOT_DIR = "IgFamily_root_4.txt";
	const string DEFAULT_INPUT_FASTA = "MOST_RECENT_20160804.fasta";
	const string DEFAULT_INPUT_FASTA_DIRECTORY = "FASTA\\" + DEFAULT_INPUT_FASTA;

	const double DENOVO_PEPTIDE_CONFIDENCE_THRESHOLD = 50;
	const double BLASTP_THRESHOLD = 10000;
	const double BLASTP_EVALUE_THRESHOLD = 20;
	const double BLASTP_PARPROP_SCALE = BLASTP_EVALUE_THRESHOLD;
	const double BLASTP_EVALUETRANSFORMED_THRESHOLD = 10;
	const double MULTINOMIAL_ELEMENT_OUTPUT_THRESHOLD = 0.1;
}

#endif
