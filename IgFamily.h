// * * IgFamily.h * *
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *


#define EXIT_FILESYTEM_CURRENT 0;

#include <string> // provides - string
#include <cstdlib> // provides - size_t

#ifndef IgFamily
#define	IgFamily

namespace IgFamily {

	using std::string;

	const string version = "v0.7.4";

	const bool FILESYSTEM_UPDATE_ALL = 1;
	const bool NOVOR_DENOVO = 1;
	const bool MAP_FOUT_BY_SCORE = 1;
	const bool MAP_FOUT_BY_DISTINCT = 1;
	const bool MAP_FOUT_PEPTIDE_SUMMARY_BY_SPECTRALCOUNT = 1;
	const bool DEBUG_MODE = 0;
	const bool OUTPUT_FASTA = 1;

	const string IGFAMILY_ROOT_DIR = "IgFamily_root_3.txt";
	const string DEFAULT_INPUT_FASTA = "LUKAH_WITHUNIPROT_20160815.fasta";
	const string DEFAULT_INPUT_FASTA_DIRECTORY = "FASTA\\" + DEFAULT_INPUT_FASTA;

	const double BLASTP_THRESHOLD = 10000;
	const double EVALUE_THRESHOLD = 20;
	const double PARPROP_SCALE = EVALUE_THRESHOLD;
	const double MULTINOMIAL_ELEMENT_OUTPUT_THRESHOLD = 0.1;
}

#endif
