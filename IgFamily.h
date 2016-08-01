// * * IgFamily.h * *
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *


#define EXIT_FILESYTEM_CURRENT 0;

#include <string> // provides - std::string
#include <vector> // provides - std::vector
#include <cstdlib> // provides - size_t

#ifndef IgFamily
#define	IgFamily

namespace IgFamily {

	typedef std::string string_type;
	typedef char char_type;
	typedef double value_type;
	typedef size_t size_type;

	string_type version = "v0.5.11";

	bool FILESYSTEM_MODE = true;
	bool FILESYSTEM_UPDATE_ALL = true;
	bool SIMPLE_SCORE = false;
	size_type MAP_FOUT_BY_SCORE = 1;
	size_type MAP_FOUT_BY_DISTINCT = 1;
	size_type MAP_FOUT_PEPTIDE_SUMMARY_BY_SPECTRALCOUNT = 1;
	size_type DEBUG_MODE = 0;
	size_type OUTPUT_FASTA = 0;

	string_type IGFAMILY_ROOT_DIR = "IgFamily_root_2.txt";
	string_type INPUT_CSV;
	string_type INPUT_FASTA = "FASTA\\FPF_V_mouse_20160620.fasta";
	string_type output = "z";

	double PARSE_THRESHOLD_IgP = 30;
	double SCORE_MEAN = 0;
	double SCORE_THRESHOLD = 10;
	size_type GLOBAL_ITERATOR = 0;
	size_type ITERATE_TRAIN_SCORE = 1;
	double BLASTP_THRESHOLD = 10000;
	double PARPROP_SCALE = 10;

	char_type ct_parse_data_separator = ':';
	char_type ct_parse_data_delimitor = ' ';
}

#endif
