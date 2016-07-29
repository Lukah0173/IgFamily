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

	string_type version = "v0.5.4";

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

	value_type alpha_sum = 6.4; // Parameter of prior on doc-topic distribution"
	value_type beta = 0.01; // Parameter of prior on topic-word distribution"
	value_type cost = 1.6; // C param in SVM
	value_type ell = 64.0; // Margin param in SVM: usually 1
	size_type num_iter = 50; // Number of burn - in iterations for MCMC
	size_type st_constcount_topic = 40; // Model size, usually called K
	size_type eval_every = 100; //Evaluate the model every N iterations
	size_type num_mh = 6; // Number of MH steps for each token
	size_type num_gibbs = 2; // Number of Gibbs sampling for classifier
	size_type top = 10; // Save top N words for each topic

	const value_type NUSQINV = 1.0;
	const size_type  TRAIN_COLLECT = 1;
	const size_type  TEST_COLLECT = 20;
	const size_type  MAX_TEST_BURNIN = 100;
	const value_type GIBBS_CONVERGED = 1e-4;
	const value_type EXP_THRESHOLD = 640.0;

	char_type ct_parse_data_separator = ':';
	char_type ct_parse_data_delimitor = ' ';
}

#endif
