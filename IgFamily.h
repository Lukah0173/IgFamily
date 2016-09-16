// * * IgFamily.h * *
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef IgFamily
#define	IgFamily

#include <string>


#include <cstdlib>
#include <future>
#include <iostream>
#include <numeric>



namespace IgFamily {

	using std::string;

	const string version = "v0.9.5f";

	const bool FILESYSTEM_MODE = 0;
	const bool FILESYSTEM_UPDATE_ALL = 1;
	const bool MAP_FOUT_BY_SCORE = 1;
	const bool MAP_FOUT_BY_DISTINCT = 1;
	const bool MAP_FOUT_PEPTIDE_SUMMARY_BY_SPECTRALCOUNT = 1;
	const bool OUTPUT_FASTA = 1;
	const bool BLASTP_BY_SELECTED_PEPTIDE = 1;

	const string IGFAMILY_ROOT_DIR = "IgFamily_root_3.txt";
	const string DEFAULT_IGFAMILY_DIRECTORY = "";
	//const string DEFAULT_IGFAMILY_DIRECTORY = "Z:\\Lukah_Dykes\\IgFamily\\";
	const string DEFAULT_FASTA_DIRECTORY = DEFAULT_IGFAMILY_DIRECTORY + "FASTA\\";
	const string DEFAULT_FASTA_MODULE_DIRECTORY = DEFAULT_FASTA_DIRECTORY + "FASTA_modules\\";
	const string DEFAULT_BLASTP_DIRECTORY = DEFAULT_IGFAMILY_DIRECTORY + "blast_directory\\";
	const string DEFAULT_MSCONVERT_DIRECTORY = DEFAULT_IGFAMILY_DIRECTORY + "msconvert\\";
	const string DEFAULT_NOVOR_DIRECTORY = DEFAULT_IGFAMILY_DIRECTORY + "novor\\win\\";
	const string DEFAULT_TRANSCRIPT_DIRECTORY = DEFAULT_IGFAMILY_DIRECTORY + "transcript_data\\";
	const string DEFAULT_GENOME_DIRECTORY = DEFAULT_IGFAMILY_DIRECTORY + "genome_data\\";
	const string DEFAULT_INPUT_FASTA = "IGHV_IGLV_IGKV_CONT_20160827.fasta";
	const string DEFAULT_PEPTIDE_ASSIGNMENT_METHOD = "PEAKS de novo";

	const int OUTPUT_FASTA_ACCESSION_WIDTH = 60;
	const double DENOVO_PEPTIDE_CONFIDENCE_THRESHOLD = 50;
	const double DENOVO_LOCAL_CONFIDENCE_MOVING_AVERAGE_THRESHOLD = 85;
	const double DENOVO_PEPTIDE_SIZE_THRESHOLD = 5;
	const double BLASTP_THRESHOLD = 10000;
	const double BLASTP_EVALUE_THRESHOLD = 20;
	const double BLASTP_PARPROP_SCALE = BLASTP_EVALUE_THRESHOLD;
	const double BLASTP_EVALUETRANSFORMED_THRESHOLD = 10;
	const string SELECT_TYPE_GENE_FAMILIES = "IG";
	const double SELECT_N_MANY_GENE_FAMILIES = 50;
	const double MULTINOMIAL_ELEMENT_OUTPUT_THRESHOLD = 0.1;
	const double REPORT_SCORE_THRESHOLD = 10;
	const double REPORT_QUERY_ALIGNMENT_TOTALSCORE_OUTPUT_THRESHOLD = 0.1;
	const double REPORT_QUERY_ALIGNMENT_PARSCORE_OUTPUT_THRESHOLD = 0.1;
	const double REPORT_QUERY_EVALUETRANSFORMED_THRESHOLD = 10;

	bool POLYMORPHISM_SELECTED = 0;
	double MULTINOMIAL_CONJUGATION_FACTOR = 0.5;
	double MULTINOMIAL_CONJUGATION_FACTOR_LIMIT_1 = 50;
	double MULTINOMIAL_CONJUGATION_FACTOR_CONVERGE_1 = 0.05;
	double MULTINOMIAL_CONJUGATION_FACTOR_LIMIT_2 = 50;
	double MULTINOMIAL_CONJUGATION_FACTOR_CONVERGE_2 = 0.01;




	int brute_calc(int i) {
		std::vector<int> v(10000000);
		std::iota(v.begin(), v.end(), 0);
		if (std::find(v.begin(), v.end(), 500) == v.end()) {
			i = -i;
		}
		return (i);
	}

	void async_test() {
		std::vector<std::future<int>> test_future{};
		for (auto i = 1; i < 5000000; ++i) {
			test_future.push_back(std::async(brute_calc, i));
		}
		for (auto& get_future : test_future) {
			std::cout << get_future.get() << std::endl;
		}
	}
}

#endif
