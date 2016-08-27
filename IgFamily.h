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

	const string version = "v0.7.13";

	const bool FILESYSTEM_UPDATE_ALL = 1;
	const bool MAP_FOUT_BY_SCORE = 1;
	const bool MAP_FOUT_BY_DISTINCT = 1;
	const bool MAP_FOUT_PEPTIDE_SUMMARY_BY_SPECTRALCOUNT = 1;
	const bool DEBUG_MODE = 0;
	const bool OUTPUT_FASTA = 1;

	const string IGFAMILY_ROOT_DIR = "IgFamily_root_4.txt";
	const string DEFAULT_INPUT_FASTA = "IGHV_IGLV_IGKV_CONT_20160827.fasta";
	const string DEFAULT_PEPTIDE_ASSIGNMENT_METHOD = "PEAKS de novo";

	const double DENOVO_PEPTIDE_CONFIDENCE_THRESHOLD = 50;
	const double BLASTP_THRESHOLD = 10000;
	const double BLASTP_EVALUE_THRESHOLD = 20;
	const double BLASTP_PARPROP_SCALE = BLASTP_EVALUE_THRESHOLD;
	const double BLASTP_EVALUETRANSFORMED_THRESHOLD = 10;
	const double MULTINOMIAL_ELEMENT_OUTPUT_THRESHOLD = 0.1;







	int brute_calc(int i) {
		std::vector<int> v(1000000);
		std::iota(v.begin(), v.end(), 0);
		if (std::find(v.begin(), v.end(), 999999) == v.end()) {
			i = -i;
		}
		return (i);
	}

	void async_test() {
		std::vector<std::future<int>> test_future{};

		for (auto i = 1; i < 50000; ++i) {
			test_future.push_back(std::async(brute_calc, i));
		}

		for (auto& get_future : test_future) {
			std::cout << get_future.get() << std::endl;
		}
	}
}

#endif
