// * * fpf_fasta.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_FASTA_CREATOR
#define	FPF_FASTA_CREATOR

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "IgFamily.h"


namespace fpf_fasta_creator {

	using std::string;
	using std::vector;

	struct fasta_creator_directory {
	public:
		vector<string> v_fasta_module;
	};

	vector<string> read_FASTA_module_directory(const string& par_FASTA_module_directory) {
		vector<string> temp_v_FASTA_module{};
		string temp_FASTA_module{};
		const string FASTA_module_directory{ IgFamily::DEFAULT_FASTA_MODULE_DIRECTORY + "_FASTA_module_directory.txt" };
		std::ifstream fin_FASTA_module_directory(FASTA_module_directory);
		char read_FASTA_module_directory{};
		while (fin_FASTA_module_directory.get(read_FASTA_module_directory)) {
			if (read_FASTA_module_directory == ',') {
				temp_v_FASTA_module.push_back(temp_FASTA_module);
				temp_FASTA_module.clear();
			}
			else {
				temp_FASTA_module += read_FASTA_module_directory;
			}
		}
		return temp_v_FASTA_module;
	}
}

#endif
