// * * IgFamily v0.12.2s * * 
//
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#include <iostream>
#include <map>
#include <omp.h>
#include <string>
#include <vector>

#include "IgFamily.h"
#include "fpf_core.h"
#include "fpf_filesystem.h"
#include "fpf_interface.h"
#include "fpf_report.h"


int main() {

	using std::string;
	using std::vector;

	using fpf_filesystem::filesystem;

	std::cout << "\n -- IgFamily " << IgFamily::version << " --\n";

	string select_FASTA{ IgFamily::DEFAULT_INPUT_FASTA };
	vector<string> v_select_peptide_assignment{ IgFamily::DEFAULT_PEPTIDE_ASSIGNMENT_METHOD };

	fpf_interface::display_settings(select_FASTA, v_select_peptide_assignment);
	fpf_interface::select_settings(select_FASTA, v_select_peptide_assignment);

	vector<string> v_root_directory{ fpf_filesystem::read_root_dir(IgFamily::IGFAMILY_ROOT_DIR) };
	vector<filesystem> v_filesystem{ fpf_filesystem::read_filesystem(v_root_directory) };

	for (auto& itr_v_filesystem : v_filesystem) {
		if (itr_v_filesystem.perform_wiff_fileconversion) {
			fpf_core::core_perform_wiff_fileconversion(itr_v_filesystem);
		}
		if (itr_v_filesystem.perform_novor_denovo) {
			fpf_core::core_perform_novor_denovo(itr_v_filesystem);
		}
		bool FASTA_exists = fpf_core::core_create_FASTA_data_structures(itr_v_filesystem, select_FASTA);
		if (FASTA_exists) {
			fpf_core::core_parse_data(itr_v_filesystem, v_select_peptide_assignment);
			for (auto& itr_v_sample_analysis : itr_v_filesystem.v_sample_analysis) {
				if (itr_v_sample_analysis.file_found) {
					fpf_core::core_create_data_structures(itr_v_filesystem, itr_v_sample_analysis);
					fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_sample_analysis, false);
					fpf_core::core_data_analysis(itr_v_filesystem, itr_v_sample_analysis, false);
					fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_sample_analysis, true);
					fpf_core::core_data_analysis(itr_v_filesystem, itr_v_sample_analysis, true);
					fpf_core::core_multinomial(itr_v_filesystem, itr_v_sample_analysis);
					fpf_core::core_report(itr_v_filesystem, itr_v_sample_analysis);
					itr_v_filesystem = {};
				}
			}
		}
	}

	string farewell;
	std::cout << "\n\n\n program complete...";
	std::cout << "\n\n\n input any key to exit...\n\n > ";
	std::cin >> farewell;

	return EXIT_SUCCESS;
}
