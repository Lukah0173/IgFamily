// * * IgFamily v0.12.0 * * 
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

	// temp

	//vector<std::multimap<double, double>> v2_combined_contaminants_list{};
	//for (auto& itr_v_filesystem : v_filesystem) {
	//	for (auto& itr_sample_analysis : itr_v_filesystem.v_sample_analysis) {
	//		v2_combined_contaminants_list.push_back(fpf_report::fout_multinomial_contaminants_list(itr_v_filesystem, itr_sample_analysis));
	//	}
	//}

	//std::multimap<double, double> v_selected_combined_contaminants_list{};
	//bool hit_sample{};
	//size_t count_sample{};
	//for (auto i = 0; i < v2_combined_contaminants_list.size(); ++i) {
	//	for (auto& itr_combined_contaminants_list_map : v2_combined_contaminants_list[i]) {
	//		for (auto j = 0; j < v2_combined_contaminants_list.size(); ++j) {
	//			if (i != j) {
	//				if (!hit_sample) {
	//					auto& temp_find = v2_combined_contaminants_list[j].find(itr_combined_contaminants_list_map.first);
	//					if (temp_find != v2_combined_contaminants_list[j].end()) {
	//						++count_sample;
	//						if (count_sample == 4) {
	//							v_selected_combined_contaminants_list.insert(itr_combined_contaminants_list_map);
	//						}
	//						hit_sample = true;
	//					}
	//				}
	//			}
	//			hit_sample = {};
	//		}
	//		count_sample = {};
	//	}
	//}

	//std::map<double, double> v_averaged_selected_combined_contaminants_list{};
	//for (const auto& itr_v_selected_combined_contaminants_list : v_selected_combined_contaminants_list) {
	//	auto& temp_find = v_averaged_selected_combined_contaminants_list.find(itr_v_selected_combined_contaminants_list.first);
	//	if (temp_find != v_averaged_selected_combined_contaminants_list.end()) {
	//		temp_find->second = (itr_v_selected_combined_contaminants_list.second + temp_find->second) / double(2);
	//	}
	//	else {
	//		v_averaged_selected_combined_contaminants_list[itr_v_selected_combined_contaminants_list.first] = itr_v_selected_combined_contaminants_list.second;
	//	}
	//}

	//for (auto& itr_v_filesystem : v_filesystem) {
	//	for (auto& itr_v_sample_analysis : itr_v_filesystem.v_sample_analysis) {
	//		fpf_report::temp_fout_multinomial_contaminants_list(itr_v_filesystem, itr_v_sample_analysis, v_selected_combined_contaminants_list);
	//	}
	//}



	string farewell;
	std::cout << "\n\n\n program complete...";
	std::cout << "\n\n\n input any key to exit...\n\n > ";
	std::cin >> farewell;

	return EXIT_SUCCESS;
}
