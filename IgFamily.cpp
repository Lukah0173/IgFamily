// * * IgFamily v0.9.3b * * 
//
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <string>
#include <utility>

#include "IgFamily.h"
#include "fpf_core.h"
#include "fpf_filesystem.h"
#include "fpf_filesystem_analysis.h"
#include "fpf_homology_analysis.h"
#include "fpf_interface.h"
#include "fpf_multinomial.h"
#include "fpf_report.h"
#include "fpf_utility.h"


int main() {

	using std::string;
	using std::vector;

	std::cout << "\n -- IgFamily " << IgFamily::version << " --\n";

	string select_FASTA{ IgFamily::DEFAULT_INPUT_FASTA };
	vector<string> v_select_peptide_assignment{ IgFamily::DEFAULT_PEPTIDE_ASSIGNMENT_METHOD };

	fpf_interface::display_settings(select_FASTA, v_select_peptide_assignment);
	fpf_interface::select_settings(select_FASTA, v_select_peptide_assignment);

	vector<string> v_root_directory = fpf_filesystem::read_root_dir(IgFamily::IGFAMILY_ROOT_DIR);
	vector<fpf_filesystem::filesystem> v_filesystem = fpf_filesystem::read_filesystem(v_root_directory);

	for (auto& itr_v_filesystem : v_filesystem) {
		fpf_core::core_perform_wiff_fileconversion(itr_v_filesystem);
		fpf_core::core_perform_novor_denovo(itr_v_filesystem);
		fpf_core::core_create_peptide_data_structures(itr_v_filesystem, v_select_peptide_assignment);
		for (auto& itr_v_sample_analysis : itr_v_filesystem.v_sample_analysis) {
			if (itr_v_sample_analysis.file_found) {
				std::ifstream fin_FASTA(select_FASTA);
				fpf_parse::custom_FASTA_output(select_FASTA);
				std::cout << "\n\n\n * parsing FASTA file... \n";
				vector<fpf_parse::FASTA_data> main_FASTA = fpf_parse::parse_FASTA(fin_FASTA);

				if (check_FASTA_file_empty(main_FASTA)) { return 1; };
				fpf_parse::output_v_FASTA_data(main_FASTA);
				fpf_parse::output_v_FASTA_data_to_blastdirectory(main_FASTA);

				std::cout << "\n\n peptides parsed - " << itr_v_sample_analysis.main_v_csv_peptides.size();
				std::cout << "\n\n\n\n creating data structures for file ";
				std::cout << itr_v_filesystem.filename;
				std::cout << "...";
				itr_v_sample_analysis.v_peptide_data = fpf_data::create_v_peptide_data(itr_v_sample_analysis.main_v_csv_peptides);
				itr_v_sample_analysis.v_peptide_data_map = fpf_data::create_v_peptide_data_map(itr_v_sample_analysis.v_peptide_data);
				itr_v_sample_analysis.v_peptide_analysis = fpf_data_analysis::create_v_peptide_analysis(itr_v_sample_analysis.v_peptide_data);
				itr_v_sample_analysis.v_peptide_analysis_map = fpf_data_analysis::create_v_peptide_analysis_map(itr_v_sample_analysis.v_peptide_analysis);
				itr_v_sample_analysis.v_protein_data = fpf_data::create_v_protein_data(main_FASTA);
				itr_v_sample_analysis.v_protein_data_map = fpf_data::create_v_protein_data_map(itr_v_sample_analysis.v_protein_data);
				std::cout << "\n ...data structures assigned";

				fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_sample_analysis, false);
				fpf_core::core_data_analysis(itr_v_filesystem, itr_v_sample_analysis);
				fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_sample_analysis, true);
				fpf_core::core_data_analysis(itr_v_filesystem, itr_v_sample_analysis);
				fpf_core::core_multinomial(itr_v_filesystem, itr_v_sample_analysis);
				fpf_core::core_report(itr_v_filesystem, itr_v_sample_analysis);
				fpf_filesystem::fout_filesystem(itr_v_filesystem);
			}
			else {
				std::cout << "\n ...no data file found";
			}
		}
	}

	string farewell;
	std::cout << "\n\n\n program complete...";
	std::cout << "\n\n\n input any key to exit...\n\n > ";
	std::cin >> farewell;

	return 0;
}
