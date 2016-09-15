// * * IgFamily v0.9.3 * * 
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
#include "fpf_data_analysis.h"
#include "fpf_dirichlet_mixture_model.h"
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
	}

	for (auto& itr_v_filesystem : v_filesystem) {
		bool filesystem_modified{};
		vector<fpf_parse::csv_data> main_v_csv_PEAKS_database_peptides;
		vector<fpf_parse::csv_data> main_v_csv_PEAKS_denovo_peptides;
		vector<fpf_parse::csv_data> main_v_csv_NOVOR_denovo_peptides;
		if ((IgFamily::FILESYSTEM_UPDATE_ALL) || (itr_v_filesystem.fileversion != IgFamily::version)) {
			std::cout << "\n\n\n\n parsing data for file ";
			std::cout << itr_v_filesystem.filename;
			std::cout << "...\n";
			main_v_csv_PEAKS_database_peptides = fpf_filesystem::parse_filesystem_PEAKS_database_peptides(fpf_filesystem::read_filesystem_PEAKS_database_peptides(itr_v_filesystem.directory));
			main_v_csv_PEAKS_denovo_peptides = fpf_filesystem::parse_filesystem_PEAKS_denovo_peptides(fpf_filesystem::read_filesystem_PEAKS_denovo_peptides(itr_v_filesystem.directory));
			main_v_csv_NOVOR_denovo_peptides = fpf_filesystem::parse_filesystem_NOVOR_denovo_peptides(fpf_filesystem::read_filesystem_NOVOR_denovo_peptides(itr_v_filesystem.directory));
			filesystem_modified = true;
		}

		bool file_found{};
		for (const auto& itr_v_select_spectra_assignment : v_select_peptide_assignment) {
			fpf_filesystem::sample_analysis temp_sample_analysis{};
			if (itr_v_select_spectra_assignment == "PEAKS database match") {
				temp_sample_analysis.PEAKS_database_exists = fpf_parse::check_csv_PEAKS_database_peptides_empty(main_v_csv_PEAKS_database_peptides, filesystem_modified);
				temp_sample_analysis.peptide_assignment_method = "PEAKS_database";
				file_found = temp_sample_analysis.PEAKS_database_exists;
			}
			if (itr_v_select_spectra_assignment == "PEAKS de novo") {
				temp_sample_analysis.PEAKS_denovo_exists = fpf_parse::check_csv_PEAKS_denovo_peptides_empty(main_v_csv_PEAKS_denovo_peptides, filesystem_modified);
				temp_sample_analysis.peptide_assignment_method = "PEAKS_denono";
				file_found = temp_sample_analysis.PEAKS_denovo_exists;
			}
			if (itr_v_select_spectra_assignment == "NOVOR de novo") {
				temp_sample_analysis.NOVOR_denovo_exists = fpf_parse::check_csv_NOVOR_denovo_peptides_empty(main_v_csv_NOVOR_denovo_peptides, filesystem_modified);
				temp_sample_analysis.peptide_assignment_method = "NOVOR_denono";
				file_found = temp_sample_analysis.NOVOR_denovo_exists;
			}
			itr_v_filesystem.v_sample_analysis.push_back(std::move(temp_sample_analysis));
		}

		if (file_found) {
			for (auto& itr_v_sample_analysis : itr_v_filesystem.v_sample_analysis) {
				std::ifstream fin_FASTA(select_FASTA);
				fpf_parse::custom_FASTA_output(select_FASTA);
				std::cout << "\n\n\n * parsing FASTA file... \n";
				vector<fpf_parse::FASTA_data> main_FASTA = fpf_parse::parse_FASTA(fin_FASTA);

				if (check_FASTA_file_empty(main_FASTA)) { return 1; };
				fpf_parse::output_v_FASTA_data(main_FASTA);
				fpf_parse::output_v_FASTA_data_to_blastdirectory(main_FASTA);

				vector<fpf_parse::csv_data> main_v_csv_peptides{};
				if (itr_v_sample_analysis.PEAKS_database_exists) {
					std::cout << "\n\n\n * parsing " << itr_v_filesystem.filename << " PEAKS database matched peptides...";
					main_v_csv_peptides = std::move(main_v_csv_PEAKS_database_peptides);
				}
				if (itr_v_sample_analysis.PEAKS_denovo_exists) {
					std::cout << "\n\n\n * parsing " << itr_v_filesystem.filename << " PEAKS de novo peptides...";
					main_v_csv_peptides = std::move(main_v_csv_PEAKS_denovo_peptides);
					if (!IgFamily::FILESYSTEM_MODE) {
						itr_v_filesystem.filename = main_v_csv_peptides.begin()->csv_sourcefile;
						std::cout << main_v_csv_peptides.begin()->csv_sourcefile;;
					}
				}
				if (itr_v_sample_analysis.NOVOR_denovo_exists) {
					std::cout << "\n\n\n * parsing " << itr_v_filesystem.filename << " NOVOR de novo matched peptides...";
					main_v_csv_peptides = std::move(main_v_csv_NOVOR_denovo_peptides);
				}

				std::cout << "\n\n peptides parsed - " << main_v_csv_peptides.size();
				std::cout << "\n\n\n\n creating data structures for file ";
				std::cout << itr_v_filesystem.filename;
				std::cout << "...";
				itr_v_sample_analysis.v_peptide_data = fpf_data::create_v_peptide_data(main_v_csv_peptides);
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
				if (itr_v_sample_analysis.PEAKS_database_exists) {
					fpf_core::core_report(itr_v_filesystem, itr_v_sample_analysis);
				}
				if (itr_v_sample_analysis.PEAKS_denovo_exists) {
					fpf_core::core_report(itr_v_filesystem, itr_v_sample_analysis);
				}
				if (itr_v_sample_analysis.NOVOR_denovo_exists) {
					fpf_core::core_report(itr_v_filesystem, itr_v_sample_analysis);
				}
				fpf_filesystem::fout_filesystem(itr_v_filesystem);
			}
		}
		else {
			std::cout << "\n ...no data file found";
		}
	}

	string farewell;
	std::cout << "\n\n\n program complete...";
	std::cout << "\n\n\n input any key to exit...\n\n > ";
	std::cin >> farewell;

	return 0;
}
