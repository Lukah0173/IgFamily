// * * IgFamily v0.8.7a * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <utility>

#include "IgFamily.h"
#include "fpf_core.h"
#include "fpf_data_analysis.h"
#include "fpf_denovo.h"
#include "fpf_dirichlet_mixture_model.h"
#include "fpf_filesystem.h"
#include "fpf_filesystem_analysis.h"
#include "fpf_homology_analysis.h"
#include "fpf_multinomial.h"
#include "fpf_report.h"
#include "fpf_utility.h"


int main() {

	using std::string;
	using std::vector;

	vector<fpf_utility::transcript> main_v_transcript = fpf_utility::parse_transcript_data();
	fpf_utility::translate_v_transcript(main_v_transcript);
	fpf_utility::fout_transcript_and_translation(main_v_transcript);

	std::cout << "\n -- IgFamily " << IgFamily::version << " --\n";

	string select_FASTA{ IgFamily::DEFAULT_INPUT_FASTA };
	vector<string> v_select_spectra_assignment{ IgFamily::DEFAULT_PEPTIDE_ASSIGNMENT_METHOD };

	fpf_filesystem::display_settings(select_FASTA, v_select_spectra_assignment);
	string menu_selection{ fpf_filesystem::display_menu() };

	bool menu_continue{};
	while (!menu_continue) {
		if (menu_selection == "F") {
			select_FASTA = fpf_filesystem::display_FASTA_menu(select_FASTA);	
			fpf_filesystem::display_settings(select_FASTA, v_select_spectra_assignment);
			menu_selection = fpf_filesystem::display_menu();
		}
		if (menu_selection == "P") {
			*v_select_spectra_assignment.begin() = fpf_filesystem::display_peptide_assignment_menu(*v_select_spectra_assignment.begin());
			fpf_filesystem::display_settings(select_FASTA, v_select_spectra_assignment);
			menu_selection = fpf_filesystem::display_menu();
		}
		if (menu_selection == "X") {
			select_FASTA = IgFamily::DEFAULT_INPUT_FASTA_DIRECTORY + select_FASTA;
			menu_continue = true;
		}
	}

	std::cout << "\n\n\n reading root directory...\n";
	vector<string> v_root_directory = fpf_filesystem::read_root_dir(IgFamily::IGFAMILY_ROOT_DIR);
	vector<fpf_filesystem::filesystem> v_filesystem = fpf_filesystem::read_filesystem(v_root_directory);
	
	for (auto itr_v_filesystem : v_filesystem) {
		if (itr_v_filesystem.fileconversion) {
			//fpf_filesystem::perform_fileconversion(itr_v_filesystem);
			//fpf_denovo::perform_novor_denovo(itr_v_filesystem);
			//string catchin{};
			//std::cin >> catchin;
		}
	}

	for (auto& itr_v_filesystem : v_filesystem) {
		bool filesystem_modified{};
		fpf_filesystem::sample_analysis temp_sample_analysis{};
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
		for (const auto& itr_v_select_spectra_assignment : v_select_spectra_assignment) {
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
			itr_v_filesystem.v_sample_analysis.push_back(temp_sample_analysis);
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
				}
				if (itr_v_sample_analysis.NOVOR_denovo_exists) {
					std::cout << "\n\n\n * parsing " << itr_v_filesystem.filename << " NOVOR de novo matched peptides...";
					main_v_csv_peptides = std::move(main_v_csv_NOVOR_denovo_peptides);
				}

				std::cout << "\n\n peptides parsed - " << main_v_csv_peptides.size();
				std::cout << "\n\n\n\n creating data structures for file ";
				std::cout << itr_v_filesystem.filename;
				std::cout << "...";
				itr_v_sample_analysis.v_peptide_data = fpf_data::create_peptide_data(main_v_csv_peptides);
				itr_v_sample_analysis.v_peptide_analysis = fpf_data::create_peptide_analysis(itr_v_sample_analysis.v_peptide_data);
				itr_v_sample_analysis.v_protein_data = fpf_data::create_protein_data(main_FASTA);
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
