// * * IgFamily v0.7.13 * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <utility>

#include "IgFamily.h"
#include "fpf_core.h"
#include "fpf_data_analysis.h"
#include "fpf_dirichlet_mixture_model.h"
#include "fpf_filesystem.h"
#include "fpf_filesystem_analysis.h"
#include "fpf_homology_analysis.h"
#include "fpf_multinomial.h"
#include "fpf_report.h"


int main() {

	using std::string;
	using std::vector;

	std::cout << "-- IgFamily " << IgFamily::version << " --\n\n\n";

	string select_FASTA{ IgFamily::DEFAULT_INPUT_FASTA };
	string select_peptide_assignment{ IgFamily::DEFAULT_PEPTIDE_ASSIGNMENT_METHOD };

	fpf_filesystem::display_settings(select_FASTA, select_peptide_assignment);
	string menu_selection{ fpf_filesystem::display_menu() };

	bool menu_continue{};
	while (!menu_continue) {
		if (menu_selection == "F") {
			select_FASTA = fpf_filesystem::display_FASTA_menu(select_FASTA);	
			fpf_filesystem::display_settings(select_FASTA, select_peptide_assignment);
			menu_selection = fpf_filesystem::display_menu();
		}
		if (menu_selection == "P") {
			select_peptide_assignment = fpf_filesystem::display_peptide_assignment_menu(select_peptide_assignment);
			fpf_filesystem::display_settings(select_FASTA, select_peptide_assignment);
			menu_selection = fpf_filesystem::display_menu();
		}
		if (menu_selection == "X") {
			menu_continue = true;
		}
	}

	std::cout << "\n\n\n reading root directory...\n";

	vector<string> v_root_directory = fpf_filesystem::read_root_dir(IgFamily::IGFAMILY_ROOT_DIR);
	vector<fpf_filesystem::filesystem> v_filesystem = fpf_filesystem::read_filesystem(v_root_directory);
	for (auto itr_v_filesystem : v_filesystem) {
		if (itr_v_filesystem.fileconversion) {
			//fpf_filesystem::perform_fileconversion(itr_v_filesystem);
		}
	}

	for (auto& itr_v_filesystem : v_filesystem) {
		std::ifstream fin_FASTA(IgFamily::DEFAULT_INPUT_FASTA_DIRECTORY);
		fpf_parse::custom_FASTA_output(IgFamily::DEFAULT_INPUT_FASTA_DIRECTORY);
		std::cout << "\n\n\n * parsing FASTA file... \n\n";
		vector<fpf_parse::FASTA_data> main_FASTA = fpf_parse::parse_FASTA(fin_FASTA);

		if (check_FASTA_file_empty(main_FASTA)) { return 1; };
		fpf_parse::output_v_FASTA_data(main_FASTA);
		fpf_parse::output_v_FASTA_data_to_blastdirectory(main_FASTA);

		bool filesystem_modified = bool();
		vector<fpf_parse::csv_data> main_v_csv_PEAKS_database_peptides;
		vector<fpf_parse::csv_data> main_v_csv_PEAKS_denovo_peptides;
		vector<fpf_parse::csv_data> main_v_csv_NOVOR_denovo_peptides;
		if ((IgFamily::FILESYSTEM_UPDATE_ALL) || (itr_v_filesystem.fileversion != IgFamily::version)) {
			std::cout << "\n\n\n\n parsing data...";
			main_v_csv_PEAKS_database_peptides = fpf_filesystem::parse_filesystem_PEAKS_database_peptides(fpf_filesystem::read_filesystem_PEAKS_database_peptides(itr_v_filesystem.directory));
			main_v_csv_PEAKS_denovo_peptides = fpf_filesystem::parse_filesystem_PEAKS_denovo_peptides(fpf_filesystem::read_filesystem_PEAKS_denovo_peptides(itr_v_filesystem.directory));
			main_v_csv_NOVOR_denovo_peptides = fpf_filesystem::parse_filesystem_NOVOR_denovo_peptides(fpf_filesystem::read_filesystem_NOVOR_denovo_peptides(itr_v_filesystem.directory));
			filesystem_modified = true;
		}

		itr_v_filesystem.PEAKS_database_exists = fpf_parse::check_csv_PEAKS_database_peptides_empty(main_v_csv_PEAKS_database_peptides, filesystem_modified);
		itr_v_filesystem.PEAKS_denovo_exists = fpf_parse::check_csv_PEAKS_denovo_peptides_empty(main_v_csv_PEAKS_denovo_peptides, filesystem_modified);
		itr_v_filesystem.NOVOR_denovo_exists = fpf_parse::check_csv_NOVOR_denovo_peptides_empty(main_v_csv_NOVOR_denovo_peptides, filesystem_modified);

		if (itr_v_filesystem.PEAKS_database_exists) {
			std::cout << "\n\n * parsing " << itr_v_filesystem.filename << " PEAKS database matched peptides...";
			std::cout << "\n\n ~ peptides parsed - " << main_v_csv_PEAKS_database_peptides.size();
			std::cout << "\n\n creating data structures...";
			fpf_filesystem::sample_analysis main_sample_analysis{};
			main_sample_analysis.v_peptide_data = fpf_data::create_peptide_data(main_v_csv_PEAKS_database_peptides);
			main_sample_analysis.v_peptide_analysis = fpf_data::create_peptide_analysis(main_sample_analysis.v_peptide_data);
			main_sample_analysis.v_protein_data = fpf_data::create_protein_data(main_FASTA);
			itr_v_filesystem.sample_PEAKS_database_analysis = main_sample_analysis;
			std::cout << "\n\n ...data structures assigned";

			fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_database_analysis, false);
			fpf_core::core_data_analysis(itr_v_filesystem.sample_PEAKS_database_analysis);
			fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_database_analysis, true);
			fpf_core::core_data_analysis(itr_v_filesystem.sample_PEAKS_database_analysis);
			fpf_core::core_multinomial(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_database_analysis);
			fpf_core::core_report(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_database_analysis, "PEAKS_database");
			fpf_filesystem::fout_filesystem(itr_v_filesystem);
		}

		if (itr_v_filesystem.PEAKS_denovo_exists) {
			std::cout << "\n\n * parsing " << itr_v_filesystem.filename << " PEAKS de novo matched peptides...";
			std::cout << "\n\n ~ peptides parsed - " << main_v_csv_PEAKS_denovo_peptides.size();
			std::cout << "\n\n creating data structures...";
			fpf_filesystem::sample_analysis main_sample_analysis{};
			main_sample_analysis.v_peptide_data = fpf_data::create_peptide_data(main_v_csv_PEAKS_denovo_peptides);
			main_sample_analysis.v_peptide_analysis = fpf_data::create_peptide_analysis(main_sample_analysis.v_peptide_data);
			main_sample_analysis.v_protein_data = fpf_data::create_protein_data(main_FASTA);
			itr_v_filesystem.sample_PEAKS_denovo_analysis = main_sample_analysis;
			std::cout << "\n\n ...data structures assigned";

			fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_denovo_analysis, false);
			fpf_core::core_data_analysis(itr_v_filesystem.sample_PEAKS_denovo_analysis);
			fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_denovo_analysis, true);
			fpf_core::core_data_analysis(itr_v_filesystem.sample_PEAKS_denovo_analysis);
			fpf_core::core_multinomial(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_denovo_analysis);
			fpf_core::core_report(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_denovo_analysis, "PEAKS_denovo");
			fpf_filesystem::fout_filesystem(itr_v_filesystem);
		}

		if (itr_v_filesystem.NOVOR_denovo_exists) {
			std::cout << "\n\n * parsing " << itr_v_filesystem.filename << " NOVOR de novo matched peptides...";
			std::cout << "\n\n ~ peptides parsed - " << main_v_csv_NOVOR_denovo_peptides.size();
			std::cout << "\n\n creating data structures...";
			fpf_filesystem::sample_analysis main_sample_analysis{};
			main_sample_analysis.v_peptide_data = fpf_data::create_peptide_data(main_v_csv_NOVOR_denovo_peptides);
			main_sample_analysis.v_peptide_analysis = fpf_data::create_peptide_analysis(main_sample_analysis.v_peptide_data);
			main_sample_analysis.v_protein_data = fpf_data::create_protein_data(main_FASTA);
			itr_v_filesystem.sample_NOVOR_denovo_analysis = main_sample_analysis;
			std::cout << "\n\n ...data structures assigned";

			fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_filesystem.sample_NOVOR_denovo_analysis, false);
			fpf_core::core_data_analysis(itr_v_filesystem.sample_NOVOR_denovo_analysis);
			fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_filesystem.sample_NOVOR_denovo_analysis, true);
			fpf_core::core_data_analysis(itr_v_filesystem.sample_NOVOR_denovo_analysis);
			fpf_core::core_multinomial(itr_v_filesystem, itr_v_filesystem.sample_NOVOR_denovo_analysis);
			fpf_core::core_report(itr_v_filesystem, itr_v_filesystem.sample_NOVOR_denovo_analysis, "NOVOR_denovo");
			fpf_filesystem::fout_filesystem(itr_v_filesystem);
		}
	}



	string farewell;
	std::cout << "\n\n\n\n program complete...";
	std::cout << "\n\n\n\n input any key to exit...\n\n > ";
	std::cin >> farewell;

	return 0;
}
