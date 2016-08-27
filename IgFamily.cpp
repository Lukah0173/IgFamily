// * * IgFamily v0.7.11 * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#include <string>			// provides - string
#include <iostream>			// provides - std::cin, std::cout
#include <iomanip>			// provides - std::setprecision
#include <fstream>			// provides - std::ifstream
#include <utility>			// provides - std::move

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

	typedef std::string string;

	std::cout << "-- IgFamily " << IgFamily::version << " --\n\n\n";

	//fpf_filesystem::display_default();

	int catch_menu_selection;
	catch_menu_selection = fpf_filesystem::perform_menu_selection(fpf_filesystem::display_menu());
	if (catch_menu_selection == 1) {
	}

	std::ifstream fin_FASTA(IgFamily::DEFAULT_INPUT_FASTA_DIRECTORY);

	fpf_parse::custom_FASTA_output(IgFamily::DEFAULT_INPUT_FASTA_DIRECTORY);

	std::cout << "\n\n\n * parsing FASTA file... \n\n";

	vector<fpf_parse::FASTA_data> main_FASTA = fpf_parse::parse_FASTA(fin_FASTA);
	if (check_FASTA_file_empty(main_FASTA)) { return 1; };
	fpf_parse::output_v_FASTA_data(main_FASTA);
	fpf_parse::output_v_FASTA_data_to_blastdirectory(main_FASTA);

	vector<fpf_parse::csv_data> main_v_csv_proteinpeptides;
	vector<fpf_parse::csv_data> main_v_csv_denovopeptides;

	std::cout << "\n\n\n reading root directory...\n";

	vector<string> v_root_directory = fpf_filesystem::read_root_dir(IgFamily::IGFAMILY_ROOT_DIR);
	vector<fpf_filesystem::filesystem> v_filesystem = fpf_filesystem::read_filesystem(v_root_directory);
	for (auto itr_v_filesystem : v_filesystem) {
		if (itr_v_filesystem.fileconversion) {
			//fpf_filesystem::perform_fileconversion(itr_v_filesystem);
		}
	}

	for (auto& itr_v_filesystem : v_filesystem) {
		bool filesystem_modified = bool();
		if ((IgFamily::FILESYSTEM_UPDATE_ALL) || (itr_v_filesystem.fileversion != IgFamily::version)) {
			std::cout << "\n\n\n\n parsing data...";
			//main_v_csv_proteinpeptides = fpf_filesystem::parse_filesystem_proteinpeptides(fpf_filesystem::read_filesystem_proteinpeptides(itr_v_filesystem.directory));
			main_v_csv_denovopeptides = fpf_filesystem::parse_filesystem_denovopeptides(fpf_filesystem::read_filesystem_denovopeptides(itr_v_filesystem.directory));
			filesystem_modified = true;
		}
		//itr_v_filesystem.proteinpeptides_exist = fpf_parse::check_csv_proteinpeptides_empty(main_v_csv_proteinpeptides, filesystem_modified);
		itr_v_filesystem.denovopeptides_exist = fpf_parse::check_csv_denovopeptides_empty(main_v_csv_denovopeptides, filesystem_modified);

		//if (fpf_filesystem::check_filesystem_current(main_v_csv_proteinpeptides, main_v_csv_denovopeptides, filesystem_modified)) {
		//	return EXIT_FILESYTEM_CURRENT;
		//}

		if (itr_v_filesystem.denovopeptides_exist) {
			std::cout << "\n\n * parsing " << itr_v_filesystem.filename;
			std::cout << "\n\n ~ peptides parsed - " << main_v_csv_denovopeptides.size();
			std::cout << "\n\n creating data structures...";
			vector<fpf_data::peptide_data> main_v_peptide_data_PEAKS_denovo = fpf_data::create_peptide_data(main_v_csv_denovopeptides);
			vector<fpf_data::peptide_analysis> main_v_peptide_analysis_PEAKS_denovo = fpf_data::create_peptide_analysis(main_v_peptide_data_PEAKS_denovo);
			std::cout << "\n\n ...data structures assigned";

			static vector<fpf_data::peptide_data> main_v_peptide_data = std::move(main_v_peptide_data_PEAKS_denovo);
			static vector<fpf_data::peptide_analysis> main_v_peptide_analysis = std::move(main_v_peptide_analysis_PEAKS_denovo);
			static vector<fpf_data::protein_data> main_v_protein_data = fpf_data::create_protein_data(main_FASTA);
			fpf_filesystem::sample_analysis main_sample_analysis{};
			main_sample_analysis.v_peptide_data = main_v_peptide_data;
			main_sample_analysis.v_peptide_analysis = main_v_peptide_analysis;
			main_sample_analysis.v_protein_data = main_v_protein_data;
			itr_v_filesystem.sample_PEAKS_denovo_analysis = main_sample_analysis;
		}
	}

	for (auto& itr_v_filesystem : v_filesystem) {
		fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_denovo_analysis, false);
		fpf_core::core_data_analysis(itr_v_filesystem.sample_PEAKS_denovo_analysis);
		fpf_core::core_homology_analysis(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_denovo_analysis, true);
		fpf_core::core_data_analysis(itr_v_filesystem.sample_PEAKS_denovo_analysis);
		fpf_core::core_multinomial(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_denovo_analysis);
		fpf_core::core_report(itr_v_filesystem, itr_v_filesystem.sample_PEAKS_denovo_analysis);
		fpf_filesystem::fout_filesystem(itr_v_filesystem);
	}



	string farewell;
	std::cout << "\n\n\n\n program complete...";
	std::cout << "\n\n\n\n input any key to exit...\n\n > ";
	std::cin >> farewell;

	return 0;
}
