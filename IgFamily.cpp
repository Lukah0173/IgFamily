// * * IgFamily v0.7.5 * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#include <string> // provides - string
#include <iostream> // provides - std::cin, std::cout
#include <iomanip> // provides - std::setprecision
#include <fstream> // provides - std::ifstream

#include "IgFamily.h"
#include "fpf_filesystem.h"
#include "fpf_filesystem_analysis.h"
#include "fpf_data_analysis.h"
#include "fpf_blastp_analysis.h"
#include "fpf_dirichlet_mixture_model.h"
#include "fpf_report.h"
#include "fpf_multinomial.h"


int main() {

	using std::string;
	using std::vector;

	typedef std::string string;

	std::cout << "-- IgFamily " << IgFamily::version << " --\n\n\n";

	fpf_filesystem::display_default();

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

	vector<fpf_data::FASTA_category> main_v_FASTA_category = fpf_data::create_FASTA_category(main_FASTA);

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

		if (itr_v_filesystem.proteinpeptides_exist) {
			std::cout << "\n\n * parsing " << main_v_csv_proteinpeptides.begin()->csv_file;
			std::cout << "\n\n ~ peptides parsed - " << main_v_csv_proteinpeptides.size();
			std::cout << "\n\n\n creating data structures...";
			vector<fpf_data::peptide_data> main_v_peptide_data = fpf_data::create_peptide_data(main_v_csv_proteinpeptides);
			std::cout << "\n\n ...data structures assigned";
		}

		if (itr_v_filesystem.denovopeptides_exist) {
			std::cout << "\n\n * parsing " << main_v_csv_denovopeptides.begin()->csv_file;
			std::cout << "\n\n ~ peptides parsed - " << main_v_csv_denovopeptides.size();
			std::cout << "\n\n creating data structures...";
			vector<fpf_data::peptide_data> main_v_peptide_data = fpf_data::create_peptide_data(main_v_csv_denovopeptides);
			std::cout << "\n\n ...data structures assigned";

			itr_v_filesystem.v_peptide_data = main_v_peptide_data;
			itr_v_filesystem.v_FASTA_category = main_v_FASTA_category;
		}
	}

	for (auto& itr_v_filesystem : v_filesystem) {
		if (itr_v_filesystem.denovopeptides_exist) {
			std::cout << "\n\n\n\n analysing homology for file ";
			std::cout << itr_v_filesystem.filename;
			std::cout << "...";
			fpf_blastp_analysis::create_blastp_input(itr_v_filesystem);
			fpf_blastp_analysis::create_blastp_database(itr_v_filesystem);
			fpf_blastp_analysis::sys_blastp(itr_v_filesystem);
			std::cout << "\n\n\n ...homology analysis complete";
			std::cout << "\n\n\n creating homology data structures for file ";
			std::cout << itr_v_filesystem.filename;
			std::cout << "...";
			fpf_blastp_analysis::create_v_blastp_data(itr_v_filesystem);
			fpf_blastp_analysis::modify_filesystem_blastp_data(itr_v_filesystem);
			fpf_blastp_analysis::associate_blastp_data_to_v_FASTA_category(itr_v_filesystem);
			fpf_blastp_analysis::associate_blastp_data_to_v_peptide_data(itr_v_filesystem);
			fpf_blastp_analysis::create_query_alignment(itr_v_filesystem);		
			fpf_blastp_analysis::normalise_blastp_data(itr_v_filesystem);
			fpf_blastp_analysis::determine_blastp_parameter_density(itr_v_filesystem);
			std::cout << "\n\n ...data structures assigned";
			std::cout << "\n\n outputting homology summary...";
			fpf_blastp_analysis::fout_blastp_summary(itr_v_filesystem);
			std::cout << "\n\n ...homology file ";
			std::cout << itr_v_filesystem.filename;
			std::cout << " output";
		}
	}

	std::cout << "\n\n\n\n scoring categories...\n";
	for (auto& itr_v_filesystem : v_filesystem) {
		if (itr_v_filesystem.denovopeptides_exist) {
			fpf_data_analysis::create_category_analysis(itr_v_filesystem);
			fpf_data_analysis::create_proteinconstruct_from_denovo(itr_v_filesystem);
			for (auto& itr_v_category_analysis : itr_v_filesystem.v_category_analysis) {
				fpf_data_analysis::sort_v_blastp_data_with_spectralcount(itr_v_category_analysis.v_blastp_data_combined_by_category);
			}
			fpf_data_analysis::sort_v_category_analysis(itr_v_filesystem.v_category_analysis);
		}
	}

	for (auto& itr_v_filesystem : v_filesystem) {
		if (itr_v_filesystem.denovopeptides_exist) {
			std::cout << "\n\n\n\n determining most-probable germline representation...\n";
			fpf_data_analysis::select_category_analysis_by_score(itr_v_filesystem);
			std::cout << "\n\n\n analysing post-hoc homology for file ";
			std::cout << itr_v_filesystem.filename;
			std::cout << "...";
			fpf_blastp_analysis::create_blastp_database_refined(itr_v_filesystem);
			fpf_blastp_analysis::sys_blastp(itr_v_filesystem);
			std::cout << "\n\n\n ...homology analysis complete";
			std::cout << "\n\n\n creating homology data structures for file ";
			std::cout << itr_v_filesystem.filename;
			std::cout << "...";
			fpf_blastp_analysis::create_v_blastp_data(itr_v_filesystem);
			fpf_blastp_analysis::modify_filesystem_blastp_data(itr_v_filesystem);
			fpf_blastp_analysis::associate_blastp_data_to_v_FASTA_category(itr_v_filesystem);
			fpf_blastp_analysis::associate_blastp_data_to_v_peptide_data(itr_v_filesystem);
			fpf_blastp_analysis::create_protein_from_category_analysis(itr_v_filesystem);
			fpf_blastp_analysis::create_query_alignment(itr_v_filesystem);
			fpf_blastp_analysis::normalise_blastp_data(itr_v_filesystem);
			fpf_blastp_analysis::determine_blastp_parameter_density(itr_v_filesystem);
			std::cout << "\n\n ...data structures assigned";
			std::cout << "\n\n outputting homology summary...";
			fpf_blastp_analysis::fout_blastp_summary(itr_v_filesystem);
			std::cout << "\n\n ...homology file ";
			std::cout << itr_v_filesystem.filename;
			std::cout << " output";
		}
	}

	std::cout << "\n\n\n\n creating multinomial data frames...\n";
	for (auto& itr_v_filesystem : v_filesystem) {
		if (itr_v_filesystem.denovopeptides_exist) {
			fpf_multinomial::create_filesystem_multinomial_data(itr_v_filesystem);
			fpf_multinomial::fout_multinomial(itr_v_filesystem);
			fpf_multinomial::fout_multinomial_element(itr_v_filesystem);
			fpf_multinomial::fout_multinomial_element_nomatch(itr_v_filesystem);
		}
	}

	std::cout << "\n\n\n\n producing summary reports...\n";
	for (auto& itr_v_filesystem : v_filesystem) {
		if (itr_v_filesystem.denovopeptides_exist) {
			std::cout << "\n\n ...generating multinomial report for " << itr_v_filesystem.filename;
			fpf_report::fout_multinomial_comparison(itr_v_filesystem);
			std::cout << "\n\n ...generating html report for " << itr_v_filesystem.filename;
			fpf_report::fout_html_report(itr_v_filesystem);
			fpf_report::fout_html_report_filtered(itr_v_filesystem);
		}
	}

	for (auto itr_v_filesystem : v_filesystem) {
		fpf_filesystem::fout_filesystem(itr_v_filesystem);
	}



	string farewell;
	std::cout << "\n\n\n\n program complete...";
	std::cout << "\n\n\n\n input any key to exit...\n\n > ";
	std::cin >> farewell;

	return 0;
}
