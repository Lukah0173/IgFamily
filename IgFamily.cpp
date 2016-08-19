// * * IgFamily v0.7.0 * * 
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
#include "fpf_blastp_analysis.h"
#include "fpf_dirichlet_mixture_model.h"
#include "fpf_report.h"
#include "fpf_multinomial.h"



int main() {

	using std::string;
	using std::vector;

	typedef size_t size_t;
	typedef string string;
	typedef char char_type;

	std::cout << "-- IgFamily " << IgFamily::version << " --\n\n\n";

	int catch_menu_selection;
	catch_menu_selection = fpf_filesystem::perform_menu_selection(fpf_filesystem::display_menu());
	if (catch_menu_selection == 1) {
	}

	std::ifstream fin_input_csv(IgFamily::INPUT_CSV);
	std::ifstream fin_input_FASTA(IgFamily::INPUT_FASTA);

	fpf_parse::custom_FASTA_output(IgFamily::INPUT_FASTA);

	std::cout << "\n\n\n * parsing FASTA file... \n\n";

	vector<fpf_parse::parse_FASTA_type> main_v_c_parse_FASTA = fpf_parse::parse_FASTA(fin_input_FASTA);
	if (b_parse_FASTA_empty(main_v_c_parse_FASTA)) { return 1; };
	fpf_parse::output_v_c_parse_FASTA(main_v_c_parse_FASTA);
	fpf_parse::output_v_c_parse_FASTA_to_blastdirectory(main_v_c_parse_FASTA);

	vector<fpf_parse::parse_peptides_csv_type> main_v_c_parse_csv_proteinpeptides_data;
	vector<fpf_parse::parse_peptides_csv_type> main_v_c_parse_csv_denovopeptides_data;

	std::cout << "\n\n\n reading root directory...\n";

	vector<string> v_str_root_dir = fpf_filesystem::read_root_dir(IgFamily::IGFAMILY_ROOT_DIR);
	vector<fpf_filesystem::filesystem> v_s_filesystem = fpf_filesystem::read_filesystem(v_str_root_dir);
	for (auto& itr_v_s_filesystem : v_s_filesystem) {
		if (itr_v_s_filesystem.fileconversion) {
			//fpf_filesystem::perform_fileconversion(itr_v_s_filesystem);
		}
	}

	vector<fpf_data::multinomial_category> main_v_s_multinomial_element_data = fpf_data::create_multinomial_category(main_v_c_parse_FASTA);

	for (auto& itr_v_s_filesystem : v_s_filesystem) {
		bool filesystem_modified = bool();
		if ((IgFamily::FILESYSTEM_UPDATE_ALL) || (itr_v_s_filesystem.fileversion != IgFamily::version)) {
			std::cout << "\n\n\n\n parsing data...";
			main_v_c_parse_csv_proteinpeptides_data = fpf_filesystem::parse_filesystem_proteinpeptides(fpf_filesystem::read_filesystem_proteinpeptides(itr_v_s_filesystem.directory));
			main_v_c_parse_csv_denovopeptides_data = fpf_filesystem::parse_filesystem_denovopeptides(fpf_filesystem::read_filesystem_denovopeptides(itr_v_s_filesystem.directory));
			filesystem_modified = true;
		}
		itr_v_s_filesystem.proteinpeptides_exist = fpf_parse::check_protein_peptides(main_v_c_parse_csv_proteinpeptides_data, filesystem_modified);
		itr_v_s_filesystem.denovopeptides_exist = fpf_parse::check_denovo_peptides(main_v_c_parse_csv_denovopeptides_data, filesystem_modified);

		//if (fpf_filesystem::check_filesystem_current(main_v_c_parse_csv_proteinpeptides_data, main_v_c_parse_csv_denovopeptides_data, filesystem_modified)) {
		//	return EXIT_FILESYTEM_CURRENT;
		//}

		if (itr_v_s_filesystem.proteinpeptides_exist) {
			std::cout << "\n\n * parsing " << main_v_c_parse_csv_proteinpeptides_data.begin()->str_parse_peptides_csv_file;
			std::cout << "\n\n ~ peptides parsed - " << main_v_c_parse_csv_proteinpeptides_data.size();
			std::cout << "\n\n\n creating data structures...";
			vector<fpf_data::peptide_data> main_v_c_proteinpeptides_data = fpf_data::create_peptide_data(main_v_c_parse_csv_proteinpeptides_data);
			std::cout << "\n\n ...data structures assigned";
		}

		if (itr_v_s_filesystem.denovopeptides_exist) {
			std::cout << "\n\n * parsing " << main_v_c_parse_csv_denovopeptides_data.begin()->str_parse_peptides_csv_file;
			std::cout << "\n\n ~ peptides parsed - " << main_v_c_parse_csv_denovopeptides_data.size();
			std::cout << "\n\n creating data structures...";
			vector<fpf_data::peptide_data> main_v_c_denovopeptides_data = fpf_data::create_peptide_data(main_v_c_parse_csv_denovopeptides_data);
			std::cout << "\n\n ...data structures assigned";

			itr_v_s_filesystem.v_peptide_data = main_v_c_denovopeptides_data;
			itr_v_s_filesystem.v_multinomial_category = main_v_s_multinomial_element_data;
		}
	}

	for (auto& itr_v_s_filesystem : v_s_filesystem) {
		if (itr_v_s_filesystem.denovopeptides_exist) {
			std::cout << "\n\n\n\n analysing homology for file ";
			std::cout << itr_v_s_filesystem.filename;
			std::cout << "...";
			fpf_blastp_analysis::create_blastp_input(itr_v_s_filesystem);
			fpf_blastp_analysis::create_blastp_database(itr_v_s_filesystem);
			fpf_blastp_analysis::sys_blastp(itr_v_s_filesystem);
			std::cout << "\n\n\n ...homology analysis complete";
			std::cout << "\n\n\n creating homology data structures...";
			fpf_blastp_analysis::create_v_s_blastp(itr_v_s_filesystem);
			fpf_blastp_analysis::modify_v_s_filesystem_blastp_data(itr_v_s_filesystem);
			fpf_blastp_analysis::create_str_protein(itr_v_s_filesystem);
			fpf_blastp_analysis::create_str_query_alignment(itr_v_s_filesystem);		
			fpf_blastp_analysis::normalise_v_s_filesystem_blastp_data(itr_v_s_filesystem);
			std::cout << "\n\n ...data structures assigned";
			std::cout << "\n\n outputting homology summary...";
			fpf_blastp_analysis::fout_blastp_summary(itr_v_s_filesystem);
			std::cout << "\n\n ...homology file ";
			std::cout << itr_v_s_filesystem.filename;
			std::cout << " output";
		}
	}

	//fpf_dirichlet_mixture_model::initialise();

	std::cout << "\n\n\n\n creating multinomial data frames...\n";
	for (auto& itr_v_s_filesystem : v_s_filesystem) {
		if (itr_v_s_filesystem.denovopeptides_exist) {
			fpf_multinomial::create_filesystem_multinomial_data(itr_v_s_filesystem);
			fpf_multinomial::fout_multinomial(itr_v_s_filesystem);
			fpf_multinomial::fout_multinomial_element(itr_v_s_filesystem);
			fpf_multinomial::fout_multinomial_element_nomatch(itr_v_s_filesystem);
		}
	}

	std::cout << "\n\n\n\n producing summary reports...\n";
	for (auto& itr_v_s_filesystem : v_s_filesystem) {
		if (itr_v_s_filesystem.denovopeptides_exist) {
			fpf_report::create_category_report(itr_v_s_filesystem);
			fpf_report::create_proteinconstruct_from_denovo(itr_v_s_filesystem);
			for (auto& itr_v_s_report : itr_v_s_filesystem.v_category_report) {
				//fpf_report::sort_v_blastp_data(itr_v_s_report.v_blastp_data);
				fpf_report::sort_v_blastp_data_with_spectralcount(itr_v_s_report.v_blastp_combined_by_category);
			}
			fpf_report::sort_v_category_report(itr_v_s_filesystem.v_category_report);
			std::cout << "\n\n ...generating html report for " << itr_v_s_filesystem.filename;
			fpf_report::fout_html_report(itr_v_s_filesystem);
		}
	}


	for (auto itr_v_s_filesystem : v_s_filesystem) {
		fpf_filesystem::fout_filesystem(itr_v_s_filesystem);
	}
	//fpf_dirichlet_mixture_model::def_s_model_data s_model_data = fpf_dirichlet_mixture_model::create_s_model_data(v_s_filesystem_blastp_data, main_v_s_multinomial_element_data);


	string farewell;
	std::cout << "\n\n\n\n program complete...";
	std::cout << "\n\n\n\n input any key to exit...\n\n > ";
	std::cin >> farewell;

	return 0;
}
