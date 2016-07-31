// * * IgFamily v0.5.8 * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *

#include <string> // provides - std::string
#include <iostream> // provides - std::cin, std::cout
#include <iomanip> // provides - std::setprecision
#include <fstream> // provides - std::ifstream

#include "IgFamily.h"
#include "fpf_filesystem.h"
#include "fpf_filesystem_analysis.h"
#include "fpf_blastp_analysis.h"
#include "fpf_dirichlet_mixture_model.h"
#include "fpf_report.h"



int main() {

	typedef size_t size_type;
	typedef std::string string_type;
	typedef char char_type;

	std::cout << "-- IgFamily " << IgFamily::version << " --\n\n\n";

	if (IgFamily::FILESYSTEM_MODE) {
		std::cout << "Enter any key to continue: \n\n > ";
		string_type filesystem_initial;
		std::cin >> filesystem_initial;
	}
	else {
		std::cout << "Enter [filename].csv of a comma separated variables file from PEAKS output: \n\n > ";
		std::cin >> IgFamily::INPUT_CSV;
		std::string output_holder = IgFamily::INPUT_CSV;
		std::string con_output = std::string();
		con_output = IgFamily::INPUT_CSV;
		IgFamily::INPUT_CSV = IgFamily::INPUT_CSV + ".csv";
	}

	std::ifstream fin_input_csv(IgFamily::INPUT_CSV);
	std::ifstream fin_input_FASTA(IgFamily::INPUT_FASTA);
	std::vector<fpf_parse::parse_FASTA_type> main_v_c_parse_FASTA = fpf_parse::parse_FASTA(fin_input_FASTA);
	if (b_parse_FASTA_empty(main_v_c_parse_FASTA)) { return 1; };

	std::vector<fpf_parse::parse_peptides_csv_type> main_v_c_parse_csv_proteinpeptides_data;
	std::vector<fpf_parse::parse_peptides_csv_type> main_v_c_parse_csv_denovopeptides_data;
	std::vector<string_type> v_str_root_dir = fpf_filesystem::read_root_dir(IgFamily::IGFAMILY_ROOT_DIR);
	std::vector<fpf_filesystem::filesystem_type> v_s_filesystem = fpf_filesystem::read_filesystem(v_str_root_dir);

	std::cout << "\n\n\n * parsing FASTA file";

	std::vector<fpf_data::multinomial_element_data_type> main_v_s_multinomial_element_data = fpf_data::create_v_s_multinomial_element_data(main_v_c_parse_FASTA);
	std::vector<fpf_data::multinomial_element_data_type> main_v_s_multinomial_element_data_distinct = fpf_data::create_v_s_multinomial_element_data_distinct(main_v_s_multinomial_element_data);

	if (!IgFamily::FILESYSTEM_MODE) {
		fpf_filesystem::filesystem_type one_iteration;
		v_s_filesystem.push_back(one_iteration);
	}

	for (auto& itr_v_s_filesystem : v_s_filesystem) {
		bool filesystem_modified = bool();
		if (IgFamily::FILESYSTEM_MODE) {
			if ((IgFamily::FILESYSTEM_UPDATE_ALL) || (itr_v_s_filesystem.str_fileversion != IgFamily::version)) {
				std::cout << "\n\n\n\n ...parsing data";
				main_v_c_parse_csv_proteinpeptides_data = fpf_filesystem::parse_filesystem_proteinpeptides(fpf_filesystem::read_filesystem_proteinpeptides(itr_v_s_filesystem.str_directory));
				main_v_c_parse_csv_denovopeptides_data = fpf_filesystem::parse_filesystem_denovopeptides(fpf_filesystem::read_filesystem_denovopeptides(itr_v_s_filesystem.str_directory));
				filesystem_modified = true;
			}
		}
		else {
			std::cout << "\n\n\n\n ...parsing data";
			main_v_c_parse_csv_proteinpeptides_data = fpf_parse::parse_proteinpeptides(fin_input_csv, IgFamily::INPUT_CSV);
		}
		itr_v_s_filesystem.b_proteinpeptides_exist = fpf_parse::check_protein_peptides(main_v_c_parse_csv_proteinpeptides_data, filesystem_modified);
		itr_v_s_filesystem.b_denovopeptides_exist = fpf_parse::check_denovo_peptides(main_v_c_parse_csv_denovopeptides_data, filesystem_modified);

		//if (fpf_filesystem::check_filesystem_current(main_v_c_parse_csv_proteinpeptides_data, main_v_c_parse_csv_denovopeptides_data, filesystem_modified)) {
		//	return EXIT_FILESYTEM_CURRENT;
		//}

		if (itr_v_s_filesystem.b_proteinpeptides_exist) {
			std::cout << "\n\n  * parsing " << main_v_c_parse_csv_proteinpeptides_data.begin()->str_parse_peptides_csv_file;
			std::cout << "\n\n peptides parsed - " << main_v_c_parse_csv_proteinpeptides_data.size();
			std::cout << "\n\n ...creating data structures";
			std::vector<fpf_data::peptide_data_type> main_v_c_proteinpeptides_data = fpf_data::create_v_s_peptide_data(main_v_c_parse_csv_proteinpeptides_data);
			std::vector<fpf_data::peptide_data_type> main_v_c_proteinpeptides_data_filtered = fpf_data::create_v_s_peptide_data_filtered(main_v_c_proteinpeptides_data);
			std::vector<fpf_data::peptide_data_type> main_v_c_proteinpeptides_data_distinct = fpf_data::create_v_s_peptide_data_distinct(main_v_c_proteinpeptides_data);
			std::vector<fpf_data::peptide_data_type> main_v_c_proteinpeptides_data_filtered_distinct = fpf_data::create_v_s_peptide_data_filtered_distinct(main_v_c_proteinpeptides_data_filtered);

			if (SIMPLE_SCORE) {
				create_global_score_mean(main_v_s_multinomial_element_data);
				std::cout << "\n\n ...assigning peptides to gene families";
				fpf_data::create_v_s_multinomial_element_data_v_s_peptide_data(main_v_s_multinomial_element_data, main_v_c_proteinpeptides_data);
				fpf_data::create_v_s_multinomial_element_data_v_s_peptide_data_distinct_filtered(main_v_s_multinomial_element_data, main_v_c_proteinpeptides_data_filtered_distinct);
				fpf_data::create_v_s_multinomial_element_data_str_alignment(main_v_s_multinomial_element_data);
				std::cout << "\n\n ...determining sequence coverage and total spectral count";
				fpf_data::create_v_s_multinomial_element_data_st_totalspectralcount(main_v_s_multinomial_element_data);
				fpf_data::create_v_s_multinomial_element_data_d_coverage(main_v_s_multinomial_element_data);
				std::cout << "\n\n ...determining peptide association co-occurence";
				fpf_data::create_v_c_peptide_v_p_peptideassociation(main_v_s_multinomial_element_data, main_v_c_proteinpeptides_data);
				fpf_data::create_v_c_peptide_v_str_peptideassociation_distinct(main_v_c_proteinpeptides_data);
				fpf_data::create_v_c_peptide_v_p_peptideassociation_distinct(main_v_s_multinomial_element_data_distinct, main_v_c_proteinpeptides_data);
				std::cout << "\n\n ...calculating score";
				fpf_data::create_v_s_multinomial_element_data_d_score(main_v_c_proteinpeptides_data, main_v_s_multinomial_element_data, main_v_s_multinomial_element_data_distinct);
				std::cout << "\n\n ...training score";
				fpf_data::train_v_s_multinomial_element_data_d_score(main_v_c_proteinpeptides_data, main_v_s_multinomial_element_data, main_v_s_multinomial_element_data_distinct);
				std::cout << "\n\n ...formatting output";
				fpf_data::sort_v_s_peptide_data_str_peptide(main_v_c_proteinpeptides_data);
				fpf_data::sort_v_s_peptide_data_str_peptide(main_v_c_proteinpeptides_data_filtered_distinct);
				fpf_data::update_v_s_multinomial_element_distinctpolymorphism_data(main_v_s_multinomial_element_data_distinct);

				std::vector<fpf_data::multinomial_element_data_type*> map_main_v_s_multinomial_element_data = map_v_s_multinomial_element_data_by_score(main_v_s_multinomial_element_data);
				std::vector<fpf_data::multinomial_element_data_type*> map_main_v_s_multinomial_element_data_distict = map_v_s_multinomial_element_data_by_score(main_v_s_multinomial_element_data_distinct);

				bool b_map_main_v_s_multinomial_element_data_distict = true;
				std::vector<fpf_data::multinomial_element_data_type*> con_map_main_v_s_multinomial_element_data;
				if (b_map_main_v_s_multinomial_element_data_distict) {
					con_map_main_v_s_multinomial_element_data = map_main_v_s_multinomial_element_data_distict;
				}
				else {
					con_map_main_v_s_multinomial_element_data = map_main_v_s_multinomial_element_data;
				}

				fpf_data::output(main_v_c_parse_csv_proteinpeptides_data.begin()->str_parse_peptides_csv_file, map_main_v_s_multinomial_element_data_distict, con_map_main_v_s_multinomial_element_data);
				std::cout << "\n\n  * file " << main_v_c_parse_csv_proteinpeptides_data.begin()->str_parse_peptides_csv_file + "_output.txt" << " output";

				IgFamily::GLOBAL_ITERATOR = size_type();
			}
		}

		if (itr_v_s_filesystem.b_denovopeptides_exist) {
			std::cout << "\n\n  * parsing " << main_v_c_parse_csv_denovopeptides_data.begin()->str_parse_peptides_csv_file;
			std::cout << "\n\n peptides parsed - " << main_v_c_parse_csv_denovopeptides_data.size();
			std::cout << "\n\n ...creating data structures";
			std::vector<fpf_data::peptide_data_type> main_v_c_denovopeptides_data = fpf_data::create_v_s_peptide_data(main_v_c_parse_csv_denovopeptides_data);
			std::vector<fpf_data::peptide_data_type> main_v_c_denovopeptides_data_filtered = fpf_data::create_v_s_peptide_data_filtered(main_v_c_denovopeptides_data);
			std::vector<fpf_data::peptide_data_type> main_v_c_denovopeptides_data_distinct = fpf_data::create_v_s_peptide_data_distinct(main_v_c_denovopeptides_data);
			std::vector<fpf_data::peptide_data_type> main_v_c_denovopeptides_data_filtered_distinct = fpf_data::create_v_s_peptide_data_filtered_distinct(main_v_c_denovopeptides_data_filtered);

			if (IgFamily::FILESYSTEM_MODE) {
				itr_v_s_filesystem.v_s_peptide_data = main_v_c_denovopeptides_data;
				itr_v_s_filesystem.v_s_peptide_data_filtered = main_v_c_denovopeptides_data_filtered;
				itr_v_s_filesystem.v_s_peptide_data_distinct = main_v_c_denovopeptides_data_distinct;
				itr_v_s_filesystem.v_s_peptide_data_filtered_distinct = main_v_c_denovopeptides_data_filtered_distinct;
				itr_v_s_filesystem.v_c_analysis_data = main_v_s_multinomial_element_data;
				itr_v_s_filesystem.v_c_analysis_distinct_data = main_v_s_multinomial_element_data_distinct;
			}
		}
	}

	if (IgFamily::FILESYSTEM_MODE == 1) {
		for (auto& itr_v_s_filesystem : v_s_filesystem) {
			if (itr_v_s_filesystem.b_denovopeptides_exist) {
				std::cout << "\n\n\n\n ...analysing homology for file ";
				std::cout << itr_v_s_filesystem.str_filename;
				fpf_blastp_analysis::create_blastp_input(itr_v_s_filesystem);
				fpf_blastp_analysis::create_blastp_database(itr_v_s_filesystem);
				fpf_blastp_analysis::sys_blastp(itr_v_s_filesystem);
				fpf_blastp_analysis::create_v_s_blastp(itr_v_s_filesystem);
				fpf_blastp_analysis::modify_v_s_filesystem_blastp_data(itr_v_s_filesystem);
				fpf_blastp_analysis::create_str_protein(itr_v_s_filesystem);
				fpf_blastp_analysis::create_str_query_alignment(itr_v_s_filesystem);
				fpf_blastp_analysis::normalise_v_s_filesystem_blastp_data(itr_v_s_filesystem);
				fpf_blastp_analysis::fout_blastp_summary(itr_v_s_filesystem, BLASTP_THRESHOLD);
				fpf_blastp_analysis::create_s_filesystem_mnom(itr_v_s_filesystem);
				fpf_blastp_analysis::fout_v_s_mnom(itr_v_s_filesystem);
			}
		}

		//fpf_dirichlet_mixture_model::initialise();

		if (IgFamily::FILESYSTEM_MODE == 1) {
			std::cout << "\n\n\n\n ...producing summary reports -\n";
			for (auto& itr_v_s_filesystem : v_s_filesystem) {
				if (itr_v_s_filesystem.b_denovopeptides_exist) {
					fpf_report::create_s_report(itr_v_s_filesystem);
					fpf_report::create_str_protein_from_denovo(itr_v_s_filesystem);
					for (auto& itr_v_s_report : itr_v_s_filesystem.v_s_report) {
						fpf_report::sort_v_s_blastp(itr_v_s_report.v_s_blastp);
					}
					fpf_report::sort_v_s_report(itr_v_s_filesystem.v_s_report);
					fpf_report::fout_s_report(itr_v_s_filesystem);
					fpf_report::fout_html_report(itr_v_s_filesystem);
				}
			}
		}


		for (auto itr_v_s_filesystem : v_s_filesystem) {
			fpf_filesystem::fout_filesystem(itr_v_s_filesystem);
		}
		//fpf_dirichlet_mixture_model::def_s_model_data s_model_data = fpf_dirichlet_mixture_model::create_s_model_data(v_s_filesystem_blastp_data, main_v_s_multinomial_element_data);
	}

	


	//DWORD WINAPI GetCurrentDirectory(
	//	_In_  DWORD  nBufferLength,
	//	_Out_ LPTSTR lpBuffer
	//	);

	std::string farewell;
	std::cout << "\n\n\n\n ...program complete";
	std::cout << "\n\n\n\n input any key to exit...\n\n > ";
	std::cin >> farewell;

	return 0;
}
