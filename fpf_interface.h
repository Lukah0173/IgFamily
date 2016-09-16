// * * fpf_interface.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_INTERFACE
#define	FPF_INTERFACE

#include <cstdlib>
#include <string>
#include <vector>

#include "IgFamily.h"
#include "fpf_genome_data.h"


namespace fpf_interface {

	using std::vector;

	void display_settings(string par_FASTA_setting, vector<string> par_v_spectra_assignment_method) {
		std::cout << "\n\n Current settings:";
		std::cout << "\n\n";
		std::cout << " FASTA file -                     ";
		std::cout << par_FASTA_setting;
		std::cout << "\n Spectra assignment method -      ";
		size_t count_itr{};
		for (const auto& itr_v_spectra_assignment_method : par_v_spectra_assignment_method) {
			++count_itr;
			std::cout << itr_v_spectra_assignment_method;
			if (count_itr != par_v_spectra_assignment_method.size()) {
				std::cout << ",\n" + string(34,' ');
			}
		}
	}

	string display_menu() {
		std::cout << "\n\n";
		std::cout << " Select FASTA file:                         [F] ";
		std::cout << "\n Select spectra assignment method:          [P]";
		std::cout << "\n Continue with current settings:            [X]";
		string menu_selection{};
		while ((menu_selection != "F") && (menu_selection != "P") && (menu_selection != "X") && (menu_selection != "!!")) {
			menu_selection.clear();
			std::cout << "\n\n --> ";
			std::cin >> menu_selection;
		}
		return menu_selection;
	}

	string perform_FASTA_selection(string par_FASTA_setting_current);

	string display_FASTA_menu(string par_FASTA_setting_current) {
		std::cout << "\n\n Current setting - " << par_FASTA_setting_current;
		std::cout << "\n";
		std::cout << "\n [0] IGHV_IGLV_IGKV_20160827.fasta";
		std::cout << "\n [1] IGHV_IGLV_IGKV_CONT_20160827.fasta";
		std::cout << "\n [2] IGHV_IGLV_IGKV_CONT_UNIPROT_20160827.fasta";
		std::cout << "\n [3] IGHV_IGLV_IGKV_mABV_CONT_20160827.fasta";
		std::cout << "\n [4] IGHV_IGLV_IGKV_MIGHV_MIGKV_MIGLV_20160827.fasta";
		std::cout << "\n [5] IGHV_IGKV_IGLV_MIGHV_MIGKV_MIGLV_CONT_20160827.fasta";
		std::cout << "\n [6] IGH_IGL_IGK_20160827.fasta";
		std::cout << "\n [7] IGH_IGL_IGK_CONT_20160827.fasta";
		std::cout << "\n [8] IGH_IGL_IGK_CONT_UNIPROT_20160827.fasta";
		std::cout << "\n [9] IGH_IGL_IGK_mABV_CONT_20160827.fasta";
		std::cout << "\n [X] Use current setting";
		return perform_FASTA_selection(par_FASTA_setting_current);
	}

	string perform_FASTA_selection(string par_FASTA_setting_current) {
		string menu_selection{};
		while ((menu_selection != ("0"))
			&& (menu_selection != ("1"))
			&& (menu_selection != ("2"))
			&& (menu_selection != ("3"))
			&& (menu_selection != ("4"))
			&& (menu_selection != ("5"))
			&& (menu_selection != ("6"))
			&& (menu_selection != ("7"))
			&& (menu_selection != ("8"))
			&& (menu_selection != ("9"))
			&& (menu_selection != ("X"))) {
			menu_selection.clear();
			std::cout << "\n\n --> ";
			std::cin >> menu_selection;
		}
		if (menu_selection == "0") {
			return "IGHV_IGLV_IGKV_20160827.fasta";
		}
		if (menu_selection == "1") {
			return "IGHV_IGLV_IGKV_CONT_20160827.fasta";
		}
		if (menu_selection == "2") {
			return "IGHV_IGLV_IGKV_CONT_UNIPROT_20160827.fasta";
		}
		if (menu_selection == "3") {
			return "IGHV_IGLV_IGKV_mABV_CONT_20160827.fasta";
		}
		if (menu_selection == "4") {
			return "IGHV_IGLV_IGKV_MIGHV_MIGKV_MIGLV_20160827.fasta";
		}
		if (menu_selection == "5") {
			return "IGHV_IGKV_IGLV_MIGHV_MIGKV_MIGLV_CONT_20160827.fasta";
		}
		if (menu_selection == "6") {
			return "IGH_IGL_IGK_20160827.fasta";
		}
		if (menu_selection == "7") {
			return "IGH_IGL_IGK_CONT_20160827.fasta";
		}
		if (menu_selection == "8") {
			return "IGH_IGL_IGK_CONT_UNIPROT_20160827.fasta";
		}
		if (menu_selection == "9") {
			return "IGH_IGL_IGK_mABV_CONT_20160827.fasta";
		}
		if (menu_selection == "X") {
			return par_FASTA_setting_current;
		}
		return par_FASTA_setting_current;
	}

	string perform_peptide_method_selection(string par_spectra_method_setting);

	//void display_peptide_assignment_current(vector<string>& par_v_spectra_method_setting) {
	//	std::cout << "\n\n Current setting -      ";
	//	size_t count_itr{};
	//	for (const auto& itr_v_spectra_assignment_method : par_v_spectra_method_setting) {
	//		++count_itr;
	//		std::cout << itr_v_spectra_assignment_method;
	//		if (count_itr != par_v_spectra_method_setting.size()) {
	//			std::cout << " + ";
	//		}
	//	}
	//	std::cout << "\n";
	//	count_itr = size_t();
	//	for (const auto& itr_v_spectra_assignment_method : par_v_spectra_method_setting) {			
	//		std::cout << "\n [" << count_itr << "]";
	//		std::cout << " Modify " << itr_v_spectra_assignment_method;
	//		++count_itr;
	//	}
	//	std::cout << "\n[" << count_itr << "]";
	//	std::cout << " Include additional spectra assignment method ";
	//	std::cout << "\n [X] Use current setting";
	//	display_peptide_assignment_menu(count_itr, )
	//}

	//int perform_peptide_method_modify(string par_spectra_method_setting, vector<string>& par_v_spectra_method_setting) {
	//	string menu_selection{};
	//	while ((menu_selection != ("0"))
	//		&& (menu_selection != ("1"))
	//		&& (menu_selection != ("2"))
	//		&& (menu_selection != ("X"))) {
	//		menu_selection.clear();
	//		std::cout << "\n\n --> ";
	//		std::cin >> menu_selection;
	//	}
	//	if (menu_selection == "0") {
	//		return "PEAKS database match";
	//	}
	//	if (menu_selection == "1") {
	//		return "PEAKS de novo";
	//	}
	//	if (menu_selection == "2") {
	//		return "NOVOR de novo";
	//	}
	//	if (menu_selection == "X") {
	//		return par_spectra_method_setting;
	//	}
	//	return par_spectra_method_setting;
	//}

	//void display_peptide_assignment_menu(string par_current_setting, vector<string>& par_spectra_method_setting) {
	//	std::cout << "\n\n Current setting - ";
	//	std::cout << "[" << par_current_index << "]";
	//	std::cout << " " << par_current_setting;
	//	std::cout << "\n";
	//	std::cout << "\n [0] Modify to PEAKS database match";
	//	std::cout << "\n [1] Modify to PEAKS de novo";
	//	std::cout << "\n [2] Modify to NOVOR de novo";
	//	std::cout << "\n [X] Return";
	//}

	string display_peptide_assignment_menu(string par_peptide_method_setting_current) {
		std::cout << "\n\n Current setting - " << par_peptide_method_setting_current;
		std::cout << "\n";
		std::cout << "\n [0] PEAKS database match";
		std::cout << "\n [1] PEAKS de novo";
		std::cout << "\n [2] NOVOR de novo";
		std::cout << "\n [X] Use current setting";
		return perform_peptide_method_selection(par_peptide_method_setting_current);
	}

	string perform_peptide_method_selection(string par_spectra_method_setting) {
		string menu_selection{};
		while ((menu_selection != ("0"))
			&& (menu_selection != ("1"))
			&& (menu_selection != ("2"))
			&& (menu_selection != ("X"))) {
			menu_selection.clear();
			std::cout << "\n\n --> ";
			std::cin >> menu_selection;
		}
		if (menu_selection == "0") {
			return "PEAKS database match";
		}
		if (menu_selection == "1") {
			return "PEAKS de novo";
		}
		if (menu_selection == "2") {
			return "NOVOR de novo";
		}
		if (menu_selection == "X") {
			return par_spectra_method_setting;
		}
		return par_spectra_method_setting;
	}

	void select_settings(string& par_select_fasta, vector<string> par_v_select_peptide_assignment) {
		string menu_selection{ display_menu() };
		bool menu_continue{};
		while (!menu_continue) {
			if (menu_selection == "F") {
				par_select_fasta = display_FASTA_menu(par_select_fasta);
				display_settings(par_select_fasta, par_v_select_peptide_assignment);
				menu_selection = display_menu();
			}
			if (menu_selection == "P") {
				*par_v_select_peptide_assignment.begin() = display_peptide_assignment_menu(*par_v_select_peptide_assignment.begin());
				display_settings(par_select_fasta, par_v_select_peptide_assignment);

			}
			if (menu_selection == "X") {
				par_select_fasta = IgFamily::DEFAULT_INPUT_FASTA_DIRECTORY + par_select_fasta;
				menu_continue = true;
			}
			if (menu_selection == "!!") {
				if (IgFamily::FILESYSTEM_MODE) {
					//vector<fpf_utility::sample_transcript_and_translation> main_v_sample_transcript_and_translation = fpf_utility::parse_transcript_data();
					//fpf_utility::translate_v_transcript(main_v_sample_transcript_and_translation);
					//fpf_utility::fout_transcript_and_translation(main_v_sample_transcript_and_translation);
					fpf_genome_data::population_genome main_population_genome{};
					fpf_genome_data::create_v_genome_directory(main_population_genome);
					for (auto& itr_v_sample_genome : main_population_genome.v_sample_genome) {
						vector<fpf_genome_data::genome_data> main_v_genomic_data = fpf_genome_data::create_v_genome_data(itr_v_sample_genome.first);
						vector<fpf_genome_data::genome_analysis> main_v_genomic_analysis = fpf_genome_data::create_v_genome_analysis(main_v_genomic_data);
						fpf_genome_data::sample_genome main_sample_genome{
							&main_v_genomic_data,
							&main_v_genomic_analysis
						};
						itr_v_sample_genome.second = &main_sample_genome;
						fpf_genome_data::fout_v_genome_data(itr_v_sample_genome.first, *itr_v_sample_genome.second);
						fpf_genome_data::fout_v_genome_analysis(itr_v_sample_genome.first, *itr_v_sample_genome.second);
						fpf_genome_data::fout_v_genome_analysis_filtered(itr_v_sample_genome.first, *itr_v_sample_genome.second);
					}
				}
				menu_selection = display_menu();
			}
		}
	}
}

#endif