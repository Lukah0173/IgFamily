// * * fpf_filesystem.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_FILESYSTEM
#define	FPF_FILESYSTEM

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include "IgFamily.h"
#include "fpf_convert.h"
#include "fpf_data.h"
#include "fpf_parse.h"


namespace fpf_filesystem {

	using std::pair;
	using std::string;
	using std::vector;

	typedef fpf_convert::fileconversion_parameters fileconversion_parameters;
	typedef fpf_data::homology_data homology_data;
	typedef fpf_data::protein_data protein_data;
	typedef fpf_data::peptide_analysis peptide_analysis;
	typedef fpf_data::peptide_data peptide_data;
	typedef fpf_data::protein_analysis protein_analysis;
	typedef fpf_data::multinomial multinomial;
	typedef fpf_parse::csv_data csv_data;
	typedef fpf_parse::FASTA_data FASTA_data;

	struct filesystem;
	struct sample_analysis;
	
	struct sample_analysis {
		bool PEAKS_database_exists;
		bool PEAKS_denovo_exists;
		bool NOVOR_denovo_exists;
		string peptide_assignment_method;
		vector<protein_data> v_protein_data;
		vector<peptide_data> v_peptide_data;
		vector<peptide_analysis> v_peptide_analysis;
		vector<homology_data> v_homology_data;
		vector<protein_analysis> v_protein_analysis;
		vector<protein_analysis> v_protein_analysis_selected_by_polymorphism;
		multinomial multinomial_data;
	};

	struct filesystem {
	public:

		/* filesystem */

		pair<string, string> filesystem_id;
		string directory;
		string filename;
		string fileversion;
		vector<pair<string, string>> v_filesystem_replicates;
		size_t filesystem_replicate_count;

		/* sample factors */

		string patientstatus;
		string enzyme;	

		/* mgf conversion parameters */

		bool fileconversion = true;
		fileconversion_parameters fileconversion_parameters;

		/* de novo parameters */

		string denono_deltamass;

		/* sample analysis */

		std::vector<sample_analysis> v_sample_analysis;
	};

	void display_settings(string par_FASTA_setting, vector<string> par_v_spectra_assignment_method) {
		std::cout << "\n\n Current settings:";
		std::cout << "\n\n\n";
		std::cout << " FASTA file -                     ";
		std::cout << par_FASTA_setting;
		std::cout << "\n Spectra assignment method -      ";
		size_t count_itr{};
		for (const auto& itr_v_spectra_assignment_method : par_v_spectra_assignment_method) {
			++count_itr;
			std::cout << itr_v_spectra_assignment_method;
			if (count_itr != par_v_spectra_assignment_method.size()) {
				std::cout << " + ";
			}
		}
	}

	string display_menu() {
		std::cout << "\n\n";
		std::cout << " Select FASTA file:                         [F] ";
		std::cout << "\n Select spectra assignment method:          [P]";
		std::cout << "\n Continue with current settings:            [X]";
		string menu_selection{};
		while ((menu_selection != ("F")) && (menu_selection != ("P")) && (menu_selection != ("X"))) {
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
		std::cout << "\n [3] IGHV_IGLV_IGKV_MIGHV_MIGKV_MIGLV_20160827.fasta";
		std::cout << "\n [4] IGHV_IGKV_IGLV_MIGHV_MIGKV_MIGLV_CONT_20160827.fasta";
		std::cout << "\n [5] IGH_IGL_IGK_20160827.fasta";
		std::cout << "\n [6] IGH_IGL_IGK_CONT_20160827.fasta";
		std::cout << "\n [7] IGH_IGL_IGK_CONT_UNIPROT_20160827.fasta";
		std::cout << "\n [8] IGH_IGL_IGK_mAB_20160827.fasta";
		std::cout << "\n [9] IGH_IGL_IGK_mAB_CONT_20160827.fasta";
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
			return "IGHV_IGLV_IGKV_MIGHV_MIGKV_MIGLV_20160827.fasta";
		}
		if (menu_selection == "4") {
			return "IGHV_IGKV_IGLV_MIGHV_MIGKV_MIGLV_CONT_20160827.fasta";
		}
		if (menu_selection == "5") {
			return "IGH_IGL_IGK_20160827.fasta";
		}
		if (menu_selection == "6") {
			return "IGH_IGL_IGK_CONT_20160827.fasta";
		}
		if (menu_selection == "7") {
			return "IGH_IGL_IGK_CONT_UNIPROT_20160827.fasta";
		}
		if (menu_selection == "8") {
			return "IGH_IGL_IGK_mAB_20160827.fasta";
		}
		if (menu_selection == "9") {
			return "IGH_IGL_IGK_mAB_CONT_20160827.fasta";
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

	vector<string> read_root_dir(string par_root_directory) {
		vector<string> temp_v_IgFamily_root{};
		string fin_IgFamily_root{};
		if (IgFamily::FILESYSTEM_MODE) {
			std::ifstream fin_input_csv(par_root_directory);
			char stream_IgFamily_root{};
			while (fin_input_csv.std::istream::get(stream_IgFamily_root)) {
				if ((stream_IgFamily_root != '\n') && (stream_IgFamily_root != ',')) {
					fin_IgFamily_root += stream_IgFamily_root;
				}
				if (stream_IgFamily_root == ',') {
					temp_v_IgFamily_root.push_back(fin_IgFamily_root);
					std::cout << "\n * " << fin_IgFamily_root;
					fin_IgFamily_root.clear();
				}
			}
		}
		return temp_v_IgFamily_root;
	}

	vector<filesystem> read_filesystem(vector<string> par_root_directory) {
		filesystem temp_filesystem{};
		vector<filesystem> temp_v_filesystem{};
		if (IgFamily::FILESYSTEM_MODE) {
			std::cout << "\n";
			for (const auto& itr_root_directory : par_root_directory) {
				string stream_filesystem = itr_root_directory + "filesystem.data";
				std::ifstream fin_input_filesystem(stream_filesystem);
				char read_fin_filesystem{};
				string stream_fin_filesystem{};
				size_t switch_fin_filesystem{};
				pair<string, string> temp_filesystem_id{};
				string temp_filesystem_date{};
				string temp_filesystem_filename{};
				string temp_filesystem_version{};
				string temp_filesystem_status{};
				string temp_filesystem_enzyme{};
				string temp_filesystem_denovo_deltamass{};
				string temp_filesystem_replicatedate{};
				size_t switch__filesystem_replicatepair{};
				vector<pair<string, string>> v_filesystem_replicates{};
				while (fin_input_filesystem.std::istream::get(read_fin_filesystem)) {
					if ((read_fin_filesystem != ';') && (read_fin_filesystem != '\n') && (read_fin_filesystem != ',')) {
						stream_fin_filesystem += read_fin_filesystem;
					}
					if (switch__filesystem_replicatepair == 2) {
						switch__filesystem_replicatepair = 0;
					}
					if ((switch_fin_filesystem == 1) && (read_fin_filesystem == ',')) {
						temp_filesystem_date = stream_fin_filesystem;
						stream_fin_filesystem.clear();
					}
					if (stream_fin_filesystem == "ID: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 1;
					}
					if (stream_fin_filesystem == "FILE: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 2;
					}
					if (stream_fin_filesystem == "VERSION: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 3;
					}
					if ((switch_fin_filesystem == 4) && (read_fin_filesystem == ',') && (switch__filesystem_replicatepair == 1)) {
						v_filesystem_replicates.push_back(std::make_pair(temp_filesystem_replicatedate, stream_fin_filesystem));
						switch__filesystem_replicatepair = 2;
						temp_filesystem_replicatedate.clear();
						stream_fin_filesystem.clear();
					}
					if ((switch_fin_filesystem == 4) && (read_fin_filesystem == ';') && (switch__filesystem_replicatepair == 0)) {
						temp_filesystem_replicatedate = stream_fin_filesystem;
						switch__filesystem_replicatepair = 1;
						stream_fin_filesystem.clear();
					}
					if (stream_fin_filesystem == "REPLICATES: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 4;
					}
					if (stream_fin_filesystem == "STATUS: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 5;
					}
					if (stream_fin_filesystem == "ENZYME: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 6;
					}
					if (stream_fin_filesystem == "DENOVO_DELTAMASS: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 7;
					}
					if (read_fin_filesystem == ';') {
						if (switch_fin_filesystem == 1) {
							temp_filesystem_id = make_pair(temp_filesystem_date, stream_fin_filesystem);
							temp_filesystem_date.clear();
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 2) {
							temp_filesystem_filename = stream_fin_filesystem;
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 3) {
							temp_filesystem_version = stream_fin_filesystem;
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 4) {
							if (stream_fin_filesystem != "") {
								v_filesystem_replicates.push_back(std::make_pair(temp_filesystem_replicatedate, stream_fin_filesystem));
							}
							switch__filesystem_replicatepair = size_t();
							temp_filesystem_replicatedate.clear();
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 5) {
							temp_filesystem_status = stream_fin_filesystem;
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 6) {
							temp_filesystem_enzyme = stream_fin_filesystem;
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 7) {
							temp_filesystem_denovo_deltamass = stream_fin_filesystem;
							stream_fin_filesystem.clear();
							temp_filesystem.directory = itr_root_directory;
							temp_filesystem.filesystem_id = temp_filesystem_id;
							temp_filesystem.filename = temp_filesystem_filename;
							temp_filesystem.fileversion = temp_filesystem_version;
							temp_filesystem.patientstatus = temp_filesystem_status;
							temp_filesystem.v_filesystem_replicates = v_filesystem_replicates;
							temp_filesystem.enzyme = temp_filesystem_enzyme;
							temp_filesystem.denono_deltamass = temp_filesystem_denovo_deltamass;
							temp_filesystem.filesystem_replicate_count = size_t{ 1 };
							v_filesystem_replicates.clear();
							temp_v_filesystem.push_back(temp_filesystem);
							std::cout << "\n * reading: " << temp_filesystem_filename << "   version - " << temp_filesystem_version;
							if (temp_filesystem_version == IgFamily::version) {
								std::cout << "   ! up to date !";
							}
						}
					}
				}
				fin_input_filesystem.clear();
				fin_input_filesystem.seekg(0, std::ios::beg);
			}
		}
		else {
			std::cout << "\n * reading: Local file";
			temp_filesystem.filename = "Local_file";
			temp_filesystem.fileversion = IgFamily::version;
			temp_v_filesystem.push_back(temp_filesystem);
		}
		return temp_v_filesystem;
	}

	void perform_fileconversion(filesystem& par_filesystem) {
		par_filesystem.fileconversion_parameters = fpf_convert::create_fileconversion_parameters(fpf_convert::prompt_defaultconversion());
		string fileconversion_command{};
		fileconversion_command += "msconvert.exe ";
		fileconversion_command += "\"Z:\\Lukah_Dykes\\IgFamily\\";
		fileconversion_command += par_filesystem.directory;
		fileconversion_command += par_filesystem.filename;
		fileconversion_command += ".wiff\"";
		////fileconversion_command += " --64";
		////fileconversion_command += " --mz64";
		//fileconversion_command += " -v";
		//fileconversion_command += " --mgf";
		//if (par_filesystem.fileconversion_parameters.peakpicking.peakpicking) {
		//	//fileconversion_command += " --filter \"peakPicking cwt ";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion.peakpicking.peakpicking_mslevel_from);
		//	//fileconversion_command += "-";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion.peakpicking.peakpicking_mslevel_to);
		//	//fileconversion_command += "\"";
		//}
		//if (par_filesystem.fileconversion_parameters.threshold.intensity_threshold) {
		//	//fileconversion_command += " --filter \"threshold absolute ";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.threshold.threshold);
		//	//fileconversion_command += " most-intense\"";
		//}
		//if (par_filesystem.fileconversion_parameters.ms2denoise.ms2denoise) {
		//	//fileconversion_command += " --filter \"MS2Denoise ";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.ms2denoise.ms2denoise_peaksinwindow);
		//	//fileconversion_command += " ";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.ms2denoise.ms2denoise_windowwidth);
		//	//fileconversion_command += " true\"";
		//}
		//if (par_filesystem.fileconversion_parameters.ms2deisotope.ms2deisotope) {
		//	//fileconversion_command += " --filter MS2Deisotope";
		//}
		//if (par_filesystem.fileconversion_parameters.chargestatepredictor.chargestatepredictor) {
		//	fileconversion_command += " --filter \"chargeStatePredictor true ";
		//	fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.chargestatepredictor.chargestatepredictor_maxcharge);
		//	fileconversion_command += " ";
		//	fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.chargestatepredictor.chargestatepredictor_mincharge);
		//	fileconversion_command += " ";
		//	fileconversion_command += "0.9";
		//	fileconversion_command += "\"";
		//}
		//fileconversion_command += " outdir=Z:\\Lukah_Dykes\\IgFamily\\";
		//fileconversion_command += par_filesystem.directory;
		fileconversion_command += " --mgf --filter \"chargeStatePredictor true 3 2 0.9\"";
		std::cout << "\n\n" << fileconversion_command;
		fpf_convert::sys_msconvert(fileconversion_command, par_filesystem.directory);
	}

	string read_filesystem_PEAKS_database_peptides(string par_root_directory) {
		string temp_root_PEAKS_database_peptides = par_root_directory + "database_peptides_PEAKS.csv";
		return temp_root_PEAKS_database_peptides;
	}

	string read_filesystem_PEAKS_denovo_peptides(string par_root_directory) {
		string temp_root_PEAKS_denovo_peptides{};
		temp_root_PEAKS_denovo_peptides = par_root_directory + "denovo_peptides_PEAKS.csv";
		return temp_root_PEAKS_denovo_peptides;
	}

	string read_filesystem_NOVOR_denovo_peptides(string par_root_directory) {
		string temp_root_NOVOR_database_peptides{};
		temp_root_NOVOR_database_peptides = par_root_directory + "denovo_peptides_NOVOR.csv";
		return temp_root_NOVOR_database_peptides;
	}

	vector<csv_data> parse_filesystem_PEAKS_database_peptides(string par_fin_root_directory) {
		vector<csv_data> temp_v_csv_proteinpeptides{};
		std::ifstream fin_input_csv(par_fin_root_directory);
		temp_v_csv_proteinpeptides = fpf_parse::parse_proteinpeptides(fin_input_csv, par_fin_root_directory);
		return temp_v_csv_proteinpeptides;
	}

	vector<csv_data> parse_filesystem_PEAKS_denovo_peptides(string par_fin_root_directory) {
		vector<csv_data> temp_csv_PEAKS_denovo_peptides{};
		std::ifstream fin_input_csv(par_fin_root_directory);
		temp_csv_PEAKS_denovo_peptides = fpf_parse::parse_csv_PEAKS_denovopeptides(fin_input_csv, par_fin_root_directory);
		return temp_csv_PEAKS_denovo_peptides;
	}

	vector<csv_data> parse_filesystem_NOVOR_denovo_peptides(string par_fin_root_directory) {
		vector<csv_data> temp_csv_NOVOR_denovo_peptides{};
		std::ifstream fin_input_csv(par_fin_root_directory);
		temp_csv_NOVOR_denovo_peptides = fpf_parse::parse_csv_NOVOR_denovopeptides(fin_input_csv, par_fin_root_directory);
		return temp_csv_NOVOR_denovo_peptides;
	}

	void fout_filesystem(filesystem par_filesystem) {
		if (IgFamily::FILESYSTEM_MODE) {
			string output_filesystem = par_filesystem.directory + "filesystem.data";
			std::ofstream fout_filesystem;
			fout_filesystem.open(output_filesystem);
			fout_filesystem << "ID: " << std::get<0>(par_filesystem.filesystem_id) << "," << std::get<1>(par_filesystem.filesystem_id) << ";\n";
			fout_filesystem << "FILE: " << par_filesystem.filename << ";\n";
			fout_filesystem << "VERSION: " << IgFamily::version << ";\n";
			fout_filesystem << "REPLICATES: ";
			for (const auto& itr_v_p_replicates : par_filesystem.v_filesystem_replicates) {
				fout_filesystem << std::get<0>(itr_v_p_replicates) << ",";
				fout_filesystem << std::get<1>(itr_v_p_replicates) << ";";
			}
			fout_filesystem << "\n";
			fout_filesystem << "STATUS: " << par_filesystem.patientstatus << ";\n";
			fout_filesystem << "ENZYME: " << par_filesystem.enzyme << ";\n";
			fout_filesystem << "DENOVO_DELTAMASS: " << par_filesystem.denono_deltamass << ";\n";
		}
	}
}

#endif