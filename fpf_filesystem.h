// * * fpf_filesystem.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_FILESYSTEM
#define	FPF_FILESYSTEM

#include <cstdlib> // provides - size_t
#include <vector> // provides - vector
#include <utility> // provides - std::pair
#include <iostream> // provides - std::istream
#include <algorithm> // provides - std::find
#include <tuple> // provides - std::tuple

#include "IgFamily.h"
#include "fpf_parse.h"
#include "fpf_data.h"
#include "fpf_convert.h"


namespace fpf_filesystem {

	using std::string;
	using std::vector;

	typedef fpf_convert::fileconversion_parameters fileconversion_parameters;
	typedef fpf_data::multinomial_category multinomial_category;
	typedef fpf_data::peptide_data peptide_data;
	typedef fpf_data::blastp_data blastp_data;
	typedef fpf_data::multinomial multinomial;
	typedef fpf_data::category_report category_report;

	struct filesystem;

	struct filesystem {
	public:
		bool proteinpeptides_exist;
		bool denovopeptides_exist;

		/* filesystem */

		string directory;
		string filename;
		string fileversion;

		/* sample factors */

		string patientstatus;
		string enzyme;	

		/* mgf conversion parameters */

		bool fileconversion = true;
		fileconversion_parameters fileconversion_parameters;

		/* de novo parameters */

		string denono_deltamass;

		vector<multinomial_category> v_multinomial_category;
		vector<peptide_data> v_peptide_data;
		std::pair<string, string> filesystem_id;
		vector<std::pair<string, string>> v_filesystem_replicates;
		size_t filesystem_replicate_count;
		vector<blastp_data> v_blastp_data;
		multinomial multinomial_data;
		vector<category_report> v_category_report;
	};

	string display_menu() {
		std::cout << " FASTA utilities:     [F] ";
		std::cout << "\n Continue:            [X]";
		string menu_selection;
		while ((menu_selection != ("F")) && (menu_selection != ("X"))) {
			menu_selection.clear();
			std::cout << "\n\n Input selection: \n\n > ";
			std::cin >> menu_selection;
		}
		return menu_selection;
	}

	int perform_menu_selection(string par_menu_selection) {
		if (par_menu_selection == "F") {
			std::cout << "\n Read FASTA format:              [R] ";
			std::cout << "\n Output custom FASTA format:     [C] ";
			string menu_selection = string();
			while ((menu_selection != ("R")) && (menu_selection != ("C"))) {
				menu_selection.clear();
				std::cout << "\n\n Input selection: \n\n > ";
				std::cin >> menu_selection;
			}
			if (menu_selection == "R") {
				fpf_parse::check_FASTA_format(IgFamily::INPUT_FASTA);
			}
			if (menu_selection == "C") {
				fpf_parse::output_custom_FASTA_format(IgFamily::INPUT_FASTA);
			}
			return 0;
		}
		if (par_menu_selection == "X") {
			return 0;
		}
		return 1;
	}

	vector<string> read_root_dir(string par_IgFamily_root_dir) {
		std::ifstream fin_input_csv(par_IgFamily_root_dir);
		vector<string> temp_v_IgFamily_root;
		string fin_IgFamily_root = string();
		char stream_IgFamily_root = char();
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
		return temp_v_IgFamily_root;
	}

	vector<filesystem> read_filesystem(vector<string> par_root_dir) {
		filesystem con_s_filesystem;
		vector<filesystem> con_v_s_filesystem;
		std::cout << "\n";
		for (vector<string>::iterator itr_root_dir = par_root_dir.begin(); itr_root_dir != par_root_dir.end(); ++itr_root_dir) {
			string stream_str_filesytem = *itr_root_dir + "filesystem.data";
			std::ifstream fin_input_filesystem(stream_str_filesytem);
			char c_stream_fin_filesystem = char();
			string str_stream_fin_filesystem = string();
			size_t sw_stream_fin_filesystem = size_t();
			std::pair<string, string> con_p_filesystem_id = std::pair<string, string>();
			string con_str_filesystem_date = string();
			string con_str_filesystem_filename = string();
			string con_str_filesystem_version = string();
			string con_str_filesystem_status = string();
			string con_str_filesystem_enzyme = string();
			string con_str_filesystem_denovo_deltamass = string();
			string con_str_filesystem_replicatedate = string();
			size_t sw_filesystem_replicatepair = size_t();
			vector<std::pair<string, string>> v_p_filesystem_replicates = vector<std::pair<string, string>>();
			while (fin_input_filesystem.std::istream::get(c_stream_fin_filesystem)) {
				if ((c_stream_fin_filesystem != ';') && (c_stream_fin_filesystem != '\n') && (c_stream_fin_filesystem != ',')) {
					str_stream_fin_filesystem += c_stream_fin_filesystem;
				}
				if (sw_filesystem_replicatepair == 2) {
					sw_filesystem_replicatepair = 0;
				}
				if ((sw_stream_fin_filesystem == 1) && (c_stream_fin_filesystem == ',')) {
					con_str_filesystem_date = str_stream_fin_filesystem;
					str_stream_fin_filesystem.clear();
				}
				if (str_stream_fin_filesystem == "ID: ") {
					str_stream_fin_filesystem.clear();
					sw_stream_fin_filesystem = 1;
				}
				if (str_stream_fin_filesystem == "FILE: ") {
					str_stream_fin_filesystem.clear();
					sw_stream_fin_filesystem = 2;
				}
				if (str_stream_fin_filesystem == "VERSION: ") {
					str_stream_fin_filesystem.clear();
					sw_stream_fin_filesystem = 3;
				}
				if ((sw_stream_fin_filesystem == 4) && (c_stream_fin_filesystem == ',') && (sw_filesystem_replicatepair == 1)) {
					v_p_filesystem_replicates.push_back(std::make_pair(con_str_filesystem_replicatedate, str_stream_fin_filesystem));
					sw_filesystem_replicatepair = 2;
					con_str_filesystem_replicatedate.clear();
					str_stream_fin_filesystem.clear();
				}
				if ((sw_stream_fin_filesystem == 4) && (c_stream_fin_filesystem == ',') && (sw_filesystem_replicatepair == 0)) {
					con_str_filesystem_replicatedate = str_stream_fin_filesystem;
					sw_filesystem_replicatepair = 1;
					str_stream_fin_filesystem.clear();
				}
				if (str_stream_fin_filesystem == "REPLICATES: ") {
					str_stream_fin_filesystem.clear();
					sw_stream_fin_filesystem = 4;
				}
				if (str_stream_fin_filesystem == "STATUS: ") {
					str_stream_fin_filesystem.clear();
					sw_stream_fin_filesystem = 5;
				}
				if (str_stream_fin_filesystem == "ENZYME: ") {
					str_stream_fin_filesystem.clear();
					sw_stream_fin_filesystem = 6;
				}
				if (str_stream_fin_filesystem == "DENOVO_DELTAMASS: ") {
					str_stream_fin_filesystem.clear();
					sw_stream_fin_filesystem = 7;
				}
				if (c_stream_fin_filesystem == ';') {
					if (sw_stream_fin_filesystem == 1) {
						con_p_filesystem_id = make_pair(con_str_filesystem_date, str_stream_fin_filesystem);
						con_str_filesystem_date.clear();
						str_stream_fin_filesystem.clear();
					}
					if (sw_stream_fin_filesystem == 2) {
						con_str_filesystem_filename = str_stream_fin_filesystem;
						str_stream_fin_filesystem.clear();
					}
					if (sw_stream_fin_filesystem == 3) {
						con_str_filesystem_version = str_stream_fin_filesystem;
						str_stream_fin_filesystem.clear();
					}
					if (sw_stream_fin_filesystem == 4) {
						if (str_stream_fin_filesystem != "") {
							v_p_filesystem_replicates.push_back(std::make_pair(con_str_filesystem_replicatedate, str_stream_fin_filesystem));
						}
						sw_filesystem_replicatepair = size_t();
						con_str_filesystem_replicatedate.clear();
						str_stream_fin_filesystem.clear();
					}
					if (sw_stream_fin_filesystem == 5) {
						con_str_filesystem_status = str_stream_fin_filesystem;
						str_stream_fin_filesystem.clear();
					}
					if (sw_stream_fin_filesystem == 6) {
						con_str_filesystem_enzyme = str_stream_fin_filesystem;
						str_stream_fin_filesystem.clear();
					}
					if (sw_stream_fin_filesystem == 7) {
						con_str_filesystem_denovo_deltamass = str_stream_fin_filesystem;
						str_stream_fin_filesystem.clear();
						con_s_filesystem.directory = *itr_root_dir;
						con_s_filesystem.filesystem_id = con_p_filesystem_id;
						con_s_filesystem.filename = con_str_filesystem_filename;
						con_s_filesystem.fileversion = con_str_filesystem_version;
						con_s_filesystem.patientstatus = con_str_filesystem_status;
						con_s_filesystem.v_filesystem_replicates = v_p_filesystem_replicates;
						con_s_filesystem.enzyme = con_str_filesystem_enzyme;
						con_s_filesystem.denono_deltamass = con_str_filesystem_denovo_deltamass;
						con_s_filesystem.filesystem_replicate_count = size_t{ 1 };
						v_p_filesystem_replicates.clear();
						con_v_s_filesystem.push_back(con_s_filesystem);
						std::cout << "\n * reading: " << con_str_filesystem_filename << "   version - " << con_str_filesystem_version;
						if (con_str_filesystem_version == IgFamily::version) {
							std::cout << "   ! up to date !";
						}
					}
				}
			}
			fin_input_filesystem.clear();
			fin_input_filesystem.seekg(0, std::ios::beg);
		}
		return con_v_s_filesystem;
	}

	void perform_fileconversion(filesystem& par_filesystem) {
		par_filesystem.fileconversion_parameters = fpf_convert::create_fileconversion_parameters(fpf_convert::prompt_defaultconversion());
		string str_fileconversion_command = string();
		str_fileconversion_command += "msconvert.exe ";
		str_fileconversion_command += "\"Z:\\Lukah_Dykes\\IgFamily\\";
		str_fileconversion_command += par_filesystem.directory;
		str_fileconversion_command += par_filesystem.filename;
		str_fileconversion_command += ".wiff\"";
		str_fileconversion_command += " --64";
		str_fileconversion_command += " --mz64";
		str_fileconversion_command += " -v";
		str_fileconversion_command += " --mgf";
		if (par_filesystem.fileconversion_parameters.s_peakpicking.b_peakpicking) {
			str_fileconversion_command += " --filter \"peakPicking cwt ";
			//str_fileconversion_command += std::to_string(par_filesystem.fileconversion.s_peakpicking.st_peakpicking_mslevel_from);
			//str_fileconversion_command += "-";
			//str_fileconversion_command += std::to_string(par_filesystem.fileconversion.s_peakpicking.st_peakpicking_mslevel_to);
			str_fileconversion_command += "\"";
		}
		if (par_filesystem.fileconversion_parameters.s_threshold.b_threshold) {
			str_fileconversion_command += " --filter \"threshold absolute ";
			str_fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.s_threshold.st_threshold);
			str_fileconversion_command += " most-intense\"";
		}
		if (par_filesystem.fileconversion_parameters.s_ms2denoise.b_ms2denoise) {
			str_fileconversion_command += " --filter \"MS2Denoise ";
			str_fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.s_ms2denoise.st_ms2denoise_peaksinwindow);
			str_fileconversion_command += " ";
			str_fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.s_ms2denoise.st_ms2denoise_windowwidth);
			str_fileconversion_command += " true\"";
		}
		if (par_filesystem.fileconversion_parameters.s_ms2deisotope.b_ms2deisotope) {
			str_fileconversion_command += " --filter MS2Deisotope";
		}
		if (par_filesystem.fileconversion_parameters.s_chargestatepredictor.b_chargestatepredictor) {
			//str_fileconversion_command += " --filter \"chargeStatePredictor true ";
			//str_fileconversion_command += std::to_string(par_filesystem.fileconversion.s_chargestatepredictor.st_chargestatepredictor_mincharge);
			//str_fileconversion_command += " ";
			//str_fileconversion_command += std::to_string(par_filesystem.fileconversion.s_chargestatepredictor.st_chargestatepredictor_maxcharge);
			//str_fileconversion_command += " ";
			//str_fileconversion_command += std::to_string(par_filesystem.fileconversion.s_chargestatepredictor.d_chargestatepredictor_chargefraction);
			//str_fileconversion_command += "\"";
		}
		str_fileconversion_command += " -o Z:\\Lukah_Dykes\\IgFamily\\";
		str_fileconversion_command += par_filesystem.directory;
		//std::cout << "\n\n" << str_fileconversion_command;
		fpf_convert::sys_msconvert(str_fileconversion_command, par_filesystem.directory);
	}

	string read_filesystem_proteinpeptides(string par_IgFamily_root_dir) {
		string con_str_root_proteinpeptides = par_IgFamily_root_dir + "protein_peptides.csv";
		return con_str_root_proteinpeptides;
	}

	string read_filesystem_denovopeptides(string par_IgFamily_root_dir) {
		string con_str_root_denovopeptides;
		if (IgFamily::NOVOR_DENOVO) {
			con_str_root_denovopeptides = par_IgFamily_root_dir + "denovo_peptides_NOVOR.csv";
		}
		else {
			con_str_root_denovopeptides = par_IgFamily_root_dir + "denovo_peptides.csv";
		}
		return con_str_root_denovopeptides;
	}

	vector<fpf_parse::parse_peptides_csv_type> parse_filesystem_proteinpeptides(string par_str_fin_root) {
		vector<fpf_parse::parse_peptides_csv_type> con_v_c_parse_csv;
		std::ifstream fin_input_csv(par_str_fin_root);
		con_v_c_parse_csv = fpf_parse::parse_proteinpeptides(fin_input_csv, par_str_fin_root);
		return con_v_c_parse_csv;
	}

	vector<fpf_parse::parse_peptides_csv_type> parse_filesystem_denovopeptides(string par_str_fin_root) {
		vector<fpf_parse::parse_peptides_csv_type> con_v_c_parse_csv;
		std::ifstream fin_input_csv(par_str_fin_root);
		if (IgFamily::NOVOR_DENOVO) {
			con_v_c_parse_csv = fpf_parse::parse_NOVOR_denovopeptides_csv(fin_input_csv, par_str_fin_root);
		}
		else {
			con_v_c_parse_csv = fpf_parse::parse_PEAKS_denovopeptides_csv(fin_input_csv, par_str_fin_root);
		}
		return con_v_c_parse_csv;
	}

	void fout_filesystem(filesystem par_filesystem) {
		string output_v_c_analysis = par_filesystem.directory + "filesystem.data";
		std::ofstream fout_v_s_filesystem;
		fout_v_s_filesystem.open(output_v_c_analysis);
		fout_v_s_filesystem << "ID: " << std::get<0>(par_filesystem.filesystem_id) << "," << std::get<1>(par_filesystem.filesystem_id) << ";\n";
		fout_v_s_filesystem << "FILE: " << par_filesystem.filename << ";\n";
		fout_v_s_filesystem << "VERSION: " << IgFamily::version << ";\n";
		fout_v_s_filesystem << "REPLICATES: ";
		for (vector<std::pair<string, string>>::iterator itr_v_p_replicates = par_filesystem.v_filesystem_replicates.begin(); itr_v_p_replicates != par_filesystem.v_filesystem_replicates.end(); ++itr_v_p_replicates) {
			fout_v_s_filesystem << std::get<0>(*itr_v_p_replicates) << "," << std::get<1>(*itr_v_p_replicates);
			if ((itr_v_p_replicates + 1) != par_filesystem.v_filesystem_replicates.end()) {
				fout_v_s_filesystem << ",";
			}
		}
		fout_v_s_filesystem << ";\n";
		fout_v_s_filesystem << "STATUS: " << par_filesystem.patientstatus << ";\n";
		fout_v_s_filesystem << "ENZYME: " << par_filesystem.enzyme << ";\n";
		fout_v_s_filesystem << "DENOVO_DELTAMASS: " << par_filesystem.denono_deltamass << ";\n";
	}
}

#endif