// * * fpf_filesystem.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_FILESYSTEM
#define	FPF_FILESYSTEM
#include <cstdlib> // provides - size_t
#include <vector> // provides - std::vector
#include <utility> // provides - std::pair
#include <iostream> // provides - std::istream
#include <algorithm> // provides - std::find
#include <tuple> // provides - std::tuple
#include "IgFamily.h"
#include "fpf_parse.h"
#include "fpf_data.h"
#include "fpf_convert.h"



namespace fpf_filesystem {

	struct filesystem_type;

	typedef std::string string_type;
	typedef size_t size_type;
	typedef fpf_convert::fileconversion_type fileconversion_type;
	typedef fpf_data::multinomial_category_data_type multinomial_category_data_type;
	typedef fpf_data::peptide_data_type peptide_data_type;
	typedef fpf_data::blastp_type blastp_type;
	typedef fpf_data::multinomial_type multinomial_type;
	typedef fpf_data::blastp_type blastp_type;
	typedef fpf_data::report_type report_type;

	struct filesystem_type {
	public:
		bool b_proteinpeptides_exist;
		bool b_denovopeptides_exist;

		/* filesystem */

		string_type str_directory;
		string_type str_filename;
		string_type str_fileversion;

		/* sample factors */

		string_type str_patientstatus;
		string_type str_enzyme;	

		/* mgf conversion parameters */

		bool b_fileconversion = true;
		fileconversion_type s_fileconversion;

		/* de novo parameters */

		string_type str_denono_deltamass;

		std::vector<multinomial_category_data_type> v_c_multinomial_catagory;
		std::vector<multinomial_category_data_type> v_c_multinomial_catagory_distinct;
		std::vector<peptide_data_type> v_s_peptide_data;
		std::pair<string_type, string_type> p_filesystemid;
		std::vector<std::pair<string_type, string_type>> v_p_replicates;
		size_type st_replicate_count;
		std::vector<blastp_type> v_s_blastp;
		std::vector<report_type> v_s_report;
		multinomial_type s_multinomial;
	};

	string_type display_menu() {
		std::cout << " FASTA utilities:     [F] ";
		std::cout << "\n Continue:            [X]";
		string_type str_menu_selection;
		while ((str_menu_selection != ("F")) && (str_menu_selection != ("X"))) {
			str_menu_selection.clear();
			std::cout << "\n\n Input selection: \n\n > ";
			std::cin >> str_menu_selection;
		}
		return str_menu_selection;
	}

	int perform_menu_selection(string_type par_str_menu_selection) {
		if (par_str_menu_selection == "F") {
			std::cout << "\n Read FASTA format:              [R] ";
			std::cout << "\n Output custom FASTA format:     [C] ";
			string_type str_menu_selection = string_type();
			while ((str_menu_selection != ("R")) && (str_menu_selection != ("C"))) {
				str_menu_selection.clear();
				std::cout << "\n\n Input selection: \n\n > ";
				std::cin >> str_menu_selection;
			}
			if (str_menu_selection == "R") {
				fpf_parse::check_FASTA_format(IgFamily::INPUT_FASTA);
			}
			if (str_menu_selection == "C") {
				fpf_parse::output_custom_FASTA_format(IgFamily::INPUT_FASTA);
			}
			return 0;
		}
		if (par_str_menu_selection == "X") {
			return 0;
		}
		return 1;
	}

	std::vector<string_type> read_root_dir(string_type par_IgFamily_root_dir) {
		std::ifstream fin_input_csv(par_IgFamily_root_dir);
		std::vector<string_type> con_v_str_fin_IgFamily_root;
		string_type str_stream_fin_IgFamily_root = string_type();
		char c_stream_fin_IgFamily_root = char();
		while (fin_input_csv.std::istream::get(c_stream_fin_IgFamily_root)) {
			if ((c_stream_fin_IgFamily_root != '\n') && (c_stream_fin_IgFamily_root != ',')) {
				str_stream_fin_IgFamily_root += c_stream_fin_IgFamily_root;
			}
			if (c_stream_fin_IgFamily_root == ',') {
				con_v_str_fin_IgFamily_root.push_back(str_stream_fin_IgFamily_root);
				std::cout << "\n * " << str_stream_fin_IgFamily_root;
				str_stream_fin_IgFamily_root.clear();
			}
		}
		return con_v_str_fin_IgFamily_root;
	}

	std::vector<filesystem_type> read_filesystem(std::vector<string_type> par_root_dir) {
		filesystem_type con_s_filesystem;
		std::vector<filesystem_type> con_v_s_filesystem;
		std::cout << "\n";
		for (std::vector<string_type>::iterator itr_root_dir = par_root_dir.begin(); itr_root_dir != par_root_dir.end(); ++itr_root_dir) {
			string_type stream_str_filesytem = *itr_root_dir + "filesystem.data";
			std::ifstream fin_input_filesystem(stream_str_filesytem);
			char c_stream_fin_filesystem = char();
			string_type str_stream_fin_filesystem = string_type();
			size_type sw_stream_fin_filesystem = size_type();
			std::pair<string_type, string_type> con_p_filesystem_id = std::pair<string_type, string_type>();
			string_type con_str_filesystem_date = string_type();
			string_type con_str_filesystem_filename = string_type();
			string_type con_str_filesystem_version = string_type();
			string_type con_str_filesystem_status = string_type();
			string_type con_str_filesystem_enzyme = string_type();
			string_type con_str_filesystem_denovo_deltamass = string_type();
			string_type con_str_filesystem_replicatedate = string_type();
			size_type sw_filesystem_replicatepair = size_type();
			std::vector<std::pair<string_type, string_type>> v_p_filesystem_replicates = std::vector<std::pair<string_type, string_type>>();
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
						sw_filesystem_replicatepair = size_type();
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
						con_s_filesystem.str_directory = *itr_root_dir;
						con_s_filesystem.p_filesystemid = con_p_filesystem_id;
						con_s_filesystem.str_filename = con_str_filesystem_filename;
						con_s_filesystem.str_fileversion = con_str_filesystem_version;
						con_s_filesystem.str_patientstatus = con_str_filesystem_status;
						con_s_filesystem.v_p_replicates = v_p_filesystem_replicates;
						con_s_filesystem.str_enzyme = con_str_filesystem_enzyme;
						con_s_filesystem.str_denono_deltamass = con_str_filesystem_denovo_deltamass;
						con_s_filesystem.st_replicate_count = size_type{ 1 };
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

	void perform_fileconversion(filesystem_type& par_s_filesystem) {
		par_s_filesystem.s_fileconversion = fpf_convert::create_s_fileconversion(fpf_convert::prompt_b_defaultconversion());
		string_type str_fileconversion_command = string_type();
		str_fileconversion_command += "msconvert.exe ";
		str_fileconversion_command += "\"Z:\\Lukah_Dykes\\IgFamily\\";
		str_fileconversion_command += par_s_filesystem.str_directory;
		str_fileconversion_command += par_s_filesystem.str_filename;
		str_fileconversion_command += ".wiff\"";
		str_fileconversion_command += " --64";
		str_fileconversion_command += " --mz64";
		str_fileconversion_command += " -v";
		str_fileconversion_command += " --mgf";
		if (par_s_filesystem.s_fileconversion.s_peakpicking.b_peakpicking) {
			str_fileconversion_command += " --filter \"peakPicking cwt ";
			//str_fileconversion_command += std::to_string(par_s_filesystem.s_fileconversion.s_peakpicking.st_peakpicking_mslevel_from);
			//str_fileconversion_command += "-";
			//str_fileconversion_command += std::to_string(par_s_filesystem.s_fileconversion.s_peakpicking.st_peakpicking_mslevel_to);
			str_fileconversion_command += "\"";
		}
		if (par_s_filesystem.s_fileconversion.s_threshold.b_threshold) {
			str_fileconversion_command += " --filter \"threshold absolute ";
			str_fileconversion_command += std::to_string(par_s_filesystem.s_fileconversion.s_threshold.st_threshold);
			str_fileconversion_command += " most-intense\"";
		}
		if (par_s_filesystem.s_fileconversion.s_ms2denoise.b_ms2denoise) {
			str_fileconversion_command += " --filter \"MS2Denoise ";
			str_fileconversion_command += std::to_string(par_s_filesystem.s_fileconversion.s_ms2denoise.st_ms2denoise_peaksinwindow);
			str_fileconversion_command += " ";
			str_fileconversion_command += std::to_string(par_s_filesystem.s_fileconversion.s_ms2denoise.st_ms2denoise_windowwidth);
			str_fileconversion_command += " true\"";
		}
		if (par_s_filesystem.s_fileconversion.s_ms2deisotope.b_ms2deisotope) {
			str_fileconversion_command += " --filter MS2Deisotope";
		}
		if (par_s_filesystem.s_fileconversion.s_chargestatepredictor.b_chargestatepredictor) {
			//str_fileconversion_command += " --filter \"chargeStatePredictor true ";
			//str_fileconversion_command += std::to_string(par_s_filesystem.s_fileconversion.s_chargestatepredictor.st_chargestatepredictor_mincharge);
			//str_fileconversion_command += " ";
			//str_fileconversion_command += std::to_string(par_s_filesystem.s_fileconversion.s_chargestatepredictor.st_chargestatepredictor_maxcharge);
			//str_fileconversion_command += " ";
			//str_fileconversion_command += std::to_string(par_s_filesystem.s_fileconversion.s_chargestatepredictor.d_chargestatepredictor_chargefraction);
			//str_fileconversion_command += "\"";
		}
		str_fileconversion_command += " -o Z:\\Lukah_Dykes\\IgFamily\\";
		str_fileconversion_command += par_s_filesystem.str_directory;
		//std::cout << "\n\n" << str_fileconversion_command;
		fpf_convert::sys_msconvert(str_fileconversion_command, par_s_filesystem.str_directory);
	}

	string_type read_filesystem_proteinpeptides(string_type par_IgFamily_root_dir) {
		string_type con_str_root_proteinpeptides = par_IgFamily_root_dir + "protein_peptides.csv";
		return con_str_root_proteinpeptides;
	}

	string_type read_filesystem_denovopeptides(string_type par_IgFamily_root_dir) {
		string_type con_str_root_denovopeptides;
		if (IgFamily::NOVOR_DENOVO) {
			con_str_root_denovopeptides = par_IgFamily_root_dir + "denovo_peptides_NOVOR.csv";
		}
		else {
			con_str_root_denovopeptides = par_IgFamily_root_dir + "denovo_peptides.csv";
		}
		return con_str_root_denovopeptides;
	}

	std::vector<fpf_parse::parse_peptides_csv_type> parse_filesystem_proteinpeptides(string_type par_str_fin_root) {
		std::vector<fpf_parse::parse_peptides_csv_type> con_v_c_parse_csv;
		std::ifstream fin_input_csv(par_str_fin_root);
		con_v_c_parse_csv = fpf_parse::parse_proteinpeptides(fin_input_csv, par_str_fin_root);
		std::cout << "\n\n ! " << con_v_c_parse_csv.size();
		return con_v_c_parse_csv;
	}

	std::vector<fpf_parse::parse_peptides_csv_type> parse_filesystem_denovopeptides(string_type par_str_fin_root) {
		std::vector<fpf_parse::parse_peptides_csv_type> con_v_c_parse_csv;
		std::ifstream fin_input_csv(par_str_fin_root);
		if (IgFamily::NOVOR_DENOVO) {
			con_v_c_parse_csv = fpf_parse::parse_NOVOR_denovopeptides_csv(fin_input_csv, par_str_fin_root);
		}
		else {
			con_v_c_parse_csv = fpf_parse::parse_PEAKS_denovopeptides_csv(fin_input_csv, par_str_fin_root);
		}
		std::cout << "\n\n ! " << con_v_c_parse_csv.size();
		return con_v_c_parse_csv;
	}

	void create_v_p_replicate_data(filesystem_type& par_s_filesystem) {
		for (std::vector<peptide_data_type>::iterator itr_v_s_peptide_data = par_s_filesystem.v_s_peptide_data.begin(); itr_v_s_peptide_data != par_s_filesystem.v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
			itr_v_s_peptide_data->v_p_peptideassociation.clear();
			itr_v_s_peptide_data->v_p_peptideassociation_distinct.clear();
			for (size_type itr_st_replicate_count = 0; itr_st_replicate_count < par_s_filesystem.st_replicate_count; ++itr_st_replicate_count) {
				bool b_peptide_replicate_found = bool();
				size_type st_peptide_replicate_spectralcount = size_type();
				size_type st_peptide_replicate_IgP = size_type();
				for (auto itr_v_s_peptide_data_2 = par_s_filesystem.v_s_peptide_data.begin(); itr_v_s_peptide_data_2 != par_s_filesystem.v_s_peptide_data.end(); ++itr_v_s_peptide_data_2) {
					if ((itr_v_s_peptide_data_2->str_peptide) == (itr_v_s_peptide_data->str_peptide) && (itr_v_s_peptide_data_2->st_filesystem_replicate == (itr_st_replicate_count + size_type{ 1 })) && (!itr_v_s_peptide_data_2->b_replicate_merged)) {
						b_peptide_replicate_found = true;
						st_peptide_replicate_spectralcount = itr_v_s_peptide_data_2->st_spectralcount;
						st_peptide_replicate_IgP = itr_v_s_peptide_data_2->st_IgP;
						if (!itr_v_s_peptide_data->v_p_replicate_data.empty()) {
							itr_v_s_peptide_data_2->b_replicate_merged = true;
						}
						for (auto i = size_type(itr_v_s_peptide_data->v_p_replicate_data.size()); i < (itr_v_s_peptide_data_2->st_filesystem_replicate); ++i) {
							if ((itr_v_s_peptide_data->v_p_replicate_data.size() + size_type(1)) == itr_v_s_peptide_data_2->st_filesystem_replicate) {
								itr_v_s_peptide_data->v_p_replicate_data.push_back(std::make_tuple(par_s_filesystem.str_filename, st_peptide_replicate_spectralcount, st_peptide_replicate_IgP));
								break;
							}
							itr_v_s_peptide_data->v_p_replicate_data.push_back(std::make_tuple(string_type(), size_type(), size_type()));
							if (i > par_s_filesystem.st_replicate_count) {
								std::cout << "\n\n\nERROR: replicate count exceeded - iterator out of range";
								string_type test;
								std::cin >> test;
								break;
							}
						}
					}
				}
				if (itr_v_s_peptide_data->v_p_replicate_data.size() != 0) {
					itr_v_s_peptide_data->st_spectralcount = (size_type(((itr_v_s_peptide_data->v_p_replicate_data.size() * itr_v_s_peptide_data->st_spectralcount) + st_peptide_replicate_spectralcount) / double(itr_v_s_peptide_data->v_p_replicate_data.size() + 1)));
					itr_v_s_peptide_data->st_IgP = (size_type(((itr_v_s_peptide_data->v_p_replicate_data.size() * itr_v_s_peptide_data->st_IgP) + st_peptide_replicate_IgP) / double(itr_v_s_peptide_data->v_p_replicate_data.size() + 1)));
				}
			}
			for (auto i = size_type(itr_v_s_peptide_data->v_p_replicate_data.size()); i < par_s_filesystem.st_replicate_count; ++i) {
				itr_v_s_peptide_data->v_p_replicate_data.push_back(std::make_tuple(string_type(), size_type(), size_type()));
			}
		}
		for (std::vector<peptide_data_type>::iterator itr_v_s_peptide_data = par_s_filesystem.v_s_peptide_data.begin(); itr_v_s_peptide_data != par_s_filesystem.v_s_peptide_data.end(); itr_v_s_peptide_data->b_replicate_merged ? itr_v_s_peptide_data : ++itr_v_s_peptide_data) {
			if (itr_v_s_peptide_data->b_replicate_merged) {
				par_s_filesystem.v_s_peptide_data.erase(itr_v_s_peptide_data);
			}
		}
	}	

	std::vector<string_type> create_v_str_patientstatus(std::vector<filesystem_type> par_v_s_filesystem) {
		std::vector<string_type> con_v_str_patientstatus;
		for (auto itr_v_s_filesystem = par_v_s_filesystem.begin(); itr_v_s_filesystem != par_v_s_filesystem.end(); ++itr_v_s_filesystem) {
			auto find_str_patientstatus = std::find(con_v_str_patientstatus.begin(), con_v_str_patientstatus.end(), itr_v_s_filesystem->str_patientstatus);
			if (find_str_patientstatus == con_v_str_patientstatus.end()) {
				con_v_str_patientstatus.push_back(itr_v_s_filesystem->str_patientstatus);
			}
		}
		return con_v_str_patientstatus;
	}

	void fout_filesystem(filesystem_type par_s_filesystem) {
		std::string output_v_c_analysis = par_s_filesystem.str_directory + "filesystem.data";
		std::ofstream fout_v_s_filesystem;
		fout_v_s_filesystem.open(output_v_c_analysis);
		fout_v_s_filesystem << "ID: " << std::get<0>(par_s_filesystem.p_filesystemid) << "," << std::get<1>(par_s_filesystem.p_filesystemid) << ";\n";
		fout_v_s_filesystem << "FILE: " << par_s_filesystem.str_filename << ";\n";
		fout_v_s_filesystem << "VERSION: " << IgFamily::version << ";\n";
		fout_v_s_filesystem << "REPLICATES: ";
		for (std::vector<std::pair<string_type, string_type>>::iterator itr_v_p_replicates = par_s_filesystem.v_p_replicates.begin(); itr_v_p_replicates != par_s_filesystem.v_p_replicates.end(); ++itr_v_p_replicates) {
			fout_v_s_filesystem << std::get<0>(*itr_v_p_replicates) << "," << std::get<1>(*itr_v_p_replicates);
			if ((itr_v_p_replicates + 1) != par_s_filesystem.v_p_replicates.end()) {
				fout_v_s_filesystem << ",";
			}
		}
		fout_v_s_filesystem << ";\n";
		fout_v_s_filesystem << "STATUS: " << par_s_filesystem.str_patientstatus << ";\n";
		fout_v_s_filesystem << "ENZYME: " << par_s_filesystem.str_enzyme << ";\n";
		//if (par_s_filesystem.b_wiffconversion) {
		//	fout_v_s_filesystem << "CHARGE_STATE_PREDICTOR: ";
		//	if (par_s_filesystem.str_chargestatepredictor == "") {
		//		fout_v_s_filesystem << par_s_filesystem.str_chargestatepredictor << ";\n";
		//	}
		//}
		fout_v_s_filesystem << "DENOVO_DELTAMASS: " << par_s_filesystem.str_denono_deltamass << ";\n";
	}
}

#endif