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
#include "fpf_data.h"



namespace fpf_filesystem {

	struct s_filesystem;

	typedef std::string string_type;
	typedef size_t size_type;
	typedef fpf_data::s_multinomial_element_data s_multinomial_element_data;
	typedef fpf_data::s_peptide_data s_peptide_data;
	typedef fpf_data::s_blastp s_blastp_type;
	typedef fpf_data::s_mnom s_mnom_type;
	typedef fpf_data::s_blastp s_blastp_type;
	typedef fpf_data::s_report s_report_type;

	struct s_filesystem {
	public:
		bool b_proteinpeptides_exist;
		bool b_denovopeptides_exist;
		string_type str_directory;
		string_type str_filename;
		string_type str_fileversion;
		string_type str_patientstatus;
		std::vector<s_multinomial_element_data> v_c_analysis_data;
		std::vector<s_multinomial_element_data> v_c_analysis_distinct_data;
		std::vector<s_peptide_data> v_s_peptide_data;
		std::vector<s_peptide_data> v_s_peptide_data_filtered;
		std::vector<s_peptide_data> v_s_peptide_data_distinct;
		std::vector<s_peptide_data> v_s_peptide_data_filtered_distinct;
		std::pair<string_type, string_type> p_filesystemid;
		std::vector<std::pair<string_type, string_type>> v_p_replicates;
		size_type st_replicate_count;
		std::vector<s_blastp_type> v_s_blastp;
		std::vector<s_mnom_type> v_s_mnom;
		std::vector<s_report_type> v_s_report;
	};

	std::vector<string_type> read_root_dir(string_type par_IgFamily_root_dir) {
		std::ifstream fin_input_csv(par_IgFamily_root_dir);
		if (IgFamily::FILESYSTEM_MODE == 1) {
			std::cout << "\n\n ...reading root directory\n";
		}
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

	std::vector<s_filesystem> read_filesystem(std::vector<string_type> par_root_dir) {
		s_filesystem con_s_filesystem;
		std::vector<s_filesystem> con_v_s_filesystem;
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
						con_s_filesystem.str_directory = *itr_root_dir;
						con_s_filesystem.p_filesystemid = con_p_filesystem_id;
						con_s_filesystem.str_filename = con_str_filesystem_filename;
						con_s_filesystem.str_fileversion = con_str_filesystem_version;
						con_s_filesystem.str_patientstatus = con_str_filesystem_status;
						con_s_filesystem.v_p_replicates = v_p_filesystem_replicates;
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
		std::cout << con_v_s_filesystem.begin()->str_filename;
		return con_v_s_filesystem;
	}

	string_type read_filesystem_proteinpeptides(string_type par_IgFamily_root_dir) {
		string_type con_str_root_proteinpeptides = par_IgFamily_root_dir + "protein_peptides.csv";
		return con_str_root_proteinpeptides;
	}

	string_type read_filesystem_denovopeptides(string_type par_IgFamily_root_dir) {
		string_type con_str_root_denovopeptides = par_IgFamily_root_dir + "denovo_peptides.csv";
		return con_str_root_denovopeptides;
	}

	std::vector<fpf_parse::s_parse_peptides_csv> parse_filesystem_proteinpeptides(string_type par_str_fin_root) {
		std::vector<fpf_parse::s_parse_peptides_csv> con_v_c_parse_csv;
		std::ifstream fin_input_csv(par_str_fin_root);
		con_v_c_parse_csv = fpf_parse::parse_proteinpeptides(fin_input_csv, par_str_fin_root);
		std::cout << "\n\n ! " << con_v_c_parse_csv.size();
		return con_v_c_parse_csv;
	}

	std::vector<fpf_parse::s_parse_peptides_csv> parse_filesystem_denovopeptides(string_type par_str_fin_root) {
		std::vector<fpf_parse::s_parse_peptides_csv> con_v_c_parse_csv;
		std::ifstream fin_input_csv(par_str_fin_root);
		con_v_c_parse_csv = fpf_parse::parse_denovopeptides_csv(fin_input_csv, par_str_fin_root);
		std::cout << "\n\n ! " << con_v_c_parse_csv.size();
		return con_v_c_parse_csv;
	}

	void create_v_p_replicate_data(s_filesystem& par_s_filesystem) {
		for (std::vector<s_peptide_data>::iterator itr_v_s_peptide_data = par_s_filesystem.v_s_peptide_data.begin(); itr_v_s_peptide_data != par_s_filesystem.v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
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
		for (std::vector<s_peptide_data>::iterator itr_v_s_peptide_data = par_s_filesystem.v_s_peptide_data.begin(); itr_v_s_peptide_data != par_s_filesystem.v_s_peptide_data.end(); itr_v_s_peptide_data->b_replicate_merged ? itr_v_s_peptide_data : ++itr_v_s_peptide_data) {
			if (itr_v_s_peptide_data->b_replicate_merged) {
				par_s_filesystem.v_s_peptide_data.erase(itr_v_s_peptide_data);
			}
		}
	}	

	std::vector<string_type> create_v_str_patientstatus(std::vector<s_filesystem> par_v_s_filesystem) {
		std::vector<string_type> con_v_str_patientstatus;
		for (auto itr_v_s_filesystem = par_v_s_filesystem.begin(); itr_v_s_filesystem != par_v_s_filesystem.end(); ++itr_v_s_filesystem) {
			auto find_str_patientstatus = std::find(con_v_str_patientstatus.begin(), con_v_str_patientstatus.end(), itr_v_s_filesystem->str_patientstatus);
			if (find_str_patientstatus == con_v_str_patientstatus.end()) {
				con_v_str_patientstatus.push_back(itr_v_s_filesystem->str_patientstatus);
			}
		}
		return con_v_str_patientstatus;
	}

	void fout_filesystem(s_filesystem par_s_filesystem) {
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
		fout_v_s_filesystem << "STATUS: " << par_s_filesystem.str_patientstatus << ";";
	}
}

#endif