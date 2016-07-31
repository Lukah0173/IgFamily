// * * fpf_parse.h * *
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#include <cstdlib> // provides - size_t
#include <string> // provides - std::string
#include <iostream> // provides - std::cin, std::cout
#include <fstream> // provides - std::ifstream, std::ofstream
#include <istream> // provides - std::istream::get
#include <vector> // provides - std::vector
#include <map> // provides - std::map
#include <set>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "IgFamily.h"

#ifndef FPF_PARSE
#define	FPF_PARSE

namespace fpf_parse {

	struct s_parse_peptides_csv;
	class c_parse_FASTA;

	typedef size_t size_type;
	typedef std::string string_type;
	typedef std::vector<s_parse_peptides_csv> v_s_parse_csv_data_type;

	struct s_parse_peptides_csv {
		s_parse_peptides_csv() {
		};

		~s_parse_peptides_csv() {
		};

	public:
		string_type str_parse_peptides_csv_file;
		string_type str_parse_peptides_csv_peptide;
		string_type str_parse_peptides_csv_spectralcount;
		string_type str_parse_csv_IgP;
		std::vector<double> v_d_denovo_localconfidence;
	};

	class c_parse_FASTA {
	public:
		c_parse_FASTA() {
		};

		~c_parse_FASTA() {
		};

		inline void set_str_parse_FASTA_accession(string_type par_str_parse_FASTA_accession) {
			str_parse_FASTA_accession = par_str_parse_FASTA_accession;
		};

		inline void set_str_parse_FASTA_protein_delimited(string_type par_str_parse_FASTA_protein_delimited) {
			str_parse_FASTA_protein_delimited = par_str_parse_FASTA_protein_delimited;
		};

		inline void set_str_parse_FASTA_genefamily(string_type par_str_parse_FASTA_genefamily) {
			str_parse_FASTA_genefamily = par_str_parse_FASTA_genefamily;
		};

		inline void set_str_parse_FASTA_genefamily_class(string_type par_str_parse_FASTA_genefamily_class) {
			str_parse_FASTA_genefamily_class = par_str_parse_FASTA_genefamily_class;
		};

		inline void set_str_parse_FASTA_species(string_type par_str_parse_FASTA_species) {
			str_parse_FASTA_species = par_str_parse_FASTA_species;
		};

		inline void set_str_parse_FASTA_protein(string_type par_str_parse_FASTA_protein) {
			str_parse_FASTA_protein = par_str_parse_FASTA_protein;
		};

		inline const string_type return_str_parse_FASTA_accession() const {
			return str_parse_FASTA_accession;
		};

		inline const string_type return_str_parse_FASTA_protein_delimited() const {
			return str_parse_FASTA_protein_delimited;
		};

		inline const string_type return_str_parse_FASTA_genefamily() const {
			return str_parse_FASTA_genefamily;
		};

		inline const string_type return_str_parse_FASTA_genefamily_class() const {
			return str_parse_FASTA_genefamily_class;
		};

		inline const string_type return_str_parse_FASTA_species() const {
			return str_parse_FASTA_species;
		};

		inline const string_type return_str_parse_FASTA_protein() const {
			return str_parse_FASTA_protein;
		};

	private:
		string_type str_parse_FASTA_accession;
		string_type str_parse_FASTA_protein_delimited;
		string_type str_parse_FASTA_genefamily;
		string_type str_parse_FASTA_genefamily_class;
		string_type str_parse_FASTA_species;
		string_type str_parse_FASTA_protein;
	};

	std::vector<s_parse_peptides_csv> parse_proteinpeptides(std::ifstream& par_fin_proteinpeptides_csv, string_type par_str_dir) {

		std::vector<s_parse_peptides_csv> con_v_c_parse_proteinpeptides_csv;
		string_type str_parse_csv;
		size_type st_count_csv = size_type();
		size_type sw_parse_csv = size_type();
		size_type sw_parse_csv_ignore_header = size_type();
		char ch_parse_csv;
		string_type str_parse_peptide;
		string_type str_parse_spectralcount;
		string_type str_parse_IgP;

		while (par_fin_proteinpeptides_csv.std::istream::get(ch_parse_csv)) {
			if (sw_parse_csv_ignore_header == 1) {
				str_parse_csv = ch_parse_csv;
				if (st_count_csv == 18) {
					if (con_v_c_parse_proteinpeptides_csv.size() == 0) {
						s_parse_peptides_csv con_c_parse_csv = s_parse_peptides_csv();
						con_c_parse_csv.str_parse_peptides_csv_file = par_str_dir;
						con_c_parse_csv.str_parse_peptides_csv_peptide = str_parse_peptide;
						con_c_parse_csv.str_parse_peptides_csv_spectralcount = str_parse_spectralcount;
						con_c_parse_csv.str_parse_csv_IgP= str_parse_IgP;
						con_v_c_parse_proteinpeptides_csv.push_back(con_c_parse_csv);
					}
					else {
						size_type count_v_c_parse_csv = size_type();
						for (auto& itr_v_c_parse_csv : con_v_c_parse_proteinpeptides_csv) {
							++count_v_c_parse_csv;
							if (itr_v_c_parse_csv.str_parse_peptides_csv_peptide == str_parse_peptide) {
								break;
							}
							if (count_v_c_parse_csv == con_v_c_parse_proteinpeptides_csv.size()) {
								s_parse_peptides_csv con_c_parse_csv = s_parse_peptides_csv();
								con_c_parse_csv.str_parse_peptides_csv_file = par_str_dir;
								con_c_parse_csv.str_parse_peptides_csv_peptide = str_parse_peptide;
								con_c_parse_csv.str_parse_peptides_csv_spectralcount = str_parse_spectralcount;
								con_c_parse_csv.str_parse_csv_IgP = str_parse_IgP;
								con_v_c_parse_proteinpeptides_csv.push_back(con_c_parse_csv);

								break;
							}
						}
					}
					str_parse_peptide.clear();
					str_parse_spectralcount.clear();
					str_parse_IgP.clear();
					st_count_csv = 0;
				}
				if (str_parse_csv == ",") {
					++st_count_csv;
				}
				if ((st_count_csv == 3) && (str_parse_csv != ",")) {
					if (str_parse_csv == "(") {
						sw_parse_csv = 2;
					}
					if (str_parse_csv == ")") {
						sw_parse_csv = 1;
					}
					if ((sw_parse_csv == 1) && ((str_parse_csv == ".") && (str_parse_peptide.length() > 2))) {
						sw_parse_csv = 0;
					}
					if ((sw_parse_csv == 1) || (sw_parse_csv == 2)) {
						str_parse_peptide += str_parse_csv;
					}
					if ((str_parse_csv == ".") && (str_parse_peptide.length() <= 2)) {
						str_parse_peptide.clear();
						sw_parse_csv = 1;
					}
					if ((sw_parse_csv != 2) && (str_parse_csv != ".")) {
						sw_parse_csv = 1;
					}
				}
				if ((st_count_csv == 4) && (str_parse_csv != ",")) {
					str_parse_csv.clear();
				}
				if ((st_count_csv == 5) && (str_parse_csv != ",")) {
					str_parse_IgP += str_parse_csv;
				}
				if ((st_count_csv == 14) && (str_parse_csv != ",")) {
					str_parse_spectralcount += str_parse_csv;
				}
			}
			if ((sw_parse_csv_ignore_header == 0) && (ch_parse_csv == ',')) {
				++st_count_csv;
				if (st_count_csv == 18) {
					sw_parse_csv_ignore_header = 1;
					st_count_csv = 0;
				}
			}
		}

		par_fin_proteinpeptides_csv.clear();
		par_fin_proteinpeptides_csv.seekg(0, std::ios::beg);

		return con_v_c_parse_proteinpeptides_csv;
	}

	std::vector<s_parse_peptides_csv> parse_denovopeptides_csv(std::ifstream& par_fin_denovopeptides_csv, string_type par_str_dir) {

		std::vector<s_parse_peptides_csv> con_v_c_parse_peptides_csv;
		string_type str_parse_csv;
		size_type st_count_csv = size_type();
		size_type sw_parse_csv = size_type();
		size_type sw_parse_csv_ignore_header = size_type();
		char ch_parse_csv;
		string_type str_parse_peptide;
		string_type con_str_denovo_localconfidence;
		std::vector<double> con_v_d_denovo_localconfidence;

		while (par_fin_denovopeptides_csv.std::istream::get(ch_parse_csv)) {
			if (sw_parse_csv_ignore_header == 1) {
				str_parse_csv = ch_parse_csv;
				if (st_count_csv == 16) {
					if (con_v_c_parse_peptides_csv.size() == 0) {
						s_parse_peptides_csv con_c_parse_csv = s_parse_peptides_csv();
						con_c_parse_csv.str_parse_peptides_csv_file = par_str_dir;
						con_c_parse_csv.str_parse_peptides_csv_peptide = str_parse_peptide;
						con_c_parse_csv.v_d_denovo_localconfidence = con_v_d_denovo_localconfidence;
						con_v_c_parse_peptides_csv.push_back(con_c_parse_csv);
						con_v_d_denovo_localconfidence.clear();
					}
					else {
						bool b_parse_peptide_found = bool();
						for (auto itr_v_c_parse_csv = con_v_c_parse_peptides_csv.begin(); itr_v_c_parse_csv != con_v_c_parse_peptides_csv.end(); ++itr_v_c_parse_csv) {
							if (itr_v_c_parse_csv->str_parse_peptides_csv_peptide == str_parse_peptide) {
								b_parse_peptide_found = true;
							}
							if (((itr_v_c_parse_csv + 1) == con_v_c_parse_peptides_csv.end()) && (!b_parse_peptide_found)) {
								s_parse_peptides_csv con_c_parse_csv = s_parse_peptides_csv();
								con_c_parse_csv.str_parse_peptides_csv_file = par_str_dir;
								con_c_parse_csv.str_parse_peptides_csv_peptide = str_parse_peptide;
								con_c_parse_csv.v_d_denovo_localconfidence = con_v_d_denovo_localconfidence;
								con_v_c_parse_peptides_csv.push_back(con_c_parse_csv);
								con_v_d_denovo_localconfidence.clear();
								break;
							}
						}
					}
					str_parse_peptide.clear();
					st_count_csv = 0;
				}
				if (str_parse_csv == ",") {
					++st_count_csv;
				}
				if ((st_count_csv == 3) && (str_parse_csv != ",")) {
					if (str_parse_csv == "(") {
						sw_parse_csv = 2;
					}
					if (str_parse_csv == ")") {
						sw_parse_csv = 1;
					}
					if ((sw_parse_csv == 1) && ((str_parse_csv == ".") && (str_parse_peptide.length() > 2))) {
						sw_parse_csv = 0;
					}
					if ((sw_parse_csv == 1) || (sw_parse_csv == 2)) {
						str_parse_peptide += str_parse_csv;
					}
					if ((str_parse_csv == ".") && (str_parse_peptide.length() <= 2)) {
						str_parse_peptide.clear();
						sw_parse_csv = 1;
					}
					if ((sw_parse_csv != 2) && (str_parse_csv != ".")) {
						sw_parse_csv = 1;
					}
				}
				if ((st_count_csv == 4) && (str_parse_csv != ",")) {
					str_parse_csv.clear();
				}
				if ((st_count_csv == 14) && (ch_parse_csv != ',')) {
					if (ch_parse_csv == ' ') {
						con_v_d_denovo_localconfidence.push_back(std::stod(con_str_denovo_localconfidence));
						con_str_denovo_localconfidence.clear();
					}
					if (ch_parse_csv != ' ') {
						con_str_denovo_localconfidence += str_parse_csv;
					}
				}
				if ((st_count_csv == 15) && (ch_parse_csv == ',')) {
					con_v_d_denovo_localconfidence.push_back(std::stod(con_str_denovo_localconfidence));
					con_str_denovo_localconfidence.clear();
				}
			}
			if ((sw_parse_csv_ignore_header == 0) && (ch_parse_csv == ',')) {
				++st_count_csv;
				if (st_count_csv == 16) {
					sw_parse_csv_ignore_header = 1;
					st_count_csv = 0;
				}
			}
		}

		par_fin_denovopeptides_csv.clear();
		par_fin_denovopeptides_csv.seekg(0, std::ios::beg);

		return con_v_c_parse_peptides_csv;
	}

	std::vector<c_parse_FASTA> parse_FASTA(std::ifstream& par_fin_input_FASTA) {

		char ch_parse_FASTA;
		size_type sw_input_FASTA = size_type();
		size_type sw_2_input_FASTA = size_type();
		string_type con_str_FASTA_genefamily = string_type();
		string_type con_str_FASTA_genefamily_class = string_type();
		string_type con_str_FASTA_protein = string_type();
		string_type con_str_parse_FASTA_accession = string_type();
		string_type con_str_FASTA_species = string_type();
		string_type con_str_parse_FASTA_protein_delimited = string_type();
		c_parse_FASTA con_c_parse_FASTA = c_parse_FASTA();
		std::vector<c_parse_FASTA> con_v_c_parse_FASTA;

		while (par_fin_input_FASTA.get(ch_parse_FASTA)) {
			if (sw_2_input_FASTA == 1) {
				if (ch_parse_FASTA != '\n') {
					con_str_FASTA_protein += ch_parse_FASTA;
				}
				con_str_parse_FASTA_protein_delimited += ch_parse_FASTA;
			}
			if ((sw_2_input_FASTA == 0) && (ch_parse_FASTA == '\n')) {
				con_str_FASTA_protein.clear();
				sw_2_input_FASTA = 1;
			}
			if ((sw_input_FASTA == 2) && (ch_parse_FASTA == '|')) {
				sw_input_FASTA = 3;
			}
			if (sw_input_FASTA == 2) {
				con_str_FASTA_species += ch_parse_FASTA;
			}
			if ((sw_input_FASTA == 1) && (ch_parse_FASTA == '|')) {
				sw_input_FASTA = 2;
			}
			if (sw_input_FASTA == 1) {
				if (ch_parse_FASTA != '_') {
					if (con_str_FASTA_genefamily == "IGHV") {
						con_str_FASTA_genefamily_class = con_str_FASTA_genefamily;
					}
					if (con_str_FASTA_genefamily == "IGLV") {
						con_str_FASTA_genefamily_class = con_str_FASTA_genefamily;
					}
					con_str_FASTA_genefamily += ch_parse_FASTA;
				}
				else {
					con_str_FASTA_genefamily += '*';
				}
			}
			if ((sw_input_FASTA == 0) && (ch_parse_FASTA == '|')) {
				con_str_FASTA_genefamily.clear();
				sw_input_FASTA = 1;
			}
			if ((sw_2_input_FASTA == 1) && (par_fin_input_FASTA.peek() == '>')) {
				con_c_parse_FASTA.set_str_parse_FASTA_accession("test");
				con_c_parse_FASTA.set_str_parse_FASTA_protein_delimited(con_str_parse_FASTA_protein_delimited);
				con_c_parse_FASTA.set_str_parse_FASTA_genefamily(con_str_FASTA_genefamily);
				con_c_parse_FASTA.set_str_parse_FASTA_genefamily_class(con_str_FASTA_genefamily_class);
				con_c_parse_FASTA.set_str_parse_FASTA_species(con_str_FASTA_species);
				con_c_parse_FASTA.set_str_parse_FASTA_protein(con_str_FASTA_protein);
				con_v_c_parse_FASTA.push_back(con_c_parse_FASTA);
				con_str_parse_FASTA_accession.clear();
				con_str_parse_FASTA_protein_delimited.clear();
				con_str_FASTA_protein.clear();
				con_str_FASTA_genefamily.clear();
				con_str_FASTA_genefamily_class.clear();
				con_str_FASTA_species.clear();
				sw_input_FASTA = 0;
				sw_2_input_FASTA = 0;
			}
		}

		par_fin_input_FASTA.clear();
		par_fin_input_FASTA.seekg(0, std::ios::beg);

		return con_v_c_parse_FASTA;
	}

	bool check_protein_peptides(std::vector<fpf_parse::s_parse_peptides_csv> par_v_c_parse_csv_proteinpeptides_data, bool par_filesystem_modified) {
		if (par_v_c_parse_csv_proteinpeptides_data.empty()) {
			if (par_filesystem_modified) {
				std::cout << "\n\n --- the protein_peptides file is empty or does not exist";
			}
			return false;
		}
		if (!par_v_c_parse_csv_proteinpeptides_data.empty()) {
			if (par_filesystem_modified) {
				std::cout << "\n\n --- protein_peptides file found";
			}
			return true;
		}
		return false;
	}

	bool check_denovo_peptides(std::vector<fpf_parse::s_parse_peptides_csv> par_v_c_parse_csv_denovopeptides_data, bool par_filesystem_modified) {
		if (par_v_c_parse_csv_denovopeptides_data.empty()) {
			if (par_filesystem_modified) {
				std::cout << "\n\n --- the denovo_peptides file is empty or does not exist";
			}
			return false;
		}
		if (!par_v_c_parse_csv_denovopeptides_data.empty()) {
			if (par_filesystem_modified) {
				std::cout << "\n\n --- denovo_peptides file found";
			}
			return true;
		}
		return false;
	}

	//if ((par_v_c_parse_csv_proteinpeptides_data.empty()) && (par_v_c_parse_csv_denovopeptides_data.empty()) && (par_filesystem_modified)) {
//	std::cout << "\n\n\n input any key to exit...\n\n > ";
//	std::string s_catch_error;
//	std::cin >> s_catch_error;
//	return true;
//}

	bool b_parse_FASTA_empty(std::vector<c_parse_FASTA> par_v_c_parse_FASTA) {
		if (par_v_c_parse_FASTA.size() == 0) {
			std::cout << "\n\n * * * FASTA file empty..\n\n * * * Is the file correctly directed?";
			std::cout << "\n\n The program will now terminate. Input any key to continue -\n\n -> ";
			string_type str_exit;
			std::cin >> str_exit;
			return true;
		}
		return false;
	}

	void output_v_c_parse_FASTA(std::vector<c_parse_FASTA> par_v_c_parse_FASTA) {
		string_type output_FASTA_filtered = "output.fasta";
		std::ofstream fout_FASTA_filtered;
		fout_FASTA_filtered.open(output_FASTA_filtered);
		for (auto itr_v_c_parse_FASTA : par_v_c_parse_FASTA) {
			if ((IgFamily::OUTPUT_FASTA == 1) && ((itr_v_c_parse_FASTA.return_str_parse_FASTA_accession().find("IGHV3-23") != std::string::npos))) {
				fout_FASTA_filtered << itr_v_c_parse_FASTA.return_str_parse_FASTA_accession();
				fout_FASTA_filtered << itr_v_c_parse_FASTA.return_str_parse_FASTA_protein_delimited();
			}
		}
	}

	//void clear_v_s_parse_csv_data

	//void IgSplicer() {
	//	string_type CMD_SPLICER = "-spl";
	//	if (CMD_SPLICER == "-spl") {
	//		string_type input_data;
	//		int debug_mode_spl = 1;

	//		std::cout << "Input?\n\n";
	//		std::cin >> input_data;
	//		input_data = input_data + ".txt";

	//		std::ifstream fin_input_data(input_data);
	//		char constructor_input_data;
	//		string_type string_input_data;
	//		string_type condition_input_data;
	//		string_type constructor_IGHJ;
	//		string_type constructor_IGHC;
	//		std::map <string_type, string_type> map_IGHJ;
	//		std::map <string_type, string_type> map_IGHC;
	//		int switch_csv = 0;
	//		int switch_JC = 0;

	//		while (fin_input_data.get(constructor_input_data)) {
	//			condition_input_data += constructor_input_data;
	//			if ((condition_input_data != ",") && (condition_input_data != "\n") && (condition_input_data != "\r\n") && (condition_input_data != "\r")) {
	//				string_input_data += constructor_input_data;
	//			}
	//			if (condition_input_data == "!") {
	//				switch_JC = 1;
	//				string_input_data.clear();
	//			}
	//			if (condition_input_data == ",") {
	//				if (switch_JC == 0) {
	//					if (switch_csv == 0) {
	//						constructor_IGHJ = string_input_data;
	//						string_input_data.clear();
	//						switch_csv = 1;
	//					}
	//					else {
	//						map_IGHJ[constructor_IGHJ] = string_input_data;
	//						string_input_data.clear();
	//						switch_csv = 0;
	//					}
	//				}
	//				if (switch_JC == 1) {
	//					if (switch_csv == 0) {
	//						constructor_IGHC = string_input_data;
	//						string_input_data.clear();
	//						switch_csv = 1;
	//					}
	//					else {
	//						map_IGHC[constructor_IGHC] = string_input_data;
	//						string_input_data.clear();
	//						switch_csv = 0;
	//					}
	//				}
	//			}
	//			condition_input_data.clear();
	//		}

	//		if (debug_mode_spl == 1) {
	//			string_type output_debug_map_IGHJ_IGHC = "z" + input_data + "_debug_map_IGHJ-IGHC.txt";
	//			std::ofstream fout_debug_map_IGHJ_IGHC;
	//			fout_debug_map_IGHJ_IGHC.open(output_debug_map_IGHJ_IGHC);
	//			fout_debug_map_IGHJ_IGHC << "-- IgSplicer " << version << " --\n\n\n" << "Input file: " << input_data << "\n\n\n";
	//			for (std::map<string_type, string_type>::const_iterator i = map_IGHJ.begin(); i != map_IGHJ.end(); ++i) {
	//				fout_debug_map_IGHJ_IGHC << i->first << " " << i->second << "\n\n";
	//			}
	//			for (std::map<string_type, string_type>::const_iterator i = map_IGHC.begin(); i != map_IGHC.end(); ++i) {
	//				fout_debug_map_IGHJ_IGHC << i->first << " " << i->second << "\n\n";
	//			}
	//		}

	//		int count_accession = 001;
	//		int count_accession_peptide = 0;
	//		string_type constructor_accession;
	//		string_type accession;
	//		string_type constructor_accession_peptide;

	//		string_type output_map_IGHJ_IGHC = input_data + "_output.txt";
	//		std::ofstream fout_map_IGHJ_IGHC;
	//		fout_map_IGHJ_IGHC.open(output_map_IGHJ_IGHC);
	//		fout_map_IGHJ_IGHC << "-- IgSplicer " << version << " --\n\n\n" << "Input file: " << input_data << "\n\n\n";

	//		std::map<string_type, string_type> map_IGHJ_IGHC;
	//		for (std::map<string_type, string_type>::const_iterator i = map_IGHJ.begin(); i != map_IGHJ.end(); ++i) {
	//			for (std::map<string_type, string_type>::const_iterator j = map_IGHC.begin(); j != map_IGHC.end(); ++j) {
	//				std::stringstream constructor_accession;
	//				constructor_accession << count_accession;
	//				accession = constructor_accession.str();
	//				fout_map_IGHJ_IGHC << ">F";
	//				for (unsigned k = 0; k < 6 - accession.length(); ++k) {
	//					fout_map_IGHJ_IGHC << "0";
	//				}
	//				fout_map_IGHJ_IGHC << accession + "|" + i->first + "-" + j->first + " \n";
	//				for (unsigned n1 = 0; n1 < (i->second + j->second).length(); ++n1) {
	//					constructor_accession_peptide += (i->second + j->second).at(n1);
	//					++count_accession_peptide;
	//					if (count_accession_peptide % 60 == 0) {
	//						constructor_accession_peptide += "\n";
	//					}
	//				}
	//				fout_map_IGHJ_IGHC << constructor_accession_peptide + "\n";
	//				constructor_accession_peptide.clear();
	//				count_accession_peptide = 0;
	//				++count_accession;
	//			}
	//		}
	//	}
	//}

}

#endif