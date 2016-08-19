// * * fpf_parse.h * *
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#include <cstdlib> // provides - size_t
#include <string> // provides - string
#include <iostream> // provides - std::cin, std::cout
#include <fstream> // provides - std::ifstream, std::ofstream
#include <istream> // provides - std::istream::get
#include <vector> // provides - vector
#include <map> // provides - std::map
#include <set>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "IgFamily.h"

#ifndef FPF_PARSE
#define	FPF_PARSE

namespace fpf_parse {

	using std::string;
	using std::vector;

	struct parse_peptides_csv_type;
	class parse_FASTA_type;

	typedef vector<parse_peptides_csv_type> v_s_parse_csv_data_type;

	struct parse_peptides_csv_type {
	public:
		string str_parse_peptides_csv_file;
		string str_parse_peptides_csv_peptide;
		string str_parse_peptides_csv_spectralcount;
		string str_parse_csv_IgP;
		vector<double> v_d_denovo_localconfidence;
	};

	class parse_FASTA_type {
	public:
		parse_FASTA_type() {
		};

		~parse_FASTA_type() {
		};

		inline void set_str_parse_FASTA_accession(string par_str_parse_FASTA_accession) {
			str_parse_FASTA_accession = par_str_parse_FASTA_accession;
		};

		inline void set_str_parse_FASTA_genefamily(string par_str_parse_FASTA_genefamily) {
			str_parse_FASTA_genefamily = par_str_parse_FASTA_genefamily;
		};

		inline void set_str_parse_FASTA_genefamily_class(string par_str_parse_FASTA_genefamily_class) {
			str_parse_FASTA_genefamily_class = par_str_parse_FASTA_genefamily_class;
		};

		inline void set_str_parse_FASTA_species(string par_str_parse_FASTA_species) {
			str_parse_FASTA_species = par_str_parse_FASTA_species;
		};

		inline void set_str_parse_FASTA_protein(string par_str_parse_FASTA_protein) {
			str_parse_FASTA_protein = par_str_parse_FASTA_protein;
		};

		inline const string return_str_parse_FASTA_accession() const {
			return str_parse_FASTA_accession;
		};

		inline const string return_str_parse_FASTA_genefamily() const {
			return str_parse_FASTA_genefamily;
		};

		inline const string return_str_parse_FASTA_genefamily_class() const {
			return str_parse_FASTA_genefamily_class;
		};

		inline const string return_str_parse_FASTA_species() const {
			return str_parse_FASTA_species;
		};

		inline const string return_str_parse_FASTA_protein() const {
			return str_parse_FASTA_protein;
		};

	private:
		string str_parse_FASTA_accession;
		string str_parse_FASTA_protein_delimited;
		string str_parse_FASTA_genefamily;
		string str_parse_FASTA_genefamily_class;
		string str_parse_FASTA_species;
		string str_parse_FASTA_protein;
	};

	vector<parse_peptides_csv_type> parse_proteinpeptides(std::ifstream& par_fin_proteinpeptides_csv, string par_str_dir) {

		vector<parse_peptides_csv_type> con_v_c_parse_proteinpeptides_csv;
		string str_parse_csv;
		size_t st_count_csv = size_t();
		size_t sw_parse_csv = size_t();
		size_t sw_parse_csv_ignore_header = size_t();
		char ch_parse_csv;
		string str_parse_peptide;
		string str_parse_spectralcount;
		string str_parse_IgP;

		while (par_fin_proteinpeptides_csv.std::istream::get(ch_parse_csv)) {
			if (sw_parse_csv_ignore_header == 1) {
				str_parse_csv = ch_parse_csv;
				if (st_count_csv == 18) {
					if (con_v_c_parse_proteinpeptides_csv.size() == 0) {
						parse_peptides_csv_type con_c_parse_csv = parse_peptides_csv_type();
						con_c_parse_csv.str_parse_peptides_csv_file = par_str_dir;
						con_c_parse_csv.str_parse_peptides_csv_peptide = str_parse_peptide;
						con_c_parse_csv.str_parse_peptides_csv_spectralcount = str_parse_spectralcount;
						con_c_parse_csv.str_parse_csv_IgP = str_parse_IgP;
						con_v_c_parse_proteinpeptides_csv.push_back(con_c_parse_csv);
					}
					else {
						size_t count_v_c_parse_csv = size_t();
						for (auto& itr_v_c_parse_csv : con_v_c_parse_proteinpeptides_csv) {
							++count_v_c_parse_csv;
							if (itr_v_c_parse_csv.str_parse_peptides_csv_peptide == str_parse_peptide) {
								break;
							}
							if (count_v_c_parse_csv == con_v_c_parse_proteinpeptides_csv.size()) {
								parse_peptides_csv_type con_c_parse_csv = parse_peptides_csv_type();
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

	vector<parse_peptides_csv_type> parse_PEAKS_denovopeptides_csv(std::ifstream& par_fin_denovopeptides_csv, string par_str_dir) {

		vector<parse_peptides_csv_type> con_v_c_parse_peptides_csv;
		string str_parse_csv;
		size_t st_count_csv = size_t();
		size_t sw_parse_csv = size_t(1);
		size_t sw_parse_csv_ignore_header = size_t();
		char ch_parse_csv;
		string str_parse_peptide;
		string con_str_denovo_localconfidence;
		vector<double> con_v_d_denovo_localconfidence;

		while (par_fin_denovopeptides_csv.std::istream::get(ch_parse_csv)) {
			if (sw_parse_csv_ignore_header == 1) {
				str_parse_csv = ch_parse_csv;
				if (st_count_csv == 16) {
					parse_peptides_csv_type con_c_parse_csv = parse_peptides_csv_type();
					con_c_parse_csv.str_parse_peptides_csv_file = par_str_dir;
					con_c_parse_csv.str_parse_peptides_csv_peptide = str_parse_peptide;
					con_c_parse_csv.v_d_denovo_localconfidence = con_v_d_denovo_localconfidence;
					con_v_c_parse_peptides_csv.push_back(con_c_parse_csv);
					con_v_d_denovo_localconfidence.clear();
					str_parse_peptide.clear();
					st_count_csv = 0;
				}
				if (str_parse_csv == ",") {
					++st_count_csv;
				}
				if ((st_count_csv == 3) && (str_parse_csv != ",") && (str_parse_csv != " ")) {
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

	vector<parse_peptides_csv_type> parse_NOVOR_denovopeptides_csv(std::ifstream& par_fin_denovopeptides_csv, string par_str_dir) {

		vector<parse_peptides_csv_type> con_v_c_parse_peptides_csv;
		string str_parse_csv;
		size_t st_count_csv = size_t();
		size_t sw_parse_csv = size_t(1);
		bool b_parse_csv_ignore_header = bool();
		string str_parse_csv_ignore_header = string();
		char ch_parse_csv;
		string str_parse_peptide;
		string con_str_denovo_localconfidence;
		vector<double> con_v_d_denovo_localconfidence;

		while (par_fin_denovopeptides_csv.std::istream::get(ch_parse_csv)) {
			if (b_parse_csv_ignore_header) {
				str_parse_csv = ch_parse_csv;
				if (str_parse_csv == ",") {
					++st_count_csv;
				}
				if ((st_count_csv == 11) && (ch_parse_csv == ',')) {
					con_v_d_denovo_localconfidence.push_back(std::stod(con_str_denovo_localconfidence));
					parse_peptides_csv_type con_c_parse_csv = parse_peptides_csv_type();
					con_c_parse_csv.str_parse_peptides_csv_file = par_str_dir;
					con_c_parse_csv.str_parse_peptides_csv_peptide = str_parse_peptide;
					con_c_parse_csv.v_d_denovo_localconfidence = con_v_d_denovo_localconfidence;
					con_v_c_parse_peptides_csv.push_back(con_c_parse_csv);
					con_v_d_denovo_localconfidence.clear();
					str_parse_peptide.clear();
					st_count_csv = 1;
				}
				if ((st_count_csv == 10) && (ch_parse_csv != ',')) {
					if (ch_parse_csv == '-') {
						con_v_d_denovo_localconfidence.push_back(std::stod(con_str_denovo_localconfidence));
						con_str_denovo_localconfidence.clear();
					}
					if ((ch_parse_csv != ' ') && (ch_parse_csv != '-')) {
						con_str_denovo_localconfidence += str_parse_csv;
					}
				}
				if ((st_count_csv == 9) && (str_parse_csv != ",") && (str_parse_csv != " ")) {
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
				if ((st_count_csv == 10) && (str_parse_csv != ",")) {
					str_parse_csv.clear();
				}
			}
			if (!b_parse_csv_ignore_header) {
				str_parse_csv_ignore_header += ch_parse_csv;
				if (str_parse_csv_ignore_header == " aaScore,") {
					b_parse_csv_ignore_header = true;
				}
				if (ch_parse_csv == ',') {
					str_parse_csv_ignore_header.clear();
				}
			}
		}

		par_fin_denovopeptides_csv.clear();
		par_fin_denovopeptides_csv.seekg(0, std::ios::beg);

		return con_v_c_parse_peptides_csv;
	}

	void check_FASTA_format(string par_INPUT_FASTA) {
		std::ifstream fin_INPUT_FASTA(par_INPUT_FASTA);
		char ch_parse_FASTA = char();
		string str_parse_FASTA = string();
		size_t count_delimit = size_t();
		bool b_read_format = bool();
		bool b_header_line = bool();
		while (fin_INPUT_FASTA.get(ch_parse_FASTA)) {
			if (ch_parse_FASTA == '>') {
				b_header_line = true;
				if (b_read_format) {
					std::cout << "\n\n Accession - \n";
					std::cout << str_parse_FASTA;
					str_parse_FASTA.clear();
					string str_format_break = string();
					std::cout << "\n\n\n continue? ( y / n ) ";					
					while (str_format_break == "") {
						std::cout << "\n\n --> ";
						std::cin >> str_format_break;
					}
					if (str_format_break == "y") {
					}
					else {
						break;
					}
				}
				count_delimit = 0;
				b_read_format = true;
			}
			if (b_header_line) {
				if ((ch_parse_FASTA == '|') || (ch_parse_FASTA == '\n') || (ch_parse_FASTA == '\r\n') || (ch_parse_FASTA == '\r')) {
					std::cout << "\n Field # ";
					std::cout << count_delimit;
					std::cout << " - \"";
					std::cout << str_parse_FASTA;
					std::cout << "\"";
					++count_delimit;
					str_parse_FASTA.clear();
				}
				if ((ch_parse_FASTA != '>') && (ch_parse_FASTA != '|') && ((ch_parse_FASTA != '\n') && (ch_parse_FASTA != '\r\n') && (ch_parse_FASTA != '\r'))) {
					str_parse_FASTA += ch_parse_FASTA;
				}
				if ((ch_parse_FASTA == '\n') || (ch_parse_FASTA == '\r\n') || (ch_parse_FASTA == '\r')) {					
					b_header_line = false;
					str_parse_FASTA.clear();
				}
			}
			if (!b_header_line) {
				str_parse_FASTA += ch_parse_FASTA;
			}
		}
		fin_INPUT_FASTA.clear();
		fin_INPUT_FASTA.seekg(0, std::ios::beg);
	}

	void output_custom_FASTA_format(string par_INPUT_FASTA) {
		std::ifstream fin_INPUT_FASTA(par_INPUT_FASTA);
		char ch_parse_FASTA = char();
		string str_parse_FASTA = string();
		size_t count_delimit = size_t();
		bool b_read_format = bool();
		bool b_header_line = bool();
		while (fin_INPUT_FASTA.get(ch_parse_FASTA)) {
			if (ch_parse_FASTA == '>') {
				b_header_line = true;
				if (b_read_format) {
					std::cout << "\n\n Accession - \n";
					std::cout << str_parse_FASTA;
					str_parse_FASTA.clear();
					std::cout << "\n\n current FASTA output format: ";
					std::cout << "\n\n >";
					for (auto i = 0; i < count_delimit; ++i) {
						std::cout << "[Field #";
						std::cout << i;
						std::cout << "]|";
					}
					std::cout << "\n [Accession]";
					std::cout << "\n\n\n Create custom output format:           [C] ";
					std::cout << "\n Output FASTA with selected format:     [X] ";
					string str_menu_selection = string();
					while ((str_menu_selection != ("C")) && (str_menu_selection != ("X"))) {
						str_menu_selection.clear();
						std::cout << "\n\n Input selection: \n\n --> ";
						std::cin >> str_menu_selection;
					}
					if (str_menu_selection == "C") {
						while (true) {
							vector<string> v_str_field_value(count_delimit);
							for (auto i = 0; i < count_delimit; ++i) {
								v_str_field_value[i] = std::to_string(i);
								std::cout << "\n ";
								std::cout << i;
								std::cout << ":     Modify Field #";
								std::cout << i;
							}
							std::cout << "\n T:     Truncate field output";
							std::cout << "\n X:     Continue";
							string str_change_field_value = "-1";
							while (!((std::stoi(str_change_field_value) < count_delimit) && (std::stoi(str_change_field_value) >= 0))) {
								std::cout << "\n\n Input selection: \n\n --> ";
								std::cin >> str_change_field_value;
								if ((str_change_field_value != "T") && (str_change_field_value != "X")){
									try {
										auto try_catch = std::stoi(str_change_field_value);
										if (try_catch == 0) {
										}
									}
									catch (const std::invalid_argument&) {
										str_change_field_value = "-1";
									}
								} 
								else {
									break;
								}
							}
							if (str_change_field_value == "T") {
								size_t st_truncate_v_str_field_value = size_t();
								std::cout << "\n\n Truncate to how many fields?: \n\n --> ";
								while (!((st_truncate_v_str_field_value > 0) && (st_truncate_v_str_field_value < count_delimit))) {
									std::cin >> st_truncate_v_str_field_value;
								}
								v_str_field_value.resize(st_truncate_v_str_field_value);
								std::cout << "\n\n current FASTA output format: ";
								std::cout << "\n\n >";
								for (auto i = 0; i < v_str_field_value.size(); ++i) {
									std::cout << "[Field #";
									std::cout << v_str_field_value[i];
									std::cout << "]|";
								}
								std::cout << "\n [Accession]";
								std::cout << "\n ";
							}
							if (str_change_field_value == "X") {
								break;
							}
							if ((str_change_field_value != "T") && (str_change_field_value != "X")) {
								std::cout << "\n\n Modifying Field #";
								std::cout << str_change_field_value;
								std::cout << "\n";
								for (auto i = 0; i < count_delimit; ++i) {
									v_str_field_value[i] = std::to_string(i);
									std::cout << "\n ";
									std::cout << i;
									std::cout << ":     Assign Field #";
									std::cout << i;
								}
								std::cout << "\n ";
								std::cout << count_delimit;
								std::cout << ":     Custom string entry";
								std::cout << "\n ";
								size_t st_assign_field_value = -1;
								while (!((st_assign_field_value <= count_delimit) && (st_assign_field_value >= 0))) {
									std::cout << "\n\n Input selection: \n\n --> ";
									std::cin >> st_assign_field_value;
									if (std::cin.fail()) {
										std::cin.clear();
										std::cin.ignore(256, '\n');
										st_assign_field_value = -1;
									}
								}
								if (st_assign_field_value == count_delimit) {
									string str_assign_field_value;
									std::cout << "\n\n Input selection: \n\n --> ";
									std::cin >> str_assign_field_value;
									v_str_field_value[std::stoi(str_change_field_value)] = str_assign_field_value;
								}
								if (st_assign_field_value < count_delimit) {
									v_str_field_value[std::stoi(str_change_field_value)] = std::to_string(st_assign_field_value);
								}
								std::cout << "\n\n current FASTA output format: ";
								std::cout << "\n\n >";
								for (auto i = 0; i < v_str_field_value.size(); ++i) {
									std::cout << "[Field #";
									std::cout << v_str_field_value[i];
									std::cout << "]|";
								}
								std::cout << "\n [Accession]";
								std::cout << "\n ";
							}
						}
					}
					if (str_menu_selection == "X") {
					}
					break;
				}
				count_delimit = 0;
				b_read_format = true;
			}
			if (b_header_line) {
				if ((ch_parse_FASTA == '|') || (ch_parse_FASTA == '\n') || (ch_parse_FASTA == '\r\n') || (ch_parse_FASTA == '\r')) {
					std::cout << "\n Field # ";
					std::cout << count_delimit;
					std::cout << " - \"";
					std::cout << str_parse_FASTA;
					std::cout << "\"";
					++count_delimit;
					str_parse_FASTA.clear();
				}
				if ((ch_parse_FASTA != '>') && (ch_parse_FASTA != '|') && ((ch_parse_FASTA != '\n') && (ch_parse_FASTA != '\r\n') && (ch_parse_FASTA != '\r'))) {
					str_parse_FASTA += ch_parse_FASTA;
				}
				if ((ch_parse_FASTA == '\n') || (ch_parse_FASTA == '\r\n') || (ch_parse_FASTA == '\r')) {
					b_header_line = false;
					str_parse_FASTA.clear();
				}
			}
			if (!b_header_line) {
				str_parse_FASTA += ch_parse_FASTA;
			}
		}
		fin_INPUT_FASTA.clear();
		fin_INPUT_FASTA.seekg(0, std::ios::beg);
	}

	void custom_FASTA_output(string par_INPUT_FASTA) {
		std::ifstream fin_INPUT_FASTA(par_INPUT_FASTA);
		string output_custom_FASTA = "FASTA\\custom_FASTA.fasta";
		std::ofstream fout_custom_FASTA;
		fout_custom_FASTA.open(output_custom_FASTA);
		char ch_parse_FASTA = char();
		string str_parse_FASTA = string();
		size_t count_delimit = size_t();
		bool b_read_format = bool();
		bool b_header_line = bool();
		while (fin_INPUT_FASTA.get(ch_parse_FASTA)) {
			if (ch_parse_FASTA == '>') {
				b_header_line = true;
				if (b_read_format) {
					fout_custom_FASTA << "Homo sapiens|";
					fout_custom_FASTA << str_parse_FASTA;
					str_parse_FASTA.clear();
				}
				count_delimit = 0;
				b_read_format = true;
			}
			if (b_header_line) {
				if ((ch_parse_FASTA == '|') || (ch_parse_FASTA == '\n') || (ch_parse_FASTA == '\r\n') || (ch_parse_FASTA == '\r')) {
					if (count_delimit == 0) {
						fout_custom_FASTA << ">";
					}
					if (count_delimit == 1) {
						fout_custom_FASTA << str_parse_FASTA;
						fout_custom_FASTA << "|";
					}
					if (count_delimit == 2) {
						fout_custom_FASTA << "UNIPROT(";
						fout_custom_FASTA << str_parse_FASTA;
						fout_custom_FASTA << ")";
						fout_custom_FASTA << "|";
					}
					if (count_delimit == 3) {
						fout_custom_FASTA << "Homo sapiens";
						fout_custom_FASTA << "|";
					}
					++count_delimit;
					str_parse_FASTA.clear();
				}
				if ((ch_parse_FASTA != '>') && (ch_parse_FASTA != '|') && ((ch_parse_FASTA != '\n') && (ch_parse_FASTA != '\r\n') && (ch_parse_FASTA != '\r'))) {
					str_parse_FASTA += ch_parse_FASTA;
				}
				if ((ch_parse_FASTA == '\n') || (ch_parse_FASTA == '\r\n') || (ch_parse_FASTA == '\r')) {
					b_header_line = false;
					str_parse_FASTA.clear();
				}
			}
			if (!b_header_line) {
				str_parse_FASTA += ch_parse_FASTA;
			}
		}
		fin_INPUT_FASTA.clear();
		fin_INPUT_FASTA.seekg(0, std::ios::beg);
	}

	vector<parse_FASTA_type> parse_FASTA(std::ifstream& par_fin_input_FASTA) {

		char ch_parse_FASTA;
		size_t sw_input_FASTA = size_t();
		size_t sw_2_input_FASTA = size_t();
		string con_str_FASTA_genefamily = string();
		string con_str_FASTA_genefamily_class = string();
		string con_str_FASTA_protein = string();
		string con_str_parse_FASTA_accession = string();
		string con_str_FASTA_species = string();
		parse_FASTA_type con_c_parse_FASTA = parse_FASTA_type();
		vector<parse_FASTA_type> con_v_c_parse_FASTA;
		size_t st_count_parse_FASTA = size_t();

		while (par_fin_input_FASTA.get(ch_parse_FASTA)) {
			if (sw_2_input_FASTA == 1) {
				if (ch_parse_FASTA != '\n') {
					con_str_FASTA_protein += ch_parse_FASTA;
				}
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
					if (con_str_FASTA_genefamily == "IGKV") {
						con_str_FASTA_genefamily_class = con_str_FASTA_genefamily;
					}
					if (con_str_FASTA_genefamily == "IGLV") {
						con_str_FASTA_genefamily_class = con_str_FASTA_genefamily;
					}
					if (con_str_FASTA_genefamily == "IGKJ") {
						con_str_FASTA_genefamily_class = "IGKJ_IGLJ_IGKC_IGLC";
					}
					if (con_str_FASTA_genefamily == "IGLJ") {
						con_str_FASTA_genefamily_class = "IGKJ_IGLJ_IGKC_IGLC";
					}
					if (con_str_FASTA_genefamily == "IGHJ") {
						con_str_FASTA_genefamily_class = "IGHJ_IGHC";
					}
					if (con_str_FASTA_genefamily == "MIGHV") {
						con_str_FASTA_genefamily_class = con_str_FASTA_genefamily;
					}
					if (con_str_FASTA_genefamily == "mA") {
						con_str_FASTA_genefamily_class = "mAB";
					}
					if (con_str_FASTA_genefamily == "CON") {
						con_str_FASTA_genefamily_class = "CONT";
					}
					if (con_str_FASTA_genefamily == "UNIPROT") {
						con_str_FASTA_genefamily_class = "UNIPROT";
					}
					if (ch_parse_FASTA == ' ') {
						con_str_FASTA_genefamily += "_";
					}
					else {
						con_str_FASTA_genefamily += ch_parse_FASTA;
					}
				}
				else {
					con_str_FASTA_genefamily += '*';
				}
			}
			if ((sw_input_FASTA == 0) && (ch_parse_FASTA == '|')) {
				con_str_FASTA_genefamily.clear();
				sw_input_FASTA = 1;
			}
			if ((sw_input_FASTA == 0) && (ch_parse_FASTA != '>')) {
				con_str_parse_FASTA_accession += ch_parse_FASTA;
			}
			if ((sw_2_input_FASTA == 1) && ((par_fin_input_FASTA.peek() == '>') || (par_fin_input_FASTA.peek() == std::ifstream::traits_type::eof()))) {
				if (con_str_FASTA_genefamily_class != "MIGHV") {
					con_c_parse_FASTA.set_str_parse_FASTA_accession(con_str_parse_FASTA_accession);
					con_c_parse_FASTA.set_str_parse_FASTA_genefamily(con_str_FASTA_genefamily);
					con_c_parse_FASTA.set_str_parse_FASTA_genefamily_class(con_str_FASTA_genefamily_class);
					con_c_parse_FASTA.set_str_parse_FASTA_species(con_str_FASTA_species);
					con_c_parse_FASTA.set_str_parse_FASTA_protein(con_str_FASTA_protein);
					con_v_c_parse_FASTA.push_back(con_c_parse_FASTA);
				}
				++st_count_parse_FASTA;
				if (st_count_parse_FASTA % 1000 == 0) {
					std::cout << "\n FASTA accession parse #: ";
					std::cout << st_count_parse_FASTA;
				}
				if ((par_fin_input_FASTA.peek() == std::ifstream::traits_type::eof())) {
					std::cout << "\n FASTA accession parse #: ";
					std::cout << st_count_parse_FASTA;
				}
				con_str_parse_FASTA_accession.clear();
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

	bool check_protein_peptides(vector<fpf_parse::parse_peptides_csv_type> par_v_c_parse_csv_proteinpeptides_data, bool par_filesystem_modified) {
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

	bool check_denovo_peptides(vector<fpf_parse::parse_peptides_csv_type> par_v_c_parse_csv_denovopeptides_data, bool par_filesystem_modified) {
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

	bool b_parse_FASTA_empty(vector<parse_FASTA_type> par_v_c_parse_FASTA) {
		if (par_v_c_parse_FASTA.size() == 0) {
			std::cout << "\n\n * * * FASTA file empty..\n\n * * * Is the file correctly directed?";
			std::cout << "\n\n The program will now terminate. Input any key to continue -\n\n -> ";
			string str_exit;
			std::cin >> str_exit;
			return true;
		}
		return false;
	}

	void output_v_c_parse_FASTA(vector<parse_FASTA_type> par_v_c_parse_FASTA) {
		string output_FASTA_filtered = "FASTA\\output.fasta";
		std::ofstream fout_FASTA_filtered;
		fout_FASTA_filtered.open(output_FASTA_filtered);
		if (IgFamily::OUTPUT_FASTA == 1) {
			for (auto itr_v_c_parse_FASTA : par_v_c_parse_FASTA) {
				fout_FASTA_filtered << ">" << itr_v_c_parse_FASTA.return_str_parse_FASTA_accession();
				fout_FASTA_filtered << "|" << itr_v_c_parse_FASTA.return_str_parse_FASTA_genefamily();
				fout_FASTA_filtered << "|" << itr_v_c_parse_FASTA.return_str_parse_FASTA_species();
				fout_FASTA_filtered << "|";
				for (auto i = 0; i < itr_v_c_parse_FASTA.return_str_parse_FASTA_protein().length(); ++i) {
					if ((i % 60 == 0) && ((i + 1) < itr_v_c_parse_FASTA.return_str_parse_FASTA_protein().length())) {
						fout_FASTA_filtered << "\n";
					}
					fout_FASTA_filtered << itr_v_c_parse_FASTA.return_str_parse_FASTA_protein().at(i);
				}
				fout_FASTA_filtered << "\n";
			}
		}
	}

	void output_v_c_parse_FASTA_to_blastdirectory(vector<parse_FASTA_type> par_v_c_parse_FASTA) {
		string output_FASTA_filtered = "blast_directory\\database.fasta";
		std::ofstream fout_FASTA_filtered;
		fout_FASTA_filtered.open(output_FASTA_filtered);
		if (IgFamily::OUTPUT_FASTA == 1) {
			for (auto itr_v_c_parse_FASTA : par_v_c_parse_FASTA) {
				fout_FASTA_filtered << ">" << itr_v_c_parse_FASTA.return_str_parse_FASTA_accession();
				fout_FASTA_filtered << "|" << itr_v_c_parse_FASTA.return_str_parse_FASTA_genefamily();
				fout_FASTA_filtered << "|" << itr_v_c_parse_FASTA.return_str_parse_FASTA_species();
				fout_FASTA_filtered << "|";
				for (auto i = 0; i < itr_v_c_parse_FASTA.return_str_parse_FASTA_protein().length(); ++i) {
					if ((i % 60 == 0) && ((i + 1) < itr_v_c_parse_FASTA.return_str_parse_FASTA_protein().length())) {
						fout_FASTA_filtered << "\n";
					}
					fout_FASTA_filtered << itr_v_c_parse_FASTA.return_str_parse_FASTA_protein().at(i);
				}
				fout_FASTA_filtered << "\n";
			}
		}
	}
}

#endif