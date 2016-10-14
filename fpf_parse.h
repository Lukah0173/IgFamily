// * * fpf_parse.h * *
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_PARSE
#define	FPF_PARSE

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <istream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "IgFamily.h"


namespace fpf_parse {

	using std::string;
	using std::vector;

	struct csv_data;
	class FASTA_data;

	typedef vector<csv_data> v_csv_data;

	struct csv_data {
	public:
		string csv_sourcefile;
		string csv_scan_ID;
		string csv_mz;
		string csv_m;
		string csv_rt;
		string csv_z;
		string csv_spectralcount;
		string csv_IgP;
		string csv_peptide;
		vector<double> v_csv_denovo_localconfidence;
	};

	class FASTA_data {
	public:
		FASTA_data() {
		};

		~FASTA_data() {
		};

		inline void set_FASTA_accession(string par_FASTA_accession) {
			FASTA_accession = par_FASTA_accession;
		};

		inline void set_FASTA_name(string par_FASTA_name) {
			FASTA_name = par_FASTA_name;
		};

		inline void set_FASTA_class(string par_FASTA_class) {
			FASTA_class = par_FASTA_class;
		};

		inline void set_FASTA_type(string par_FASTA_type) {
			FASTA_type = par_FASTA_type;
		};

		inline void set_FASTA_species(string par_FASTA_species) {
			FASTA_species = par_FASTA_species;
		};

		inline void set_FASTA_protein(string par_protein_data) {
			FASTA_protein = par_protein_data;
		};

		inline const string return_FASTA_accession() const {
			return FASTA_accession;
		};

		inline const string return_FASTA_name() const {
			return FASTA_name;
		};

		inline const string return_FASTA_class() const {
			return FASTA_class;
		};

		inline const string return_FASTA_type() const {
			return FASTA_type;
		};

		inline const string return_FASTA_species() const {
			return FASTA_species;
		};

		inline const string return_FASTA_protein() const {
			return FASTA_protein;
		};

	private:
		string FASTA_accession;
		string FASTA_name;
		string FASTA_class;
		string FASTA_type;
		string FASTA_species;
		string FASTA_protein;
	};

	vector<csv_data> parse_csv_PEAKS_database_peptides(const string& par_directory) {
		std::ifstream fin_input_csv(par_directory);
		vector<csv_data> temp_v_csv_data{};
		csv_data temp_csv_data{};
		char csv_read{};
		size_t csv_count_delimit{};
		size_t csv_count_delimit_width{};
		size_t csv_peptide_parse_condition{};
		bool csv_header_parsed{};
		while (fin_input_csv.get(csv_read)) {
			if (csv_header_parsed) {
				if (csv_read != ',') {
					if (csv_count_delimit % csv_count_delimit_width == 3) {
						if (csv_read == '(') {
							csv_peptide_parse_condition = 2;
						}
						if (csv_read == ')') {
							csv_peptide_parse_condition = 1;
						}
						if ((csv_peptide_parse_condition == 1) && ((csv_read == '.') && (temp_csv_data.csv_peptide.length() > 2))) {
							csv_peptide_parse_condition = 0;
						}
						if ((csv_peptide_parse_condition == 1) || (csv_peptide_parse_condition == 2)) {
							temp_csv_data.csv_peptide += csv_read;
						}
						if ((csv_read == '.') && (temp_csv_data.csv_peptide.length() <= 2)) {
							temp_csv_data.csv_peptide.clear();
							csv_peptide_parse_condition = 1;
						}
						if ((csv_peptide_parse_condition != 2) && (csv_read != '.')) {
							csv_peptide_parse_condition = 1;
						}
					}
					if (csv_count_delimit % csv_count_delimit_width == 5) {
						temp_csv_data.csv_IgP += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 6) {
						temp_csv_data.csv_m += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 9) {
						temp_csv_data.csv_mz += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 10) {
						temp_csv_data.csv_z += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 11) {
						temp_csv_data.csv_rt += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 12) {
						temp_csv_data.csv_scan_ID += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 13) {
						temp_csv_data.csv_sourcefile += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 14) {
						temp_csv_data.csv_spectralcount += csv_read;
					}
				}
				if (csv_read == ',') {
					++csv_count_delimit;
				}
				if (csv_read == '\n') {
					for (const auto& itr_csv_peptide : temp_csv_data.csv_peptide) {
						temp_csv_data.v_csv_denovo_localconfidence.push_back(100);
					}
					for (auto i = 0; i < std::stoi(temp_csv_data.csv_spectralcount); ++i) {
						temp_v_csv_data.push_back(temp_csv_data);
					}
					temp_csv_data = csv_data();
				}
			}
			if (!csv_header_parsed) {
				if (csv_read == ',') {
					++csv_count_delimit;
				}
				if (csv_read == '\n') {
					csv_count_delimit_width = csv_count_delimit;
					csv_header_parsed = true;
				}
			}
		}

		fin_input_csv.clear();
		fin_input_csv.seekg(0, std::ios::beg);

		return temp_v_csv_data;
	}

	vector<csv_data> parse_csv_PEAKS_denovo_peptides(string par_directory) {
		std::ifstream fin_input_csv(par_directory);
		vector<csv_data> temp_v_csv_data{};
		csv_data temp_csv_data{};
		char csv_read{};
		size_t csv_count_delimit{};
		size_t csv_count_delimit_width{};
		size_t csv_peptide_parse_condition{};
		bool csv_header_parsed{};
		bool read_source{ true };
		string temp_csv_denovo_localconfidence{};
		while (fin_input_csv.get(csv_read)) {
			if (csv_header_parsed) {
				if (csv_read == ',') {
					if (csv_count_delimit % csv_count_delimit_width == 0) {
						if (!temp_csv_data.csv_scan_ID.empty()) {
							temp_v_csv_data.push_back(temp_csv_data);
						}
						temp_csv_data = {};
					}
					++csv_count_delimit;
				}
				if (csv_read != ',') {
					if (csv_count_delimit % csv_count_delimit_width == 1) {
						temp_csv_data.csv_scan_ID += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 2) {
						if (read_source) {
							if (csv_read == '.') {
								read_source = false;
							}
							if (!IgFamily::FILESYSTEM_MODE) {
								temp_csv_data.csv_sourcefile += csv_read;
							}
						}
					}
					if ((csv_count_delimit % csv_count_delimit_width == 3) && (csv_read != ' ')) {
						if (csv_read == '(') {
							csv_peptide_parse_condition = 2;
						}
						if (csv_read == ')') {
							csv_peptide_parse_condition = 1;
						}
						if ((csv_peptide_parse_condition == 1) && ((csv_read == '.') && (temp_csv_data.csv_peptide.length() > 2))) {
							csv_peptide_parse_condition = 0;
						}
						if ((csv_peptide_parse_condition == 1) || (csv_peptide_parse_condition == 2)) {
							temp_csv_data.csv_peptide += csv_read;
						}
						if ((csv_read == '.') && (temp_csv_data.csv_peptide.length() <= 2)) {
							temp_csv_data.csv_peptide.clear();
							csv_peptide_parse_condition = 1;
						}
						if ((csv_peptide_parse_condition != 2) && (csv_read != '.')) {
							csv_peptide_parse_condition = 1;
						}
					}
					if (csv_count_delimit % csv_count_delimit_width == 7) {
						temp_csv_data.csv_mz += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 8) {
						temp_csv_data.csv_z += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 9) {
						temp_csv_data.csv_rt += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 11) {
						temp_csv_data.csv_m += csv_read;
					}
					if (csv_count_delimit % csv_count_delimit_width == 14) {
						if (csv_read == ' ') {
							temp_csv_data.v_csv_denovo_localconfidence.push_back(std::stod(temp_csv_denovo_localconfidence));
							temp_csv_denovo_localconfidence.clear();
						}
						if (csv_read != ' ') {
							temp_csv_denovo_localconfidence += csv_read;
						}
						if (fin_input_csv.peek() == ',') {
							temp_csv_data.v_csv_denovo_localconfidence.push_back(std::stod(temp_csv_denovo_localconfidence));
							temp_csv_denovo_localconfidence.clear();
						}
					}
				}
			}
			if (!csv_header_parsed) {
				if (csv_read == ',') {
					++csv_count_delimit;
				}
				if (csv_read == '\n') {
					csv_count_delimit_width = csv_count_delimit;
					csv_header_parsed = true;
				}
			}
		}

		fin_input_csv.clear();
		fin_input_csv.seekg(0, std::ios::beg);

		return temp_v_csv_data;
	}

	vector<csv_data> parse_csv_NOVOR_denovo_peptides(string par_directory) {
		std::ifstream fin_input_csv(par_directory);
		vector<csv_data> temp_v_csv_data{};
		csv_data temp_csv_data{};
		char csv_read{};
		size_t csv_count_delimit{};
		size_t csv_count_delimit_width{};
		size_t csv_peptide_parse_condition{};
		bool csv_header_parsed{};
		bool read_source{ true };
		string temp_csv_denovo_localconfidence{};
		while (fin_input_csv.get(csv_read)) {
			if (csv_header_parsed) {
				if (csv_read == ',') {
					temp_v_csv_data.push_back(temp_csv_data);
					temp_csv_data = csv_data();
					++csv_count_delimit;
				}
				if (csv_read != ',') {
					if (csv_count_delimit % csv_count_delimit_width == 10) {
						if (csv_read == '-') {
							temp_csv_data.v_csv_denovo_localconfidence.push_back(std::stod(temp_csv_denovo_localconfidence));
							temp_csv_denovo_localconfidence.clear();
						}
						if ((csv_read != ' ') && (csv_read != '-')) {
							temp_csv_denovo_localconfidence += csv_read;
						}
						if (fin_input_csv.peek() == ',') {
							temp_csv_data.v_csv_denovo_localconfidence.push_back(std::stod(temp_csv_denovo_localconfidence));
							temp_csv_denovo_localconfidence.clear();
						}
					}
					if ((csv_count_delimit % csv_count_delimit_width == 9) && (csv_read != ' ')) {
						if (csv_read == '(') {
							csv_peptide_parse_condition = 2;
						}
						if (csv_read == ')') {
							csv_peptide_parse_condition = 1;
						}
						if ((csv_peptide_parse_condition == 1) && ((csv_read == '.') && (temp_csv_data.csv_peptide.length() > 2))) {
							csv_peptide_parse_condition = 0;
						}
						if ((csv_peptide_parse_condition == 1) || (csv_peptide_parse_condition == 2)) {
							temp_csv_data.csv_peptide += csv_read;
						}
						if ((csv_read == '.') && (temp_csv_data.csv_peptide.length() <= 2)) {
							temp_csv_data.csv_peptide.clear();
							csv_peptide_parse_condition = 1;
						}
						if ((csv_peptide_parse_condition != 2) && (csv_read != '.')) {
							csv_peptide_parse_condition = 1;
						}
					}
				}
			}
			if (!csv_header_parsed) {
				if (csv_read == ',') {
					++csv_count_delimit;
				}
				if (csv_read == '\n') {
					csv_count_delimit_width = csv_count_delimit;
					csv_header_parsed = true;
				}
			}
		}

		fin_input_csv.clear();
		fin_input_csv.seekg(0, std::ios::beg);

		return temp_v_csv_data;
	}

	void check_FASTA_format(string par_INPUT_FASTA) {

		std::ifstream fin_INPUT_FASTA(par_INPUT_FASTA);
		char FASTA_read{};
		string FASTA_parse{};
		size_t FASTA_count_delimit{};
		bool read_format{};
		bool header_line{};

		while (fin_INPUT_FASTA.get(FASTA_read)) {
			if (FASTA_read == '>') {
				header_line = true;
				if (read_format) {
					std::cout << "\n\n Accession - \n";
					std::cout << FASTA_parse;
					FASTA_parse.clear();
					string str_format_break{};
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
				FASTA_count_delimit = 0;
				read_format = true;
			}
			if (header_line) {
				if ((FASTA_read == '|') || (FASTA_read == '\n') || (FASTA_read == '\r\n') || (FASTA_read == '\r')) {
					std::cout << "\n Field # ";
					std::cout << FASTA_count_delimit;
					std::cout << " - \"";
					std::cout << FASTA_parse;
					std::cout << "\"";
					++FASTA_count_delimit;
					FASTA_parse.clear();
				}
				if ((FASTA_read != '>') && (FASTA_read != '|') && ((FASTA_read != '\n') && (FASTA_read != '\r\n') && (FASTA_read != '\r'))) {
					FASTA_parse += FASTA_read;
				}
				if ((FASTA_read == '\n') || (FASTA_read == '\r\n') || (FASTA_read == '\r')) {
					header_line = false;
					FASTA_parse.clear();
				}
			}
			if (!header_line) {
				FASTA_parse += FASTA_read;
			}
		}
		fin_INPUT_FASTA.clear();
		fin_INPUT_FASTA.seekg(0, std::ios::beg);
	}

	void output_custom_FASTA_format(string par_INPUT_FASTA) {

		std::ifstream fin_INPUT_FASTA(par_INPUT_FASTA);
		char FASTA_read{};
		string FASTA_parse{};
		size_t FASTA_count_delimit{};
		bool read_format{};
		bool header_line{};

		while (fin_INPUT_FASTA.get(FASTA_read)) {
			if (FASTA_read == '>') {
				header_line = true;
				if (read_format) {
					std::cout << "\n\n Accession - \n";
					std::cout << FASTA_parse;
					FASTA_parse.clear();
					std::cout << "\n\n current FASTA output format: ";
					std::cout << "\n\n >";
					for (auto i = 0; i < FASTA_count_delimit; ++i) {
						std::cout << "[Field #";
						std::cout << i;
						std::cout << "]|";
					}
					std::cout << "\n [Accession]";
					std::cout << "\n\n\n Create custom output format:           [C] ";
					std::cout << "\n Output FASTA with selected format:     [X] ";
					string menu_selection{};
					while ((menu_selection != ("C")) && (menu_selection != ("X"))) {
						menu_selection.clear();
						std::cout << "\n\n Input selection: \n\n --> ";
						std::cin >> menu_selection;
					}
					if (menu_selection == "C") {
						while (true) {
							vector<string> v_field_value(FASTA_count_delimit);
							for (auto i = 0; i < FASTA_count_delimit; ++i) {
								v_field_value[i] = std::to_string(i);
								std::cout << "\n ";
								std::cout << i;
								std::cout << ":     Modify Field #";
								std::cout << i;
							}
							std::cout << "\n T:     Truncate field output";
							std::cout << "\n X:     Continue";
							string change_field_value = "-1";
							while (!((std::stoi(change_field_value) < FASTA_count_delimit) && (std::stoi(change_field_value) >= 0))) {
								std::cout << "\n\n Input selection: \n\n --> ";
								std::cin >> change_field_value;
								if ((change_field_value != "T") && (change_field_value != "X")) {
									try {
										auto try_catch = std::stoi(change_field_value);
										if (try_catch == 0) {
										}
									}
									catch (const std::invalid_argument&) {
										change_field_value = "-1";
									}
								}
								else {
									break;
								}
							}
							if (change_field_value == "T") {
								size_t st_truncate_v_str_field_value = size_t();
								std::cout << "\n\n Truncate to how many fields?: \n\n --> ";
								while (!((st_truncate_v_str_field_value > 0) && (st_truncate_v_str_field_value < FASTA_count_delimit))) {
									std::cin >> st_truncate_v_str_field_value;
								}
								v_field_value.resize(st_truncate_v_str_field_value);
								std::cout << "\n\n current FASTA output format: ";
								std::cout << "\n\n >";
								for (auto i = 0; i < v_field_value.size(); ++i) {
									std::cout << "[Field #";
									std::cout << v_field_value[i];
									std::cout << "]|";
								}
								std::cout << "\n [Accession]";
								std::cout << "\n ";
							}
							if (change_field_value == "X") {
								break;
							}
							if ((change_field_value != "T") && (change_field_value != "X")) {
								std::cout << "\n\n Modifying Field #";
								std::cout << change_field_value;
								std::cout << "\n";
								for (auto i = 0; i < FASTA_count_delimit; ++i) {
									v_field_value[i] = std::to_string(i);
									std::cout << "\n ";
									std::cout << i;
									std::cout << ":     Assign Field #";
									std::cout << i;
								}
								std::cout << "\n ";
								std::cout << FASTA_count_delimit;
								std::cout << ":     Custom string entry";
								std::cout << "\n ";
								size_t assign_field_value = -1;
								while (!((assign_field_value <= FASTA_count_delimit) && (assign_field_value >= 0))) {
									std::cout << "\n\n Input selection: \n\n --> ";
									std::cin >> assign_field_value;
									if (std::cin.fail()) {
										std::cin.clear();
										std::cin.ignore(256, '\n');
										assign_field_value = -1;
									}
								}
								if (assign_field_value == FASTA_count_delimit) {
									string assign_field_value;
									std::cout << "\n\n Input selection: \n\n --> ";
									std::cin >> assign_field_value;
									v_field_value[std::stoi(change_field_value)] = assign_field_value;
								}
								if (assign_field_value < FASTA_count_delimit) {
									v_field_value[std::stoi(change_field_value)] = std::to_string(assign_field_value);
								}
								std::cout << "\n\n current FASTA output format: ";
								std::cout << "\n\n >";
								for (auto i = 0; i < v_field_value.size(); ++i) {
									std::cout << "[Field #";
									std::cout << v_field_value[i];
									std::cout << "]|";
								}
								std::cout << "\n [Accession]";
								std::cout << "\n ";
							}
						}
					}
					if (menu_selection == "X") {
					}
					break;
				}
				FASTA_count_delimit = 0;
				read_format = true;
			}
			if (header_line) {
				if ((FASTA_read == '|') || (FASTA_read == '\n') || (FASTA_read == '\r\n') || (FASTA_read == '\r')) {
					std::cout << "\n Field # ";
					std::cout << FASTA_count_delimit;
					std::cout << " - \"";
					std::cout << FASTA_parse;
					std::cout << "\"";
					++FASTA_count_delimit;
					FASTA_parse.clear();
				}
				if ((FASTA_read != '>') && (FASTA_read != '|') && ((FASTA_read != '\n') && (FASTA_read != '\r\n') && (FASTA_read != '\r'))) {
					FASTA_parse += FASTA_read;
				}
				if ((FASTA_read == '\n') || (FASTA_read == '\r\n') || (FASTA_read == '\r')) {
					header_line = false;
					FASTA_parse.clear();
				}
			}
			if (!header_line) {
				FASTA_parse += FASTA_read;
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
		char FASTA_get{};
		string FASTA_parse{};
		size_t FASTA_count_delimit{};
		bool read_format{};
		bool header_line{};

		while (fin_INPUT_FASTA.get(FASTA_get)) {
			if (FASTA_get == '>') {
				header_line = true;
				if (read_format) {
					fout_custom_FASTA << "Homo sapiens|";
					fout_custom_FASTA << FASTA_parse;
					FASTA_parse.clear();
				}
				FASTA_count_delimit = 0;
				read_format = true;
			}
			if (header_line) {
				if ((FASTA_get == '|') || (FASTA_get == '\n') || (FASTA_get == '\r\n') || (FASTA_get == '\r')) {
					if (FASTA_count_delimit == 0) {
						fout_custom_FASTA << ">";
					}
					if (FASTA_count_delimit == 1) {
						fout_custom_FASTA << FASTA_parse;
						fout_custom_FASTA << "|";
					}
					if (FASTA_count_delimit == 2) {
						fout_custom_FASTA << "UNIPROT(";
						fout_custom_FASTA << FASTA_parse;
						fout_custom_FASTA << ")";
						fout_custom_FASTA << "|";
					}
					if (FASTA_count_delimit == 3) {
						fout_custom_FASTA << "Homo sapiens";
						fout_custom_FASTA << "|";
					}
					++FASTA_count_delimit;
					FASTA_parse.clear();
				}
				if ((FASTA_get != '>') && (FASTA_get != '|') && ((FASTA_get != '\n') && (FASTA_get != '\r\n') && (FASTA_get != '\r'))) {
					FASTA_parse += FASTA_get;
				}
				if ((FASTA_get == '\n') || (FASTA_get == '\r\n') || (FASTA_get == '\r')) {
					header_line = false;
					FASTA_parse.clear();
				}
			}
			if (!header_line) {
				FASTA_parse += FASTA_get;
			}
		}
		fin_INPUT_FASTA.clear();
		fin_INPUT_FASTA.seekg(0, std::ios::beg);
	}

	vector<FASTA_data> parse_FASTA(const string& par_select_FASTA) {
		vector<FASTA_data> temp_v_FASTA_data{};
		FASTA_data temp_FASTA_data{};
		std::ifstream fin_FASTA(par_select_FASTA);
		std::cout << par_select_FASTA << "\n";
		char read_FASTA{};
		string parse_FASTA{};
		size_t count_FASTA_delimit{};
		size_t count_FASTA_delimit_width{ 4 };
		size_t count_FASTA_accessions{};
		bool reading_accession{};
		while (fin_FASTA.get(read_FASTA)) {
			if (read_FASTA == '>') {
				reading_accession = true;
				++count_FASTA_accessions;
			}
			if (reading_accession) {
				if (read_FASTA == '\n') {
					reading_accession = false;
				}
				if (read_FASTA == '|') {
					if (count_FASTA_delimit % count_FASTA_delimit_width == 0) {
						temp_FASTA_data.set_FASTA_accession(parse_FASTA);
						parse_FASTA.clear();
					}
					if (count_FASTA_delimit % count_FASTA_delimit_width == 1) {
						temp_FASTA_data.set_FASTA_name(parse_FASTA);
						parse_FASTA.clear();
					}
					if (count_FASTA_delimit % count_FASTA_delimit_width == 2) {
						if (parse_FASTA == "IMGT_V") {
							temp_FASTA_data.set_FASTA_type("IGV");
						}
						if (parse_FASTA == "CONT") {
							temp_FASTA_data.set_FASTA_type("CONT");
						}
						if (parse_FASTA == "UNIPROT") {
							temp_FASTA_data.set_FASTA_type("UNIPROT");
						}
						if (parse_FASTA == "mAB") {
							temp_FASTA_data.set_FASTA_type("IGV");
						}
						parse_FASTA.clear();
					}
					if (count_FASTA_delimit % count_FASTA_delimit_width == 3) {
						temp_FASTA_data.set_FASTA_species(parse_FASTA);
						parse_FASTA.clear();
					}
					++count_FASTA_delimit;
				}
				if ((read_FASTA != '|') && (read_FASTA != '>') && (read_FASTA != '\n')) {
					parse_FASTA += read_FASTA;
				}
			}
			if (!reading_accession) {
				if (read_FASTA != '\n') {
					parse_FASTA += read_FASTA;
				}
				if ((fin_FASTA.peek() == '>') || (fin_FASTA.peek() == std::ifstream::traits_type::eof())) {
					temp_FASTA_data.set_FASTA_protein(parse_FASTA);
					parse_FASTA.clear();
					if (!((temp_FASTA_data.return_FASTA_type() == "UNIPROT")
						&& ((temp_FASTA_data.return_FASTA_name().find("Ig") != std::string::npos)
							|| (temp_FASTA_data.return_FASTA_name().find("Immuno") != std::string::npos)
							|| (temp_FASTA_data.return_FASTA_name().find("immuno") != std::string::npos)
							|| (temp_FASTA_data.return_FASTA_name().find("HLA") != std::string::npos)))) {
						temp_v_FASTA_data.push_back(temp_FASTA_data);
					}
				}
			}
		}

		if (fin_FASTA.peek() == std::ifstream::traits_type::eof()) {
			std::cout << "\n FASTA accessions parsed - ";
			std::cout << count_FASTA_accessions;
		}	
		
		fin_FASTA.clear();
		fin_FASTA.seekg(0, std::ios::beg);

		return temp_v_FASTA_data;
	}

	inline bool check_csv_PEAKS_database_peptides_empty(vector<fpf_parse::csv_data> par_v_csv_data) {
		if (!par_v_csv_data.empty()) {
			std::cout << "\n --- database matched peptides file found";
			return true;
		}
		return false;
	}

	inline bool check_csv_PEAKS_denovo_peptides_empty(vector<fpf_parse::csv_data> par_v_csv_data) {
		if (!par_v_csv_data.empty()) {
			std::cout << "\n --- PEAKS de novo peptides file found";
			return true;
		}
		return false;
	}

	inline bool check_csv_NOVOR_denovo_peptides_empty(vector<fpf_parse::csv_data> par_v_csv_data) {
		if (!par_v_csv_data.empty()) {
			std::cout << "\n --- NOVOR de novo peptides file found";
			return true;
		}
		return false;
	}

	inline bool check_FASTA_file_exists(vector<FASTA_data> par_v_FASTA_data) {
		if (par_v_FASTA_data.empty()) {
			std::cout << "\n\n * * * FASTA file empty..\n\n * * * Is the file correctly directed?";
			std::cout << "\n\n The program will now terminate. Input any key to continue -\n\n -> ";
			string program_exit;
			std::cin >> program_exit;
			return true;
		}
		return false;
	}

	void output_v_FASTA_data(vector<FASTA_data> par_v_FASTA_data) {
		string output_FASTA_filtered = "FASTA\\output.fasta";
		std::ofstream fout_FASTA_filtered;
		fout_FASTA_filtered.open(output_FASTA_filtered);
		if (IgFamily::OUTPUT_FASTA == 1) {
			for (const auto& itr_v_FASTA_data : par_v_FASTA_data) {
				fout_FASTA_filtered << ">" << itr_v_FASTA_data.return_FASTA_accession();
				fout_FASTA_filtered << "|" << itr_v_FASTA_data.return_FASTA_name();
				fout_FASTA_filtered << "|" << itr_v_FASTA_data.return_FASTA_species();
				fout_FASTA_filtered << "|";
				for (auto i = 0; i < itr_v_FASTA_data.return_FASTA_protein().length(); ++i) {
					if ((i % 60 == 0) && ((i + 1) < itr_v_FASTA_data.return_FASTA_protein().length())) {
						fout_FASTA_filtered << "\n";
					}
					fout_FASTA_filtered << itr_v_FASTA_data.return_FASTA_protein().at(i);
				}
				fout_FASTA_filtered << "\n";
			}
		}
	}
}

#endif