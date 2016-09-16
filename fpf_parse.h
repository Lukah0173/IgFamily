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
		string csv_scan_ID;
		string csv_sourcefile;
		string csv_peptide;
		string csv_spectralcount;
		string csv_IgP;
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

		inline void set_protein_data(string par_protein_data) {
			protein_data = par_protein_data;
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

		inline const string return_protein_data() const {
			return protein_data;
		};

	private:
		string FASTA_accession;
		string protein_data_delimited;
		string FASTA_name;
		string FASTA_class;
		string FASTA_type;
		string FASTA_species;
		string protein_data;
	};

	vector<csv_data> parse_proteinpeptides(const string& par_directory) {
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
				if (csv_read == ',') {
					if (csv_count_delimit % csv_count_delimit_width == 0) {
						temp_v_csv_data.push_back(temp_csv_data);
						temp_csv_data = csv_data();
					}
					++csv_count_delimit;
				}
				if (csv_read != ',') {
					if (csv_count_delimit % csv_count_delimit_width == 1) {
						temp_csv_data.csv_scan_ID += csv_read;
					}
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
					if ((csv_count_delimit % csv_count_delimit_width == 5) && (csv_read != ',')) {
						temp_csv_data.csv_IgP += csv_read;
					}
					if ((csv_count_delimit % csv_count_delimit_width == 14) && (csv_read != ',')) {
						temp_csv_data.csv_spectralcount += csv_read;
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

	vector<csv_data> parse_csv_PEAKS_denovopeptides(string par_directory) {
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
						temp_v_csv_data.push_back(temp_csv_data);
						temp_csv_data = csv_data();				
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

	vector<csv_data> parse_csv_NOVOR_denovopeptides(string par_directory) {
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

	vector<FASTA_data> parse_FASTA(std::ifstream& par_fin_input_FASTA) {
		vector<FASTA_data> temp_v_FASTA_data{};
		FASTA_data temp_FASTA_data{};
		char read_FASTA{};
		size_t FASTA_condition_switch{};
		size_t FASTA_condition_switch_2{};
		size_t count_FASTA_delimit{};
		size_t count_FASTA_delimit_width{ 4 };
		bool reading_accession{};
		while (par_fin_input_FASTA.get(read_FASTA)) {
			if (read_FASTA == '>') {
				reading_accession = true;
			}
			if (reading_accession) {
				if (read_FASTA == '\n') {
					reading_accession = false;
				}
				if (read_FASTA == '|') {
					if (read_FASTA)
					++count_FASTA_delimit;
				}
			}




			if (FASTA_condition_switch_2 == 1) {
				if (read_FASTA != '\n') {
					temp_FASTA_element += read_FASTA;
				}
			}
			if ((FASTA_condition_switch_2 == 0) && (read_FASTA == '\n')) {
				temp_FASTA_element.clear();
				FASTA_condition_switch_2 = 1;
			}
			if ((FASTA_condition_switch == 2) && (read_FASTA == '|')) {
				FASTA_condition_switch = 3;
			}
			if (FASTA_condition_switch == 2) {
				tenp_FASTA_species += read_FASTA;
			}
			if ((FASTA_condition_switch == 1) && (read_FASTA == '|')) {
				FASTA_condition_switch = 2;
			}
			if (FASTA_condition_switch == 1) {
				if (read_FASTA != '_') {
					if (temp_FASTA_name == "IGHV") {
						temp_FASTA_class = "IGHV";
						temp_FASTA_type = "IG";
					}
					if (temp_FASTA_name == "IGKV") {
						temp_FASTA_class = "IGKV";
						temp_FASTA_type = "IG";
					}
					if (temp_FASTA_name == "IGLV") {
						temp_FASTA_class = "IGLV";
						temp_FASTA_type = "IG";
					}
					if (temp_FASTA_name == "IGKJ") {
						temp_FASTA_class = "IGKJ_IGLJ_IGKC_IGLC";
					}
					if (temp_FASTA_name == "IGLJ") {
						temp_FASTA_class = "IGKJ_IGLJ_IGKC_IGLC";
					}
					if (temp_FASTA_name == "IGHJ") {
						temp_FASTA_class = "IGHJ_IGHC";
					}
					if (temp_FASTA_name == "MIGHV") {
						temp_FASTA_class = "MIGHV";
					}
					if (temp_FASTA_name == "mA") {
						temp_FASTA_class = "mAB";
						temp_FASTA_type = "IG";
					}
					if (temp_FASTA_name == "CON") {
						temp_FASTA_class = "CONT";
					}
					if (temp_FASTA_name == "UNIPROT") {
						temp_FASTA_class = "UNIPROT";
					}
					if (read_FASTA == ' ') {
						temp_FASTA_name += "_";
					}
					else {
						temp_FASTA_name += read_FASTA;
					}
				}
				else {
					temp_FASTA_name += '*';
				}
			}
			if ((FASTA_condition_switch == 0) && (read_FASTA == '|')) {
				temp_FASTA_name.clear();
				FASTA_condition_switch = 1;
			}
			if ((FASTA_condition_switch == 0) && (read_FASTA != '>')) {
				temp_FASTA_accession += read_FASTA;
			}
			if ((FASTA_condition_switch_2 == 1) && ((par_fin_input_FASTA.peek() == '>') || (par_fin_input_FASTA.peek() == std::ifstream::traits_type::eof()))) {

				if (temp_FASTA_class != "MIGHV") {
					if (!((temp_FASTA_class == "UNIPROT") && ((temp_FASTA_name.find("Ig") != std::string::npos) || (temp_FASTA_name.find("Immunoglobulin") != std::string::npos)))) {
						temp_FASTA_data.set_FASTA_accession(temp_FASTA_accession);
						temp_FASTA_data.set_FASTA_name(temp_FASTA_name);
						temp_FASTA_data.set_FASTA_class(temp_FASTA_class);
						temp_FASTA_data.set_FASTA_type(temp_FASTA_type);
						temp_FASTA_data.set_FASTA_species(tenp_FASTA_species);
						temp_FASTA_data.set_protein_data(temp_FASTA_element);
						temp_v_FASTA_data.push_back(temp_FASTA_data);
					}
				}
				++count_FASTA_delimit;
				if ((par_fin_input_FASTA.peek() == std::ifstream::traits_type::eof())) {
					std::cout << "\n FASTA accessions parsed - ";
					std::cout << count_FASTA_delimit;
				}
				FASTA_condition_switch = 0;
				FASTA_condition_switch_2 = 0;
				temp_FASTA_accession.clear();
				temp_FASTA_element.clear();
				temp_FASTA_name.clear();
				temp_FASTA_class.clear();
				temp_FASTA_type.clear();
				tenp_FASTA_species.clear();
			}
		}

		par_fin_input_FASTA.clear();
		par_fin_input_FASTA.seekg(0, std::ios::beg);

		return temp_v_FASTA_data;
	}

	bool check_csv_PEAKS_database_peptides_empty(vector<fpf_parse::csv_data> par_v_csv_data) {
		if (!par_v_csv_data.empty()) {
			std::cout << "\n --- database matched peptides file found";
			return true;
		}
		return false;
	}

	bool check_csv_PEAKS_denovo_peptides_empty(vector<fpf_parse::csv_data> par_v_csv_data) {
		if (!par_v_csv_data.empty()) {
			std::cout << "\n --- PEAKS de novo peptides file found";
			return true;
		}
		return false;
	}

	bool check_csv_NOVOR_denovo_peptides_empty(vector<fpf_parse::csv_data> par_v_csv_data) {
		if (!par_v_csv_data.empty()) {
			std::cout << "\n --- NOVOR de novo peptides file found";
			return true;
		}
		return false;
	}

	bool check_FASTA_file_exists(vector<FASTA_data> par_v_FASTA_data) {
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
				for (auto i = 0; i < itr_v_FASTA_data.return_protein_data().length(); ++i) {
					if ((i % 60 == 0) && ((i + 1) < itr_v_FASTA_data.return_protein_data().length())) {
						fout_FASTA_filtered << "\n";
					}
					fout_FASTA_filtered << itr_v_FASTA_data.return_protein_data().at(i);
				}
				fout_FASTA_filtered << "\n";
			}
		}
	}

	void output_v_FASTA_data_to_blastdirectory(vector<FASTA_data> par_v_FASTA_data) {
		string output_FASTA_filtered = "blast_directory\\database.fasta";
		std::ofstream fout_FASTA_filtered;
		fout_FASTA_filtered.open(output_FASTA_filtered);
		if (IgFamily::OUTPUT_FASTA) {
			for (const auto& itr_v_FASTA_DATA : par_v_FASTA_data) {
				fout_FASTA_filtered << ">" << itr_v_FASTA_DATA.return_FASTA_accession();
				fout_FASTA_filtered << "|" << itr_v_FASTA_DATA.return_FASTA_name();
				fout_FASTA_filtered << "|" << itr_v_FASTA_DATA.return_FASTA_species();
				fout_FASTA_filtered << "|";
				for (auto i = 0; i < itr_v_FASTA_DATA.return_protein_data().length(); ++i) {
					if ((i % IgFamily::OUTPUT_FASTA_ACCESSION_WIDTH == 0) && ((i + 1) < itr_v_FASTA_DATA.return_protein_data().length())) {
						fout_FASTA_filtered << "\n";
					}
					fout_FASTA_filtered << itr_v_FASTA_DATA.return_protein_data().at(i);
				}
				fout_FASTA_filtered << "\n";
			}
		}
	}
}

#endif