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

	struct parse_peptides_csv_type;
	class parse_FASTA_type;

	typedef size_t size_type;
	typedef std::string string_type;
	typedef std::vector<parse_peptides_csv_type> v_s_parse_csv_data_type;

	struct parse_peptides_csv_type {
	public:
		string_type str_parse_peptides_csv_file;
		string_type str_parse_peptides_csv_peptide;
		string_type str_parse_peptides_csv_spectralcount;
		string_type str_parse_csv_IgP;
		std::vector<double> v_d_denovo_localconfidence;
	};

	class parse_FASTA_type {
	public:
		parse_FASTA_type() {
		};

		~parse_FASTA_type() {
		};

		inline void set_str_parse_FASTA_accession(string_type par_str_parse_FASTA_accession) {
			str_parse_FASTA_accession = par_str_parse_FASTA_accession;
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

	std::vector<parse_peptides_csv_type> parse_proteinpeptides(std::ifstream& par_fin_proteinpeptides_csv, string_type par_str_dir) {

		std::vector<parse_peptides_csv_type> con_v_c_parse_proteinpeptides_csv;
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
						parse_peptides_csv_type con_c_parse_csv = parse_peptides_csv_type();
						con_c_parse_csv.str_parse_peptides_csv_file = par_str_dir;
						con_c_parse_csv.str_parse_peptides_csv_peptide = str_parse_peptide;
						con_c_parse_csv.str_parse_peptides_csv_spectralcount = str_parse_spectralcount;
						con_c_parse_csv.str_parse_csv_IgP = str_parse_IgP;
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

	std::vector<parse_peptides_csv_type> parse_denovopeptides_csv(std::ifstream& par_fin_denovopeptides_csv, string_type par_str_dir) {

		std::vector<parse_peptides_csv_type> con_v_c_parse_peptides_csv;
		string_type str_parse_csv;
		size_type st_count_csv = size_type();
		size_type sw_parse_csv = size_type(1);
		size_type sw_parse_csv_ignore_header = size_type();
		char ch_parse_csv;
		string_type str_parse_peptide;
		string_type con_str_denovo_localconfidence;
		std::vector<double> con_v_d_denovo_localconfidence;

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

	std::vector<parse_FASTA_type> parse_FASTA(std::ifstream& par_fin_input_FASTA) {

		char ch_parse_FASTA;
		size_type sw_input_FASTA = size_type();
		size_type sw_2_input_FASTA = size_type();
		string_type con_str_FASTA_genefamily = string_type();
		string_type con_str_FASTA_genefamily_class = string_type();
		string_type con_str_FASTA_protein = string_type();
		string_type con_str_parse_FASTA_accession = string_type();
		string_type con_str_FASTA_species = string_type();
		parse_FASTA_type con_c_parse_FASTA = parse_FASTA_type();
		std::vector<parse_FASTA_type> con_v_c_parse_FASTA;

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
			if ((sw_input_FASTA == 0) && (ch_parse_FASTA != '>')) {
				con_str_parse_FASTA_accession += ch_parse_FASTA;
			}
			if ((sw_2_input_FASTA == 1) && ((par_fin_input_FASTA.peek() == '>') || (par_fin_input_FASTA.peek() == std::ifstream::traits_type::eof()))) {
				con_c_parse_FASTA.set_str_parse_FASTA_accession(con_str_parse_FASTA_accession);
				con_c_parse_FASTA.set_str_parse_FASTA_genefamily(con_str_FASTA_genefamily);
				con_c_parse_FASTA.set_str_parse_FASTA_genefamily_class(con_str_FASTA_genefamily_class);
				con_c_parse_FASTA.set_str_parse_FASTA_species(con_str_FASTA_species);
				con_c_parse_FASTA.set_str_parse_FASTA_protein(con_str_FASTA_protein);
				con_v_c_parse_FASTA.push_back(con_c_parse_FASTA);
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

	bool check_protein_peptides(std::vector<fpf_parse::parse_peptides_csv_type> par_v_c_parse_csv_proteinpeptides_data, bool par_filesystem_modified) {
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

	bool check_denovo_peptides(std::vector<fpf_parse::parse_peptides_csv_type> par_v_c_parse_csv_denovopeptides_data, bool par_filesystem_modified) {
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

	bool b_parse_FASTA_empty(std::vector<parse_FASTA_type> par_v_c_parse_FASTA) {
		if (par_v_c_parse_FASTA.size() == 0) {
			std::cout << "\n\n * * * FASTA file empty..\n\n * * * Is the file correctly directed?";
			std::cout << "\n\n The program will now terminate. Input any key to continue -\n\n -> ";
			string_type str_exit;
			std::cin >> str_exit;
			return true;
		}
		return false;
	}

	void output_v_c_parse_FASTA(std::vector<parse_FASTA_type> par_v_c_parse_FASTA) {
		string_type output_FASTA_filtered = "FASTA\\output.fasta";
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