// * * fpf_blastp_analysis.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_BLASTP_ANALYSIS 
#define	FPF_BLASTP_ANALYSIS
#include <string> // provides - std::string
#include <vector> // provides - std::vector
#include <fstream> // provides - std::ifstream
#include <iostream> // provides - std::ofstream
#include <utility> // provides - std::pair
#include <math.h> // provides - log10
#include <algorithm> // provides - std::find_if
#include "fpf_data.h"
#include "fpf_filesystem.h"
#include "fpf_filesystem_analysis.h"

namespace fpf_blastp_analysis {

	typedef std::string string_type;
	typedef fpf_data::multinomial_category_data_type multinomial_category_data_type;
	typedef fpf_data::blastp_type blastp_type;
	typedef std::vector<fpf_data::blastp_type> v_s_blastp_tpye;
	typedef fpf_data::multinomial_type s_multinomial_type;
	typedef fpf_filesystem::filesystem_type filesystem_type;
	typedef std::vector<fpf_filesystem::filesystem_type> v_filesystem_type;

	void create_blastp_input(filesystem_type par_s_filesystem) {
		std::string output_blastp_FASTA = "blast_directory\\";
		output_blastp_FASTA += par_s_filesystem.str_filename;
		output_blastp_FASTA += "_blastp_input.fasta";
		std::ofstream fout_blastp_input_FASTA;
		fout_blastp_input_FASTA.open(output_blastp_FASTA);
		size_type st_count_blastp_FASTA = size_type();
		for (auto itr_v_s_peptide_data : par_s_filesystem.v_s_peptide_data) {
			++st_count_blastp_FASTA;
			fout_blastp_input_FASTA << ">" << itr_v_s_peptide_data.str_peptide_filtered << "\n";
			fout_blastp_input_FASTA << itr_v_s_peptide_data.str_peptide_filtered << "\n";
		}
		fout_blastp_input_FASTA.close();
	}

	void create_blastp_database(filesystem_type par_s_filesystem) {
		string_type output_blastp_database_FASTA = "blast_directory\\blastp_database.fasta";
		std::ofstream fout_blastp_database_FASTA;
		fout_blastp_database_FASTA.open(output_blastp_database_FASTA);
		for (auto itr_v_c_multinomial_catagory_distinct : par_s_filesystem.v_c_multinomial_catagory_distinct) {
			fout_blastp_database_FASTA << ">" << itr_v_c_multinomial_catagory_distinct.str_multinomial_category_name << "\n";
			for (size_type i = 0; i < itr_v_c_multinomial_catagory_distinct.str_protein.length(); ++i) {
				if ((i % 60 == 0) && (i != 0)) {
					fout_blastp_database_FASTA << "\n";
				}
				fout_blastp_database_FASTA << itr_v_c_multinomial_catagory_distinct.str_protein.at(i);
			}
			fout_blastp_database_FASTA << "\n";
		}
		fout_blastp_database_FASTA.close();
	}

	void sys_blastp(filesystem_type par_s_filesystem) {
		std::cout << "\n\n";
		string_type string_system = "CD Z:\\Lukah_Dykes\\IgFamily\\blast_directory\\";
		string_system += " && makeblastdb.exe -in ";
		string_system += "database.fasta";
		string_system += " -dbtype prot -out FPF_blastpdb";
		system(string_system.c_str());
		string_system = "CD Z:\\Lukah_Dykes\\IgFamily\\blast_directory\\";
		string_system += " && blastp.exe -query ";
		string_system += par_s_filesystem.str_filename;
		string_system += "_blastp_input.fasta -db FPF_blastpdb -evalue ";
		string_system += std::to_string(EVALUE_THRESHOLD);
		string_system += " -max_target_seqs 200 -out ";
		string_system += par_s_filesystem.str_filename;
		string_system += "_blastp_output.csv -outfmt \"10 qacc qseq sseq sacc qstart sstart qlen pident ppos mismatch evalue\"";
		system(string_system.c_str());
		string_system = "EXIT";
		system(string_system.c_str());
	}

	void create_v_s_blastp(filesystem_type& par_s_filesystem) {
		if (par_s_filesystem.b_denovopeptides_exist) {
			size_type st_csv_blastp = size_type(1);
			string_type str_blastpoutput = "blast_directory\\" + par_s_filesystem.str_filename + "_blastp_output.csv";
			std::ifstream fin_input_blastp(str_blastpoutput);
			char ch_parse_blastp = char();
			blastp_type con_s_blastp = blastp_type();
			string_type con_str_parse_blastp = string_type();
			while (fin_input_blastp.get(ch_parse_blastp)) {
				if ((ch_parse_blastp != ',') && (ch_parse_blastp != '\n')) {
					con_str_parse_blastp += ch_parse_blastp;
				}
				if (ch_parse_blastp == ',') {
					if (st_csv_blastp % 10 == 1) {
						con_s_blastp.str_blastp_query = con_str_parse_blastp;
					}
					if (st_csv_blastp % 10 == 2) {
						con_s_blastp.str_blastp_query_aligned = con_str_parse_blastp;
					}
					if (st_csv_blastp % 10 == 3) {
						con_s_blastp.str_blastp_subject = con_str_parse_blastp;
					}
					if (st_csv_blastp % 10 == 4) {
						con_s_blastp.str_blastp_subject_accession = con_str_parse_blastp;
					}
					if (st_csv_blastp % 10 == 5) {
						con_s_blastp.st_blastp_query_alignment_index = std::stoi(con_str_parse_blastp);
					}
					if (st_csv_blastp % 10 == 6) {
						con_s_blastp.st_blastp_subject_alignment_index = std::stoi(con_str_parse_blastp);
					}
					con_str_parse_blastp.clear();
					++st_csv_blastp;
				}
				if (ch_parse_blastp == '\n') {
					con_s_blastp.d_blastp_evalue = std::stod(con_str_parse_blastp);
					par_s_filesystem.v_s_blastp.push_back(con_s_blastp);
					con_str_parse_blastp.clear();
				}
			}
		}
	}

	void create_str_protein(filesystem_type& par_s_filesystem) {
		for (auto& itr_v_s_blastp : par_s_filesystem.v_s_blastp) {
			auto find_str_genefamily = std::find_if(par_s_filesystem.v_c_multinomial_catagory.begin(), par_s_filesystem.v_c_multinomial_catagory.end(),
				[itr_v_s_blastp](multinomial_category_data_type par_s_multinomial_element_data) {
				return par_s_multinomial_element_data.str_multinomial_category_name == itr_v_s_blastp.str_blastp_subject_accession;
			});
			if (find_str_genefamily == par_s_filesystem.v_c_multinomial_catagory.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_v_s_blastp.str_blastp_subject_accession;
				string_type catch_error;
				std::cin >> catch_error;
			}
			itr_v_s_blastp.str_protein = find_str_genefamily->str_protein;
			itr_v_s_blastp.str_blastp_subject_accession_class = find_str_genefamily->str_multinomial_category_class;

		}
	}

	void create_str_query_alignment(filesystem_type& par_s_filesystem) {
		string_type con_str_query_alignment = string_type();
		size_type st_index_match = 1;
		for (auto& itr_v_s_blastp : par_s_filesystem.v_s_blastp) {
			for (auto itr_str_protein : itr_v_s_blastp.str_protein) {
				if (st_index_match == itr_v_s_blastp.st_blastp_subject_alignment_index) {
					if (con_str_query_alignment.length() >= itr_v_s_blastp.st_blastp_query_alignment_index) {
						con_str_query_alignment.resize(con_str_query_alignment.length() - itr_v_s_blastp.st_blastp_query_alignment_index + 1);
					}
					else {
						con_str_query_alignment.resize(0);
					}
					con_str_query_alignment += itr_v_s_blastp.str_blastp_query;
					for (size_type i = 0; i < (itr_v_s_blastp.st_blastp_query_alignment_index - 1); ++i) {
						con_str_query_alignment += ".";
					}
				}
				else {
					con_str_query_alignment += ".";
				}
				++st_index_match;
			}
			if (con_str_query_alignment.length() > itr_v_s_blastp.str_blastp_query.length() + 1) {
				con_str_query_alignment.resize(con_str_query_alignment.length() - itr_v_s_blastp.str_blastp_query.length() + 1);
			} else {
				std::cout << "\n\n ~~~ possible bad query: " << itr_v_s_blastp.str_protein;
			}
			itr_v_s_blastp.str_blastp_query_alignment = con_str_query_alignment;
			con_str_query_alignment.clear();
			st_index_match = 1;
		}
	}

	void modify_v_s_filesystem_blastp_data(filesystem_type& par_s_filesystem) {
		for (auto& itr_v_c_blastp : par_s_filesystem.v_s_blastp) {
			bool b_parse_query_accession = bool();
			string_type con_str_parse_query_accession = string_type();
			for (auto ch_blastp_query_accession : itr_v_c_blastp.str_blastp_subject_accession) {
				if ((ch_blastp_query_accession == '|') && (b_parse_query_accession)) {
					b_parse_query_accession = false;
					itr_v_c_blastp.str_blastp_subject_accession = con_str_parse_query_accession;
					con_str_parse_query_accession.clear();
					break;
				}
				if ((ch_blastp_query_accession == '|') && (!b_parse_query_accession)) {
					b_parse_query_accession = true;
				}
				if ((b_parse_query_accession) && (ch_blastp_query_accession != '|')) {
					con_str_parse_query_accession += ch_blastp_query_accession;
				}
			}
		}
	}

	inline double log_base(double d, double base) {
		return (log(d) / log(base));
	}

	void normalise_v_s_filesystem_blastp_data(filesystem_type& par_s_filesystem) {
		std::vector<blastp_type> con_v_s_blastp;
		std::vector<blastp_type> hold_v_s_blastp;
		string_type hold_str_blastp_subject = string_type();
		double hold_min_str_blastp_evalue = double();
		double hold_sum_str_blastp_evalue = double();
		bool b_hold = bool();
		for (auto itr_v_s_blastp : par_s_filesystem.v_s_blastp) {
			if (itr_v_s_blastp.str_blastp_query != hold_str_blastp_subject) {
				if (b_hold) {
					for (auto& itr_hold_v_s_blastp : hold_v_s_blastp) {
						itr_hold_v_s_blastp.d_blastp_par_prop = log_base(((double(4) * PARPROP_SCALE) / itr_hold_v_s_blastp.d_blastp_evalue), 15);
						con_v_s_blastp.push_back(itr_hold_v_s_blastp);
					}
				}
				hold_v_s_blastp.clear();
				hold_sum_str_blastp_evalue = 0;
				hold_str_blastp_subject = itr_v_s_blastp.str_blastp_query;
				hold_min_str_blastp_evalue = itr_v_s_blastp.d_blastp_evalue;
				b_hold = true;
			}
			if (b_hold && (itr_v_s_blastp.d_blastp_evalue < BLASTP_THRESHOLD)) {
				hold_v_s_blastp.push_back(itr_v_s_blastp);
				hold_sum_str_blastp_evalue += (1 / itr_v_s_blastp.d_blastp_evalue);
			}
		}
		par_s_filesystem.v_s_blastp = con_v_s_blastp;
	}

	void fout_blastp_summary(filesystem_type par_s_filesystem) {
		std::string output_blastp_summary = "blast_directory\\";
		output_blastp_summary += par_s_filesystem.str_filename;
		output_blastp_summary += "_blastp_summary.csv";
		std::ofstream fout_blastp_summary;
		fout_blastp_summary.open(output_blastp_summary);
		fout_blastp_summary << "subject,";
		fout_blastp_summary << "query_accession,";
		fout_blastp_summary << "e_value,";
		fout_blastp_summary << "par_prop,";
		fout_blastp_summary << "\n";
		for (auto itr_v_s_blastp : par_s_filesystem.v_s_blastp) {
			fout_blastp_summary << itr_v_s_blastp.str_blastp_query << ",";
			fout_blastp_summary << itr_v_s_blastp.str_blastp_subject_accession << ",";
			fout_blastp_summary << itr_v_s_blastp.d_blastp_evalue << ",";
			fout_blastp_summary << itr_v_s_blastp.d_blastp_par_prop << ",";
			fout_blastp_summary << "\n";
		}
		fout_blastp_summary.close();
	}
}

#endif