// * * fpf_blastp_analysis.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_BLASTP_ANALYSIS 
#define	FPF_BLASTP_ANALYSIS

#include <cstdlib> // provides - size_t
#include <string> // provides - string
#include <vector> // provides - vector
#include <fstream> // provides - std::ifstream
#include <iostream> // provides - std::ofstream
#include <utility> // provides - pair
#include <math.h> // provides - log10
#include <algorithm> // provides - std::find_if

#include "fpf_filesystem.h"
#include "fpf_filesystem_analysis.h"
#include "fpf_data.h"


namespace fpf_blastp_analysis {

	using std::string;
	using std::vector;

	typedef fpf_data::FASTA_category FASTA_category;
	typedef fpf_data::blastp_data blastp_data;
	typedef vector<fpf_data::blastp_data> v_s_blastp_tpye;
	typedef fpf_data::multinomial s_multinomial_type;
	typedef fpf_data::category_analysis category_analysis;
	typedef fpf_filesystem::filesystem filesystem;
	typedef vector<fpf_filesystem::filesystem> v_filesystem_type;

	void create_blastp_input(filesystem par_filesystem) {
		string output_blastp_FASTA = "blast_directory\\";
		output_blastp_FASTA += par_filesystem.filename;
		output_blastp_FASTA += "_blastp_input.fasta";
		std::ofstream fout_blastp_input_FASTA;
		fout_blastp_input_FASTA.open(output_blastp_FASTA);
		size_t st_count_blastp_FASTA{};
		for (auto itr_v_s_peptide_data : par_filesystem.v_peptide_data) {
			++st_count_blastp_FASTA;
			fout_blastp_input_FASTA << ">" << itr_v_s_peptide_data.peptide_filtered << "\n";
			fout_blastp_input_FASTA << itr_v_s_peptide_data.peptide_filtered << "\n";
		}
		fout_blastp_input_FASTA.close();
	}

	void create_blastp_database(filesystem par_filesystem) {
		string output_blastp_database_FASTA = "blast_directory\\blastp_database.fasta";
		std::ofstream fout_blastp_database_FASTA;
		fout_blastp_database_FASTA.open(output_blastp_database_FASTA);
		for (auto itr_v_multinomial_category : par_filesystem.v_multinomial_category) {
			fout_blastp_database_FASTA << ">" << itr_v_multinomial_category.category_name << "\n";
			for (size_t i = 0; i < itr_v_multinomial_category.category_protein.length(); ++i) {
				if ((i % 60 == 0) && (i != 0)) {
					fout_blastp_database_FASTA << "\n";
				}
				fout_blastp_database_FASTA << itr_v_multinomial_category.category_protein.at(i);
			}
			fout_blastp_database_FASTA << "\n";
		}
		fout_blastp_database_FASTA.close();
	}

	void create_blastp_database_refined(filesystem par_filesystem) {
		string output_blastp_database_FASTA = "blast_directory\\blastp_database.fasta";
		std::ofstream fout_blastp_database_FASTA;
		fout_blastp_database_FASTA.open(output_blastp_database_FASTA);
		for (auto itr_v_category_analysis : par_filesystem.v_category_analysis_selected_by_polymorphism) {
			fout_blastp_database_FASTA << ">" << itr_v_category_analysis.p_FASTA_category->category_name << "\n";
			for (size_t i = 0; i < itr_v_category_analysis.p_FASTA_category->category_protein.length(); ++i) {
				if ((i % 60 == 0) && (i != 0)) {
					fout_blastp_database_FASTA << "\n";
				}
				fout_blastp_database_FASTA << itr_v_category_analysis.p_FASTA_category->category_protein.at(i);
			}
			fout_blastp_database_FASTA << "\n";
		}
		fout_blastp_database_FASTA.close();
	}

	void sys_blastp(filesystem par_filesystem) {
		std::cout << "\n\n";
		string string_system = "CD C:\\Users\\LJ\\IgFamily\\blast_directory\\";
		string_system += " && makeblastdb.exe -in ";
		string_system += "blastp_database.fasta";
		string_system += " -dbtype prot -out FPF_blastpdb";
		system(string_system.c_str());
		string_system = "CD C:\\Users\\LJ\\IgFamily\\blast_directory\\";
		string_system += " && blastp.exe -query ";
		string_system += par_filesystem.filename;
		string_system += "_blastp_input.fasta -db FPF_blastpdb -evalue ";
		string_system += std::to_string(EVALUE_THRESHOLD);
		string_system += " -max_target_seqs 200 -out ";
		string_system += par_filesystem.filename;
		string_system += "_blastp_output.csv -outfmt \"10 qacc qseq sseq sacc qstart sstart qlen pident ppos mismatch evalue\"";
		system(string_system.c_str());
		string_system = "EXIT";
		system(string_system.c_str());
	}

	void create_v_s_blastp(filesystem& par_filesystem) {
		if (par_filesystem.denovopeptides_exist) {
			size_t blastp_data_count_delimit = size_t(1);
			string blastp_data_output = "blast_directory\\" + par_filesystem.filename + "_blastp_output.csv";
			std::ifstream fin_input_blastp(blastp_data_output);
			char blastp_data_read{};
			vector<blastp_data> temp_v_blastp_data{};
			blastp_data temp_blastp_data{};
			string con_str_parse_blastp{};
			while (fin_input_blastp.get(blastp_data_read)) {
				if ((blastp_data_read != ',') && (blastp_data_read != '\n')) {
					con_str_parse_blastp += blastp_data_read;
				}
				if (blastp_data_read == ',') {
					if (blastp_data_count_delimit % 10 == 1) {
						temp_blastp_data.blastp_query = con_str_parse_blastp;
					}
					if (blastp_data_count_delimit % 10 == 2) {
						temp_blastp_data.blastp_query_aligned = con_str_parse_blastp;
					}
					if (blastp_data_count_delimit % 10 == 3) {
						temp_blastp_data.blastp_subject = con_str_parse_blastp;
					}
					if (blastp_data_count_delimit % 10 == 4) {
						temp_blastp_data.blastp_subject_accession = con_str_parse_blastp;
					}
					if (blastp_data_count_delimit % 10 == 5) {
						temp_blastp_data.blastp_query_alignment_index = std::stoi(con_str_parse_blastp);
					}
					if (blastp_data_count_delimit % 10 == 6) {
						temp_blastp_data.blastp_subject_alignment_index = std::stoi(con_str_parse_blastp);
					}
					con_str_parse_blastp.clear();
					++blastp_data_count_delimit;
				}
				if (blastp_data_read == '\n') {
					temp_blastp_data.blastp_evalue = std::stod(con_str_parse_blastp);
					temp_v_blastp_data.push_back(temp_blastp_data);
					con_str_parse_blastp.clear();
				}
			}
			par_filesystem.v_blastp_data = std::move(temp_v_blastp_data);
		}
	}

	void create_str_protein(filesystem& par_filesystem) {
		for (auto& itr_v_s_blastp : par_filesystem.v_blastp_data) {
			auto& find_str_genefamily = std::find_if(par_filesystem.v_multinomial_category.begin(), par_filesystem.v_multinomial_category.end(),
				[itr_v_s_blastp](FASTA_category par_s_multinomial_element_data) {
				return par_s_multinomial_element_data.category_name == itr_v_s_blastp.blastp_subject_accession;
			});
			if (find_str_genefamily == par_filesystem.v_multinomial_category.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_v_s_blastp.blastp_subject_accession;
				string catch_error;
				std::cin >> catch_error;
			}
			itr_v_s_blastp.p_FASTA_category = &(*find_str_genefamily);
			itr_v_s_blastp.blastp_subject_accession_class = find_str_genefamily->category_class;

		}
	}

	void create_str_protein_from_category_analysis(filesystem& par_filesystem) {
		for (auto& itr_v_blastp_data : par_filesystem.v_blastp_data) {
			auto find_category_name = std::find_if(par_filesystem.v_category_analysis_selected_by_polymorphism.begin(), par_filesystem.v_category_analysis_selected_by_polymorphism.end(),
				[itr_v_blastp_data](category_analysis par_category_analysis) {
				return par_category_analysis.p_FASTA_category->category_name == itr_v_blastp_data.blastp_subject_accession;
			});
			if (find_category_name == par_filesystem.v_category_analysis_selected_by_polymorphism.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_v_blastp_data.blastp_subject_accession;
				string catch_error{};
				std::cin >> catch_error;
			}
			itr_v_blastp_data.p_FASTA_category = find_category_name->p_FASTA_category;

		}
	}

	void create_str_query_alignment(filesystem& par_filesystem) {
		string con_str_query_alignment{};
		size_t st_index_match = 1;
		for (auto& itr_v_s_blastp : par_filesystem.v_blastp_data) {
			for (auto itr_str_protein : itr_v_s_blastp.p_FASTA_category->category_protein) {
				if (st_index_match == itr_v_s_blastp.blastp_subject_alignment_index) {
					if (con_str_query_alignment.length() >= itr_v_s_blastp.blastp_query_alignment_index) {
						con_str_query_alignment.resize(con_str_query_alignment.length() - itr_v_s_blastp.blastp_query_alignment_index + 1);
					}
					else {
						con_str_query_alignment.resize(0);
					}
					con_str_query_alignment += itr_v_s_blastp.blastp_query;
					for (size_t i = 0; i < (itr_v_s_blastp.blastp_query_alignment_index - 1); ++i) {
						con_str_query_alignment += ".";
					}
				}
				else {
					con_str_query_alignment += ".";
				}
				++st_index_match;
			}
			if (con_str_query_alignment.length() > itr_v_s_blastp.blastp_query.length() + 1) {
				con_str_query_alignment.resize(con_str_query_alignment.length() - itr_v_s_blastp.blastp_query.length() + 1);
			} else {
				std::cout << "\n\n ~~~ possible bad query: " << itr_v_s_blastp.p_FASTA_category->category_protein;
			}
			itr_v_s_blastp.query_alignment = con_str_query_alignment;
			con_str_query_alignment.clear();
			st_index_match = 1;
		}
	}

	void modify_v_s_filesystem_blastp_data(filesystem& par_filesystem) {
		for (auto& itr_v_c_blastp : par_filesystem.v_blastp_data) {
			bool b_parse_query_accession{};
			string con_str_parse_query_accession{};
			for (auto ch_blastp_query_accession : itr_v_c_blastp.blastp_subject_accession) {
				if ((ch_blastp_query_accession == '|') && (b_parse_query_accession)) {
					b_parse_query_accession = false;
					itr_v_c_blastp.blastp_subject_accession = con_str_parse_query_accession;
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

	void normalise_v_s_filesystem_blastp_data(filesystem& par_filesystem) {
		vector<blastp_data> temp_v_blastp_data{};
		vector<blastp_data> hold_v_blastp_data{};
		string hold_blastp_subject{};
		double hold_min_blastp_evalue{};
		double hold_sum_blastp_evalue{};
		bool is_hold{};
		for (auto itr_v_s_blastp : par_filesystem.v_blastp_data) {
			if (itr_v_s_blastp.blastp_query != hold_blastp_subject) {
				if (is_hold) {
					for (auto& itr_hold_v_s_blastp : hold_v_blastp_data) {
						itr_hold_v_s_blastp.blastp_evalue_transformed = log_base(((double(3) * PARPROP_SCALE) / itr_hold_v_s_blastp.blastp_evalue), 1.4);
						//itr_hold_v_s_blastp.blastp_evalue_transformed = (itr_hold_v_s_blastp.blastp_evalue + 0.01);
						temp_v_blastp_data.push_back(itr_hold_v_s_blastp);
					}
				}
				hold_v_blastp_data.clear();
				hold_sum_blastp_evalue = 0;
				hold_blastp_subject = itr_v_s_blastp.blastp_query;
				hold_min_blastp_evalue = itr_v_s_blastp.blastp_evalue;
				is_hold = true;
			}
			if (is_hold && (itr_v_s_blastp.blastp_evalue < BLASTP_THRESHOLD)) {
				hold_v_blastp_data.push_back(itr_v_s_blastp);
				hold_sum_blastp_evalue += (1 / itr_v_s_blastp.blastp_evalue);
			}
		}
		par_filesystem.v_blastp_data = temp_v_blastp_data;
	}

	void determine_blastp_parameter_density(filesystem& par_filesystem) {
		for (auto& itr_v_blastp : par_filesystem.v_blastp_data) {
			double temp_evalue_transform_sum{};
			for (const auto& itr_v_blastp_2 : par_filesystem.v_blastp_data) {
				if (itr_v_blastp.blastp_query == itr_v_blastp_2.blastp_query) {
					temp_evalue_transform_sum += itr_v_blastp_2.blastp_evalue_transformed;
				}
			}
			itr_v_blastp.blastp_parameter_density = (itr_v_blastp.blastp_evalue_transformed / temp_evalue_transform_sum);
			itr_v_blastp.blastp_parameter_score = (itr_v_blastp.blastp_parameter_density * itr_v_blastp.blastp_evalue_transformed);
		}
	}

	void fout_blastp_summary(filesystem par_filesystem) {
		string output_blastp_summary = "blast_directory\\";
		output_blastp_summary += par_filesystem.filename;
		output_blastp_summary += "_blastp_summary.csv";
		std::ofstream fout_blastp_summary;
		fout_blastp_summary.open(output_blastp_summary);
		fout_blastp_summary << "subject,";
		fout_blastp_summary << "query_accession,";
		fout_blastp_summary << "e_value,";
		fout_blastp_summary << "par_prop,";
		fout_blastp_summary << "\n";
		for (auto itr_v_s_blastp : par_filesystem.v_blastp_data) {
			fout_blastp_summary << itr_v_s_blastp.blastp_query << ",";
			fout_blastp_summary << itr_v_s_blastp.blastp_subject_accession << ",";
			fout_blastp_summary << itr_v_s_blastp.blastp_evalue << ",";
			fout_blastp_summary << itr_v_s_blastp.blastp_evalue_transformed << ",";
			fout_blastp_summary << itr_v_s_blastp.blastp_parameter_density << ",";
			fout_blastp_summary << itr_v_s_blastp.blastp_parameter_score << ",";
			fout_blastp_summary << "\n";
		}
		fout_blastp_summary.close();
	}
}

#endif