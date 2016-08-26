// * * fpf_blastp_analysis.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_BLASTP_ANALYSIS 
#define	FPF_BLASTP_ANALYSIS

#include <algorithm> // provides - std::find_if
#include <cstdlib> // provides - size_t
#include <fstream> // provides - std::ifstream
#include <iostream> // provides - std::ofstream
#include <math.h> // provides - std::log10
#include <string> // provides - std::string
#include <utility> // provides - std::pair
#include <vector> // provides - std::vector

#include "fpf_data.h"
#include "fpf_filesystem.h"
#include "fpf_filesystem_analysis.h"


namespace fpf_blastp_analysis {

	using std::string;
	using std::vector;

	typedef fpf_data::blastp_data blastp_data;
	typedef vector<fpf_data::blastp_data> v_blastp_type;
	typedef fpf_data::multinomial s_multinomial_type;
	typedef fpf_data::peptide_analysis peptide_analysis;
	typedef fpf_data::peptide_data peptide_data;
	typedef fpf_data::protein_data protein_data;
	typedef fpf_data::protein_analysis protein_analysis;
	typedef fpf_filesystem::filesystem filesystem;
	typedef vector<fpf_filesystem::filesystem> v_filesystem_type;

	void create_blastp_input(filesystem par_filesystem) {
		string output_blastp_FASTA = "blast_directory\\";
		output_blastp_FASTA += par_filesystem.filename;
		output_blastp_FASTA += "_blastp_input.fasta";
		std::ofstream fout_blastp_input_FASTA;
		fout_blastp_input_FASTA.open(output_blastp_FASTA);
		size_t st_count_blastp_FASTA{};
		for (auto itr_v_peptide_analysis : par_filesystem.v_peptide_analysis) {
			++st_count_blastp_FASTA;
			fout_blastp_input_FASTA << ">" << itr_v_peptide_analysis.peptide_filtered << "\n";
			fout_blastp_input_FASTA << itr_v_peptide_analysis.peptide_filtered << "\n";
		}
		fout_blastp_input_FASTA.close();
	}

	void create_blastp_database(filesystem par_filesystem) {
		string output_blastp_database_FASTA = "blast_directory\\blastp_database.fasta";
		std::ofstream fout_blastp_database_FASTA;
		fout_blastp_database_FASTA.open(output_blastp_database_FASTA);
		for (auto itr_v_multinomial_protein : par_filesystem.v_protein_data) {
			fout_blastp_database_FASTA << ">" << itr_v_multinomial_protein.protein_name << "\n";
			for (size_t i = 0; i < itr_v_multinomial_protein.protein_protein.length(); ++i) {
				if ((i % 60 == 0) && (i != 0)) {
					fout_blastp_database_FASTA << "\n";
				}
				fout_blastp_database_FASTA << itr_v_multinomial_protein.protein_protein.at(i);
			}
			fout_blastp_database_FASTA << "\n";
		}
		fout_blastp_database_FASTA.close();
	}

	void create_blastp_database_refined(filesystem par_filesystem) {
		string output_blastp_database_FASTA = "blast_directory\\blastp_database.fasta";
		std::ofstream fout_blastp_database_FASTA;
		fout_blastp_database_FASTA.open(output_blastp_database_FASTA);
		for (auto itr_v_protein_analysis : par_filesystem.v_protein_analysis_selected_by_polymorphism) {
			fout_blastp_database_FASTA << ">" << itr_v_protein_analysis.p_protein_data->protein_name << "\n";
			for (size_t i = 0; i < itr_v_protein_analysis.p_protein_data->protein_protein.length(); ++i) {
				if ((i % 60 == 0) && (i != 0)) {
					fout_blastp_database_FASTA << "\n";
				}
				fout_blastp_database_FASTA << itr_v_protein_analysis.p_protein_data->protein_protein.at(i);
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
		string_system += std::to_string(BLASTP_EVALUE_THRESHOLD);
		string_system += " -max_target_seqs 200 -out ";
		string_system += par_filesystem.filename;
		string_system += "_blastp_output.csv -outfmt \"10 qacc qseq sseq sacc qstart sstart qlen pident ppos mismatch evalue\"";
		system(string_system.c_str());
		string_system = "EXIT";
		system(string_system.c_str());
	}

	void create_v_blastp_data(filesystem& par_filesystem) {
		if (par_filesystem.denovopeptides_exist) {
			size_t blastp_data_count_delimit = size_t(1);
			string blastp_data_output = "blast_directory\\" + par_filesystem.filename + "_blastp_output.csv";
			std::ifstream fin_input_blastp(blastp_data_output);
			char blastp_data_read{};
			vector<blastp_data> temp_v_blastp_data{};
			blastp_data temp_blastp_data{};
			string temp_parse_blastp{};
			while (fin_input_blastp.get(blastp_data_read)) {
				if ((blastp_data_read != ',') && (blastp_data_read != '\n')) {
					temp_parse_blastp += blastp_data_read;
				}
				if (blastp_data_read == ',') {
					if (blastp_data_count_delimit % 10 == 1) {
						temp_blastp_data.blastp_query = temp_parse_blastp;
					}
					if (blastp_data_count_delimit % 10 == 2) {
						temp_blastp_data.blastp_query_aligned = temp_parse_blastp;
					}
					if (blastp_data_count_delimit % 10 == 3) {
						temp_blastp_data.blastp_subject = temp_parse_blastp;
					}
					if (blastp_data_count_delimit % 10 == 4) {
						temp_blastp_data.blastp_subject_accession = temp_parse_blastp;
					}
					if (blastp_data_count_delimit % 10 == 5) {
						temp_blastp_data.blastp_query_alignment_index = std::stoi(temp_parse_blastp);
					}
					if (blastp_data_count_delimit % 10 == 6) {
						temp_blastp_data.blastp_subject_alignment_index = std::stoi(temp_parse_blastp);
					}
					temp_parse_blastp.clear();
					++blastp_data_count_delimit;
				}
				if (blastp_data_read == '\n') {
					temp_blastp_data.blastp_evalue = std::stod(temp_parse_blastp);
					temp_v_blastp_data.push_back(temp_blastp_data);
					temp_parse_blastp.clear();
				}
			}
			par_filesystem.v_blastp_data = std::move(temp_v_blastp_data);
		}
	}

	void associate_blastp_data_to_v_protein_data(filesystem& par_filesystem) {
		for (auto& itr_v_blastp_data : par_filesystem.v_blastp_data) {
			auto& find_protein_name = std::find_if(par_filesystem.v_protein_data.begin(), par_filesystem.v_protein_data.end(),
				[itr_v_blastp_data](protein_data par_protein_data) {
				return par_protein_data.protein_name == itr_v_blastp_data.blastp_subject_accession;
			});
			if (find_protein_name == par_filesystem.v_protein_data.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_v_blastp_data.blastp_subject_accession;
				string catch_error{};
				std::cin >> catch_error;
			}
			itr_v_blastp_data.p_protein_data = &(*find_protein_name);
		}
	}

	void associate_blastp_data_to_v_peptide_data(filesystem& par_filesystem) {
		for (auto& itr_v_blastp_data : par_filesystem.v_blastp_data) {
			auto& find_peptide_analysis = std::find_if(par_filesystem.v_peptide_analysis.begin(), par_filesystem.v_peptide_analysis.end(),
				[itr_v_blastp_data](peptide_analysis par_peptide_analysis) {
				return par_peptide_analysis.peptide_filtered == itr_v_blastp_data.blastp_query;
			});
			if (find_peptide_analysis == par_filesystem.v_peptide_analysis.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_v_blastp_data.blastp_subject_accession;
				string catch_error{};
				std::cin >> catch_error;
			}
			itr_v_blastp_data.p_peptide_analysis = &(*find_peptide_analysis);
		}
	}

	void create_protein_from_protein_analysis(filesystem& par_filesystem) {
		for (auto& itr_v_blastp_data : par_filesystem.v_blastp_data) {
			auto find_protein_name = std::find_if(par_filesystem.v_protein_analysis_selected_by_polymorphism.begin(), par_filesystem.v_protein_analysis_selected_by_polymorphism.end(),
				[itr_v_blastp_data](protein_analysis par_protein_analysis) {
				return par_protein_analysis.p_protein_data->protein_name == itr_v_blastp_data.blastp_subject_accession;
			});
			if (find_protein_name == par_filesystem.v_protein_analysis_selected_by_polymorphism.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_v_blastp_data.blastp_subject_accession;
				string catch_error{};
				std::cin >> catch_error;
			}
			itr_v_blastp_data.p_protein_data = find_protein_name->p_protein_data;
		}
	}

	void create_query_alignment(filesystem& par_filesystem) {
		string temp_query_alignment{};
		size_t index_match = 1;
		for (auto& itr_v_blastp_data : par_filesystem.v_blastp_data) {
			for (auto itr_protein_protein : itr_v_blastp_data.p_protein_data->protein_protein) {
				if (index_match == itr_v_blastp_data.blastp_subject_alignment_index) {
					if (temp_query_alignment.length() >= itr_v_blastp_data.blastp_query_alignment_index) {
						temp_query_alignment.resize(temp_query_alignment.length() - itr_v_blastp_data.blastp_query_alignment_index + 1);
					}
					else {
						temp_query_alignment.resize(0);
					}
					temp_query_alignment += itr_v_blastp_data.blastp_query;
					for (size_t i = 0; i < (itr_v_blastp_data.blastp_query_alignment_index - 1); ++i) {
						temp_query_alignment += ".";
					}
				}
				else {
					temp_query_alignment += ".";
				}
				++index_match;
			}
			if (temp_query_alignment.length() > itr_v_blastp_data.blastp_query.length() + 1) {
				temp_query_alignment.resize(temp_query_alignment.length() - itr_v_blastp_data.blastp_query.length() + 1);
			} else {
				std::cout << "\n\n ~~~ possible bad query: " << itr_v_blastp_data.p_protein_data->protein_protein;
			}
			itr_v_blastp_data.query_alignment = temp_query_alignment;
			temp_query_alignment.clear();
			index_match = 1;
		}
	}

	void modify_filesystem_blastp_data(filesystem& par_filesystem) {
		for (auto& itr_v_c_blastp : par_filesystem.v_blastp_data) {
			bool b_parse_query_accession{};
			string temp_parse_query_accession{};
			for (auto ch_blastp_query_accession : itr_v_c_blastp.blastp_subject_accession) {
				if ((ch_blastp_query_accession == '|') && (b_parse_query_accession)) {
					b_parse_query_accession = false;
					itr_v_c_blastp.blastp_subject_accession = temp_parse_query_accession;
					temp_parse_query_accession.clear();
					break;
				}
				if ((ch_blastp_query_accession == '|') && (!b_parse_query_accession)) {
					b_parse_query_accession = true;
				}
				if ((b_parse_query_accession) && (ch_blastp_query_accession != '|')) {
					temp_parse_query_accession += ch_blastp_query_accession;
				}
			}
		}
	}

	inline double log_base(double d, double base) {
		return (log(d) / log(base));
	}

	void normalise_blastp_data(filesystem& par_filesystem) {
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
						itr_hold_v_s_blastp.blastp_evalue_transformed = log_base(((double(3) * BLASTP_PARPROP_SCALE) / itr_hold_v_s_blastp.blastp_evalue), 1.2);
						//itr_hold_v_s_blastp.blastp_evalue_transformed = 1 / (itr_hold_v_s_blastp.blastp_evalue + 0.01);
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
			itr_v_blastp.blastp_parameter_score = (itr_v_blastp.blastp_parameter_density * itr_v_blastp.blastp_evalue_transformed * (itr_v_blastp.p_peptide_analysis->v_denovo_peptide_averagescore / 100));
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