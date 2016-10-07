// * * fpf_homology_analysis.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_HOMOLOGY_ANALYSIS
#define	FPF_HOMOLOGY_ANALYSIS

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <utility>
#include <vector>

#include "fpf_data.h"
#include "fpf_filesystem.h"


namespace fpf_homology_analysis {

	using std::string;
	using std::vector;

	using fpf_data::homology_data;
	using fpf_data::multinomial;
	using fpf_data::peptide_analysis;
	using fpf_data::peptide_data;
	using fpf_data::protein_data;
	using fpf_data::protein_analysis;
	using fpf_filesystem::filesystem;
	using fpf_filesystem::sample_analysis;

	typedef vector<homology_data> v_homology_data;
	typedef vector<filesystem> v_filesystem_type;

	void create_blastp_input(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		string output_blastp_FASTA = DEFAULT_BLASTP_DIRECTORY;
		output_blastp_FASTA += par_filesystem.filename;
		output_blastp_FASTA += "_blastp_input.fasta";
		std::ofstream fout_blastp_input_FASTA;
		fout_blastp_input_FASTA.open(output_blastp_FASTA);
		size_t st_count_blastp_FASTA{};
		for (const auto& itr_v_peptide_analysis_map : par_sample_analysis.v_peptide_analysis_map) {
			++st_count_blastp_FASTA;
			fout_blastp_input_FASTA << ">" << itr_v_peptide_analysis_map.second->key_peptide_analysis;
			fout_blastp_input_FASTA << "|" << itr_v_peptide_analysis_map.second->peptide_filtered;
			fout_blastp_input_FASTA << "\n" << itr_v_peptide_analysis_map.second->peptide_filtered;
			fout_blastp_input_FASTA << "\n";
		}
		fout_blastp_input_FASTA.close();
	}

	void create_blastp_database(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		string output_homology_database_FASTA = DEFAULT_BLASTP_DIRECTORY + "blastp_database.fasta";
		std::ofstream fout_homology_database_FASTA;
		fout_homology_database_FASTA.open(output_homology_database_FASTA);
		for (auto itr_v_protein_data : par_sample_analysis.v_protein_data) {
			fout_homology_database_FASTA << ">" << itr_v_protein_data.key_protein_data;
			fout_homology_database_FASTA << "|" << itr_v_protein_data.protein_name;
			fout_homology_database_FASTA << "\n";
			for (size_t i = 0; i < itr_v_protein_data.protein_protein.length(); ++i) {
				if ((i % 60 == 0) && (i != 0)) {
					fout_homology_database_FASTA << "\n";
				}
				fout_homology_database_FASTA << itr_v_protein_data.protein_protein.at(i);
			}
			fout_homology_database_FASTA << "\n";
		}
		fout_homology_database_FASTA.close();
	}

	void create_blastp_database_refined(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		string output_homology_database_FASTA = DEFAULT_BLASTP_DIRECTORY + "blastp_database.fasta";
		std::ofstream fout_homology_database_FASTA;
		fout_homology_database_FASTA.open(output_homology_database_FASTA);
		for (auto itr_v_protein_analysis : par_sample_analysis.v_protein_analysis_selected_by_polymorphism) {
			fout_homology_database_FASTA << ">" << itr_v_protein_analysis.p_protein_data->key_protein_data;
			fout_homology_database_FASTA << "|" << itr_v_protein_analysis.p_protein_data->protein_name;
			fout_homology_database_FASTA << "\n";
			for (size_t i = 0; i < itr_v_protein_analysis.p_protein_data->protein_protein.length(); ++i) {
				if ((i % 60 == 0) && (i != 0)) {
					fout_homology_database_FASTA << "\n";
				}
				fout_homology_database_FASTA << itr_v_protein_analysis.p_protein_data->protein_protein.at(i);
			}
			fout_homology_database_FASTA << "\n";
		}
		fout_homology_database_FASTA.close();
	}

	void systemcall_blastp(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::cout << "\n\n";
		string string_system = "CD " + IgFamily::DEFAULT_BLASTP_DIRECTORY;
		string_system += " && makeblastdb.exe -in ";
		string_system += "blastp_database.fasta";
		string_system += " -dbtype prot -out FPF_blastpdb";
		system(string_system.c_str());
		string_system = "CD " + IgFamily::DEFAULT_BLASTP_DIRECTORY;
		string_system += " && blastp.exe -query ";
		string_system += par_filesystem.filename;
		string_system += "_blastp_input.fasta -db FPF_blastpdb -max_target_seqs 100 -out ";
		string_system += par_filesystem.filename;
		string_system += "_blastp_output.csv -outfmt \"10 qacc qseq sseq sacc qstart sstart qlen pident ppos mismatch qcovhsp score\"";
		system(string_system.c_str());
		string_system = "EXIT";
		system(string_system.c_str());
	}

	void parse_homology_data(filesystem& par_filesytem, sample_analysis& par_sample_analysis) {
			size_t homology_data_count_delimit = size_t(1);
			string homology_data_input = DEFAULT_BLASTP_DIRECTORY + par_filesytem.filename + "_blastp_output.csv";
			std::ifstream fin_input_blastp(homology_data_input);
			char homology_data_read{};
			vector<homology_data> temp_v_homology_data{};
			homology_data temp_homology_data{};
			string temp_parse_blastp{};
			while (fin_input_blastp.get(homology_data_read)) {
				if ((homology_data_read != ',') && (homology_data_read != '\n')) {
					temp_parse_blastp += homology_data_read;
				}
				if (homology_data_read == ',') {
					if (homology_data_count_delimit % 11 == 1) {
						temp_homology_data.blastp_query.clear();
						bool read_key_blastp_query{};
						string read_key{};
						for (const auto& itr_parse_blastp : temp_parse_blastp) {
							if (itr_parse_blastp == '|') {
								temp_homology_data.key_blastp_query = std::stoi(read_key);
								read_key_blastp_query = true;
							}
							else {
								if (read_key_blastp_query) {
									temp_homology_data.blastp_query += itr_parse_blastp;
								}
								if (!read_key_blastp_query) {
									read_key += itr_parse_blastp;
								}
							}
						}
					}
					if (homology_data_count_delimit % 11 == 2) {
						temp_homology_data.blastp_query_aligned = temp_parse_blastp;
					}
					if (homology_data_count_delimit % 11 == 3) {
						temp_homology_data.blastp_subject = temp_parse_blastp;
					}
					if (homology_data_count_delimit % 11 == 4) {
						temp_homology_data.blastp_subject_accession.clear();
						bool read_key_blastp_subject_accession{};
						string read_key{};
						for (const auto& itr_parse_blastp : temp_parse_blastp) {
							if (itr_parse_blastp == '|') {
								temp_homology_data.key_blastp_subject_accession = std::stoi(read_key);
								read_key_blastp_subject_accession = true;
							}
							else {
								if (read_key_blastp_subject_accession) {
									temp_homology_data.blastp_subject_accession += itr_parse_blastp;
								}
								if (!read_key_blastp_subject_accession) {
									read_key += itr_parse_blastp;
								}
							}
						}						
					}
					if (homology_data_count_delimit % 11 == 5) {
						temp_homology_data.blastp_query_alignment_index = std::stoi(temp_parse_blastp);
					}
					if (homology_data_count_delimit % 11 == 6) {
						temp_homology_data.blastp_subject_alignment_index = std::stoi(temp_parse_blastp);
					}
					if (homology_data_count_delimit % 11 == 10) {
						temp_homology_data.blastp_mismatch_count = std::stoi(temp_parse_blastp);
					}
					if (homology_data_count_delimit % 11 == 0) {
						temp_homology_data.blastp_query_alignment_coverage = std::stod(temp_parse_blastp);
					}
					temp_parse_blastp.clear();
					++homology_data_count_delimit;
				}
				if (homology_data_read == '\n') {
					temp_homology_data.blastp_score = std::stod(temp_parse_blastp);
					if (temp_homology_data.blastp_query_alignment_coverage >= IgFamily::HOMOLOGY_QUERY_ALIGNMENT_COVERAGE_THRESHOLD) {
						temp_v_homology_data.push_back(temp_homology_data);
					}
					temp_parse_blastp.clear();
				}
			}
			par_sample_analysis.v_homology_data = std::move(temp_v_homology_data);
	}

	void associate_homology_data_to_v_protein_data(sample_analysis& par_sample_analysis) {
		for (auto& itr_v_homology_data : par_sample_analysis.v_homology_data) {
			auto& find_protein_name = std::find_if(par_sample_analysis.v_protein_data.begin(), par_sample_analysis.v_protein_data.end(),
				[itr_v_homology_data](const protein_data& par_protein_data) {				
				return par_protein_data.key_protein_data == itr_v_homology_data.key_blastp_subject_accession;
			});
			if (find_protein_name == par_sample_analysis.v_protein_data.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_v_homology_data.blastp_subject_accession;
				string catch_error{};
				std::cin >> catch_error;
			}
			itr_v_homology_data.p_protein_data = &(*find_protein_name);
		}
	}

	void associate_homology_data_to_v_peptide_data(sample_analysis& par_sample_analysis) {
		for (auto& itr_v_homology_data : par_sample_analysis.v_homology_data) {
			auto& find_peptide_analysis = par_sample_analysis.v_peptide_analysis_map.find(itr_v_homology_data.blastp_query);
			if (find_peptide_analysis == par_sample_analysis.v_peptide_analysis_map.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_v_homology_data.blastp_subject_accession;
				string catch_error{};
				std::cin >> catch_error;
			}
			itr_v_homology_data.p_peptide_analysis = find_peptide_analysis->second;
		}
	}

	void create_protein_construct_from_protein_analysis(sample_analysis& par_sample_analysis) {
		for (auto& itr_v_homology_data : par_sample_analysis.v_homology_data) {
			auto find_protein_name = std::find_if(par_sample_analysis.v_protein_analysis_selected_by_polymorphism.begin(), par_sample_analysis.v_protein_analysis_selected_by_polymorphism.end(),
				[itr_v_homology_data](const protein_analysis& par_protein_analysis) {
				return par_protein_analysis.p_protein_data->protein_name == itr_v_homology_data.blastp_subject_accession;
			});
			if (find_protein_name == par_sample_analysis.v_protein_analysis_selected_by_polymorphism.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_v_homology_data.blastp_subject_accession;
				string catch_error{};
				std::cin >> catch_error;
			}
			itr_v_homology_data.p_protein_data = find_protein_name->p_protein_data;
		}
	}

	void create_blastp_query_alignment(sample_analysis& par_sample_analysis) {
		string temp_query_alignment{};
		size_t index_match = 1;
		for (auto& itr_v_homology_data : par_sample_analysis.v_homology_data) {
			for (auto itr_protein_protein : itr_v_homology_data.p_protein_data->protein_protein) {
				if (index_match == itr_v_homology_data.blastp_subject_alignment_index) {
					if (temp_query_alignment.length() >= itr_v_homology_data.blastp_query_alignment_index) {
						temp_query_alignment.resize(temp_query_alignment.length() - itr_v_homology_data.blastp_query_alignment_index + 1);
					}
					else {
						temp_query_alignment.resize(0);
					}
					temp_query_alignment += itr_v_homology_data.blastp_query;
					for (size_t i = 0; i < (itr_v_homology_data.blastp_query_alignment_index - 1); ++i) {
						temp_query_alignment += ".";
					}
				}
				else {
					temp_query_alignment += ".";
				}
				++index_match;
			}
			if (temp_query_alignment.length() > itr_v_homology_data.blastp_query.length() + 1) {
				temp_query_alignment.resize(temp_query_alignment.length() - itr_v_homology_data.blastp_query.length() + 1);
			} else {
				std::cout << "\n\n ~~~ possible bad query: " << itr_v_homology_data.p_protein_data->protein_protein;
			}
			itr_v_homology_data.query_alignment = temp_query_alignment;
			temp_query_alignment.clear();
			index_match = 1;
		}
	}

	void modify_homology_data(sample_analysis& par_sample_analysis) {
		for (auto& itr_v_c_blastp : par_sample_analysis.v_homology_data) {
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

	void transform_homology_data(sample_analysis& par_sample_analysis) {
		vector<homology_data> temp_v_homology_data{};
		vector<homology_data> hold_v_homology_data{};
		string hold_blastp_query{};
		double hold_min_blastp_evalue{};
		double hold_sum_blastp_evalue{};
		bool is_hold{};
		for (auto itr_v_homology_data : par_sample_analysis.v_homology_data) {
			if (itr_v_homology_data.blastp_query != hold_blastp_query) {
				if (is_hold) {
					for (auto& itr_hold_v_homology_data : hold_v_homology_data) {
						itr_hold_v_homology_data.blastp_score_transformed = std::pow(itr_hold_v_homology_data.blastp_score, IgFamily::PARAMETER_HOMOLOGY_WEIGHT);
						itr_hold_v_homology_data.blastp_score_transformed_conjugated = itr_hold_v_homology_data.blastp_score_transformed;
						temp_v_homology_data.push_back(itr_hold_v_homology_data);
					}
				}
				hold_v_homology_data.clear();
				hold_sum_blastp_evalue = 0;
				hold_blastp_query = itr_v_homology_data.blastp_query;
				hold_min_blastp_evalue = itr_v_homology_data.blastp_score;
				is_hold = true;
			}
			if (is_hold) {
				hold_v_homology_data.push_back(itr_v_homology_data);
				hold_sum_blastp_evalue += (1 / itr_v_homology_data.blastp_score);
			}
		}
		par_sample_analysis.v_homology_data = temp_v_homology_data;
	}

	void determine_homology_parameter_density(sample_analysis& par_sample_analysis, bool par_conjugated) {
		for (auto& itr_v_homology_data : par_sample_analysis.v_homology_data) {
			double temp_score_transform_conjugated_sum{};
			for (const auto& itr_v_homology_data_2 : par_sample_analysis.v_homology_data) {
				if (itr_v_homology_data.key_blastp_query == itr_v_homology_data_2.key_blastp_query) {
					temp_score_transform_conjugated_sum += itr_v_homology_data_2.blastp_score_transformed_conjugated;
				}
			}
			itr_v_homology_data.blastp_parameter_density_conjugated = (itr_v_homology_data.blastp_score_transformed_conjugated / temp_score_transform_conjugated_sum);
			itr_v_homology_data.blastp_parameter_score = (std::pow(itr_v_homology_data.blastp_parameter_density_conjugated, IgFamily::PARAMETER_SCORE_CONJUGATION_WEIGHT) * itr_v_homology_data.blastp_score_transformed_conjugated * std::pow(IgFamily::PARAMETER_SCORE_MISMATCH_WEIGHT, itr_v_homology_data.blastp_mismatch_count));
			if (!par_conjugated) {
				itr_v_homology_data.blastp_parameter_density = itr_v_homology_data.blastp_parameter_density_conjugated;
			}
		}
	}
}

#endif