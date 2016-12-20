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

	using fpf_data::denovo_peptide;
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
		for (auto itr_v_protein_analysis : par_sample_analysis.v_protein_analysis_with_selected_polymorphism) {
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
						temp_homology_data.alignment_coverage = std::stod(temp_parse_blastp);
					}
					temp_parse_blastp.clear();
					++homology_data_count_delimit;
				}
				if (homology_data_read == '\n') {
					temp_homology_data.blastp_homology = (IgFamily::PARAMETER_HOMOLOGY_SCALING_FACTOR * std::stod(temp_parse_blastp));
					temp_homology_data.alignment_coverage_delta = temp_homology_data.blastp_query.length() - temp_homology_data.blastp_query_aligned.length();
					temp_v_homology_data.push_back(temp_homology_data);
					temp_parse_blastp.clear();
				}
			}
			par_sample_analysis.v_homology_data = std::move(temp_v_homology_data);
	}

	void create_blastp_query_alignment(sample_analysis& par_sample_analysis) {
		string temp_query_alignment{};
		size_t index_match = 1;
		for (auto& itr_v_homology_data : par_sample_analysis.v_homology_data) {
			for (const auto& itr_protein_protein : itr_v_homology_data.p_protein_data->protein_protein) {
				if (index_match == itr_v_homology_data.blastp_subject_alignment_index) {
					if (temp_query_alignment.length() >= itr_v_homology_data.blastp_query_alignment_index) 
					{
						temp_query_alignment.resize(temp_query_alignment.length() - itr_v_homology_data.blastp_query_alignment_index + 1);
					}
					else 
					{
						temp_query_alignment.resize(0);
					}
					temp_query_alignment += itr_v_homology_data.blastp_query;
					for (size_t i = 0; i < (itr_v_homology_data.blastp_query_alignment_index - 1); ++i) 
					{
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
			}
			else {
				std::cout << "\n\n ~~~ possible bad query: " << itr_v_homology_data.p_protein_data->protein_protein;
			}
			itr_v_homology_data.alignment = temp_query_alignment;
			temp_query_alignment.clear();
			index_match = 1;
			size_t count_align{};
			size_t align_mismatch{};
			size_t count_mismatch{};
			for (const auto& itr_alignment : itr_v_homology_data.alignment) {
				if (itr_alignment != '.') {
					if (itr_v_homology_data.p_protein_data->protein_protein[count_align] != itr_v_homology_data.blastp_query[align_mismatch]) {
						if (!((itr_v_homology_data.p_protein_data->protein_protein[count_align] == 'I')
							&& (itr_v_homology_data.blastp_query[align_mismatch] == 'L')
							|| ((itr_v_homology_data.p_protein_data->protein_protein[count_align] == 'L')
								&& (itr_v_homology_data.blastp_query[align_mismatch] == 'I')))) {
							++count_mismatch;
						}
					}
					++align_mismatch;
				}
				++count_align;
			}
			itr_v_homology_data.blastp_mismatch_count = count_mismatch;
		}
	}

	void determine_homology_data_uniquenesss(sample_analysis& par_sample_analysis)
	{
		for (auto& itr_v_homology_data : par_sample_analysis.v_homology_data)
		{
			for (auto& itr_v_homology_data_2 : par_sample_analysis.v_homology_data)
			{
				if ((itr_v_homology_data.blastp_query == itr_v_homology_data_2.blastp_query)
					&& (itr_v_homology_data.blastp_mismatch_count == 0))
				{
					++itr_v_homology_data.blastp_sharedidentical_count;
				}
			}
		}
	}

	void associate_homology_data_to_protein_data(sample_analysis& par_sample_analysis) {
		for (auto& itr_v_homology_data : par_sample_analysis.v_homology_data) {
			auto& find_protein_data = std::find_if(par_sample_analysis.v_protein_data.begin(), par_sample_analysis.v_protein_data.end(),
				[itr_v_homology_data](const protein_data& par_protein_data) {				
				return par_protein_data.key_protein_data == itr_v_homology_data.key_blastp_subject_accession;
			});
			if (find_protein_data == par_sample_analysis.v_protein_data.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_v_homology_data.blastp_subject_accession;
				string catch_error{};
				std::cin >> catch_error;
			}
			itr_v_homology_data.p_protein_data = &(*find_protein_data);
		}
	}

	void associate_homology_data_to_v_peptide_data(sample_analysis& par_sample_analysis) {
		denovo_peptide default_denovo_peptide{};
		for (auto& itr_homology_data : par_sample_analysis.v_homology_data) {
			auto& find_peptide_analysis = par_sample_analysis.v_peptide_analysis_map.find(itr_homology_data.blastp_query);
			if (find_peptide_analysis == par_sample_analysis.v_peptide_analysis_map.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_homology_data.blastp_subject_accession;
				string catch_error{};
				std::cin >> catch_error;
			}
			itr_homology_data.p_peptide_analysis = find_peptide_analysis->second;
			itr_homology_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence = &default_denovo_peptide;
			for (auto& itr_peptide_data : find_peptide_analysis->second->v_peptide_data) {
				++itr_homology_data.denovo_replicate_count;
				if (itr_homology_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average < itr_peptide_data->denovo_peptide_data_filtered.localconfidence_average) {
					itr_homology_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence = &itr_peptide_data->denovo_peptide_data_filtered;
				}
			}
		}
	}

	void associate_homology_data_to_protein_analysis_refined(sample_analysis& par_sample_analysis) {
		for (auto& itr_v_homology_data : par_sample_analysis.v_homology_data) {
			auto find_protein_data = std::find_if(par_sample_analysis.v_protein_analysis_with_selected_polymorphism.begin(), par_sample_analysis.v_protein_analysis_with_selected_polymorphism.end(),
				[itr_v_homology_data](const protein_analysis& par_protein_analysis) {
				return par_protein_analysis.p_protein_data->protein_name == itr_v_homology_data.blastp_subject_accession;
			});
			if (find_protein_data == par_sample_analysis.v_protein_analysis_with_selected_polymorphism.end()) {
				std::cout << "\n\n error - std::find_if returns nullptr";
				std::cout << "\n\n" << itr_v_homology_data.blastp_subject_accession;
				string catch_error{};
				std::cin >> catch_error;
			}
			itr_v_homology_data.p_protein_data = find_protein_data->p_protein_data;
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

	void transform_homology_data(sample_analysis& par_sample_analysis) {
		for (auto& itr_v_homology_data : par_sample_analysis.v_homology_data) {
			itr_v_homology_data.blastp_homology_transformed = itr_v_homology_data.blastp_homology * std::pow(IgFamily::PARAMETER_HOMOLOGY_DELTA_ALIGNMENT_WEIGHT, itr_v_homology_data.alignment_coverage_delta);
			itr_v_homology_data.blastp_homology_transformed = itr_v_homology_data.blastp_homology_transformed * std::pow(IgFamily::PARAMETER_HOMOLOGY_MISMATCH_WEIGHT, itr_v_homology_data.blastp_mismatch_count);
			if (itr_v_homology_data.blastp_homology == 0) {
				itr_v_homology_data.blastp_homology = double(1);
			}
			itr_v_homology_data.blastp_homology_transformed = std::pow(itr_v_homology_data.blastp_homology_transformed, IgFamily::PARAMETER_HOMOLOGY_WEIGHT);
			itr_v_homology_data.blastp_homology_transformed_conjugated = itr_v_homology_data.blastp_homology_transformed;
		}
	}

	void determine_HomologyDataParameters(sample_analysis& par_sample_analysis, bool par_conjugated) {
		if (par_conjugated) {
			for (auto& itr_homology_data : par_sample_analysis.v_homology_data) {
				itr_homology_data.blastp_homology_transformed_conjugated *= std::pow(itr_homology_data.blastp_homology_transformed, double(IgFamily::PARAMETER_PRIOR_DISTRIBUTION_WEIGHT));
			}
		}
		for (auto& itr_homology_data : par_sample_analysis.v_homology_data) {
			double temp_score_transform_conjugated_sum{};
			for (auto& itr_homology_data_2 : par_sample_analysis.v_homology_data) {
				if (itr_homology_data.key_blastp_query == itr_homology_data_2.key_blastp_query) {
					temp_score_transform_conjugated_sum += itr_homology_data_2.blastp_homology_transformed_conjugated;
				}
			}
			if (temp_score_transform_conjugated_sum != double(0)) {
				itr_homology_data.blastp_homology_density_conjugated = (itr_homology_data.blastp_homology_transformed_conjugated / temp_score_transform_conjugated_sum);
			}
			if (!par_conjugated) {
				itr_homology_data.blastp_homology_density = itr_homology_data.blastp_homology_density_conjugated;
			}
			itr_homology_data.score = (itr_homology_data.blastp_homology_density_conjugated * itr_homology_data.denovo_replicate_count);
		}
	}

	inline bool predicate_v_homology_data_by_homology_distribution(const homology_data* i, const homology_data* j)
	{
		return ((i->blastp_homology_density_conjugated) > (j->blastp_homology_density_conjugated));
	}
	 
	void sort_v_homology_data_by_homology_distribution(vector<homology_data*>& par_v_homology_data) {
		std::sort(par_v_homology_data.begin(), par_v_homology_data.end(), predicate_v_homology_data_by_homology_distribution);
	}

	void aggregate_v_homology_data_by_homology_distribution(sample_analysis& par_sample_analysis) {
		for (auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis) {
			for (auto& itr_homology_data : itr_protein_analysis.v_homology_data_combined_by_protein) {
				itr_homology_data.v_homology_data_aggregated_by_homology_distribution.push_back(&itr_homology_data);
				for (auto& itr_protein_analysis_2 : par_sample_analysis.v_protein_analysis) {
					for (auto& itr_homology_data_2 : itr_protein_analysis_2.v_homology_data_combined_by_protein) {
						if (itr_protein_analysis_2.p_protein_data->protein_name != itr_protein_analysis.p_protein_data->protein_name) {
							if (itr_homology_data.key_blastp_query == itr_homology_data_2.key_blastp_query) {
								if (itr_homology_data_2.blastp_homology_density_conjugated >= IgFamily::REPORT_V_HOMOLOGY_DATA_AGGREGATED_BY_PROTEIN_CONJUGATED_DENSITY_THRESHOLD) {
									itr_homology_data.v_homology_data_aggregated_by_homology_distribution.push_back(&itr_homology_data_2);
								}
							}
						}
					}
				}
				sort_v_homology_data_by_homology_distribution(itr_homology_data.v_homology_data_aggregated_by_homology_distribution);
			}			
		}
	}

	void normalise_v_HomologyData(sample_analysis& par_sample_analysis) {
		for (auto& itr_homology_data : par_sample_analysis.v_homology_data) {
			itr_homology_data.blastp_homology_transformed = std::pow(itr_homology_data.blastp_homology_transformed, (double(1) / IgFamily::PARAMETER_HOMOLOGY_WEIGHT));
			itr_homology_data.blastp_homology_transformed_conjugated = std::pow(itr_homology_data.blastp_homology_transformed_conjugated, (double(1) / IgFamily::PARAMETER_HOMOLOGY_WEIGHT));
		}
		for (auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis) {
			for (auto& itr_homology_data : itr_protein_analysis.v_homology_data_combined_by_protein) {
				itr_homology_data.blastp_homology_transformed = std::pow(itr_homology_data.blastp_homology_transformed, (double(1) / IgFamily::PARAMETER_HOMOLOGY_WEIGHT));
				itr_homology_data.blastp_homology_transformed_conjugated = std::pow(itr_homology_data.blastp_homology_transformed_conjugated, (double(1) / IgFamily::PARAMETER_HOMOLOGY_WEIGHT));
			}
		}
	}
}

#endif