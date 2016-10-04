// * * fpf_data_analysis.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_DATA_ANALYSIS
#define	FPF_DATA_ANALYSIS

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include <string>
#include <utility>
#include <vector>

#include "fpf_data.h"
#include "fpf_filesystem.h"
#include "fpf_homology_analysis.h"
#include "fpf_report.h"


namespace fpf_data_analysis {

	using std::map;
	using std::string;
	using std::vector;

	using fpf_data::homology_data;
	using fpf_data::denovo_peptide;
	using fpf_data::protein_data;
	using fpf_data::peptide_analysis;
	using fpf_data::peptide_data;
	using fpf_data::protein_analysis;
	using fpf_data::multinomial;
	using fpf_data::proteinconstruct_aminoacid;
	using fpf_filesystem::filesystem;
	using fpf_filesystem::sample_analysis;
	using fpf_parse::csv_data;
	using fpf_parse::FASTA_data;

	map<string, peptide_analysis*> create_v_peptide_analysis_map(vector<peptide_analysis>& par_v_peptide_analysis) {
		map<string, peptide_analysis*> temp_v_peptide_analysis_map{};
		for (auto& itr_v_peptide_analysis : par_v_peptide_analysis) {
			temp_v_peptide_analysis_map[itr_v_peptide_analysis.peptide_filtered] = &itr_v_peptide_analysis;
		}
		return temp_v_peptide_analysis_map;
	}

	vector<peptide_analysis> create_v_peptide_analysis(vector<peptide_data>& par_v_peptide_data) {
		vector<peptide_analysis> temp_v_peptide_analysis{};
		size_t temp_key_peptide_analysis{};
		for (auto& itr_v_peptide_data : par_v_peptide_data) {
			peptide_analysis temp_peptide_analysis{};
			auto& find_peptide_analysis = std::find_if(temp_v_peptide_analysis.begin(), temp_v_peptide_analysis.end(),
				[itr_v_peptide_data](const peptide_analysis& par_peptide_analysis) {
				return par_peptide_analysis.peptide_filtered == itr_v_peptide_data.peptide_filtered;
			});
			if (find_peptide_analysis == temp_v_peptide_analysis.end()) {
				temp_peptide_analysis.peptide_filtered = itr_v_peptide_data.peptide_filtered;
				++temp_peptide_analysis.replicate_count;
				temp_peptide_analysis.v_peptide_data.push_back(&itr_v_peptide_data);
				temp_peptide_analysis.v_peptide_withoutmod_mz.push_back(itr_v_peptide_data.peptide_mz);
				temp_peptide_analysis.v_peptide_withoutmod_z.push_back(itr_v_peptide_data.peptide_z);
				temp_peptide_analysis.v_peptide_withoutmod_rt.push_back(itr_v_peptide_data.peptide_rt);
				temp_peptide_analysis.v_peptide_withoutmod_m.push_back(itr_v_peptide_data.peptide_m);
				temp_peptide_analysis.v_denovo_peptide_averagescore = itr_v_peptide_data.denovo_peptide_data_filtered.localconfidence_average;
				temp_peptide_analysis.key_peptide_analysis = temp_key_peptide_analysis;
				temp_v_peptide_analysis.push_back(temp_peptide_analysis);
				++temp_key_peptide_analysis;
			}
			else {
				++find_peptide_analysis->replicate_count;
				find_peptide_analysis->v_peptide_data.push_back(&itr_v_peptide_data);
				find_peptide_analysis->v_peptide_withoutmod_mz.push_back(itr_v_peptide_data.peptide_mz);
				find_peptide_analysis->v_peptide_withoutmod_z.push_back(itr_v_peptide_data.peptide_z);
				find_peptide_analysis->v_peptide_withoutmod_rt.push_back(itr_v_peptide_data.peptide_rt);
				find_peptide_analysis->v_peptide_withoutmod_m.push_back(itr_v_peptide_data.peptide_m);
				find_peptide_analysis->v_denovo_peptide_averagescore
					= ((find_peptide_analysis->v_denovo_peptide_averagescore * (find_peptide_analysis->replicate_count - 1)) + itr_v_peptide_data.denovo_peptide_data_filtered.localconfidence_average) / find_peptide_analysis->replicate_count;
			}
		}
		return temp_v_peptide_analysis;
	}

	map<string, protein_analysis*> create_v_protein_analysis_map(vector<protein_analysis>& par_v_protein_analysis) {
		map<string, protein_analysis*> temp_v_protein_analysis_map{};
		for (auto& itr_v_peptide_analysis : par_v_protein_analysis) {
			temp_v_protein_analysis_map[itr_v_peptide_analysis.p_protein_data->protein_name] = &itr_v_peptide_analysis;
		}
		return temp_v_protein_analysis_map;
	}

	void create_v_protein_analysis(sample_analysis& par_sample_analysis, const size_t& par_iteration, const bool& par_refined) {
		protein_analysis temp_protein_analysis{};
		size_t temp_key_protein_analysis{};
		vector<protein_analysis> temp_v_protein_analysis{};
		for (const auto& itr_homology_data : par_sample_analysis.v_homology_data) {
			auto find_protein_analysis = std::find_if(temp_v_protein_analysis.begin(), temp_v_protein_analysis.end(),
				[itr_homology_data](const protein_analysis& par_protein_analysis) {
				return par_protein_analysis.p_protein_data->key_protein_data == itr_homology_data.key_blastp_subject_accession;
			});
			if (find_protein_analysis != temp_v_protein_analysis.end()) {
				find_protein_analysis->v_homology_data_combined_by_protein.push_back(itr_homology_data);
			}
			else {
				double temp_protein_score{};
				temp_protein_analysis.p_protein_data = itr_homology_data.p_protein_data;
				if ((par_refined) && (par_iteration > 10) && (itr_homology_data.blastp_parameter_density_conjugated > IgFamily::HOMOLOGY_SCORE_THRESHOLD)) {
					temp_protein_analysis.v_homology_data_combined_by_protein.push_back(itr_homology_data);
				}
				temp_protein_analysis.protein_score = temp_protein_score;
				temp_protein_analysis.key_protein_analysis = temp_key_protein_analysis;
				temp_v_protein_analysis.push_back(temp_protein_analysis);
				++temp_key_protein_analysis;
				temp_protein_analysis.v_homology_data_combined_by_protein.clear();
			}
		}
		denovo_peptide default_denovo_peptide{};
		for (auto& itr_protein_analysis : temp_v_protein_analysis) {
			for (auto& itr_homology_data : itr_protein_analysis.v_homology_data_combined_by_protein) {
				auto& find_peptide_analysis = par_sample_analysis.v_peptide_analysis_map.find(itr_homology_data.blastp_query);
				itr_homology_data.p_peptide_analysis = find_peptide_analysis->second;
				itr_homology_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence = &default_denovo_peptide;
				for (auto& itr_peptide_data : find_peptide_analysis->second->v_peptide_data) {
					++itr_homology_data.denovo_replicate_count;
					if (itr_homology_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average < itr_peptide_data->denovo_peptide_data_filtered.localconfidence_average) {
						itr_homology_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence = &itr_peptide_data->denovo_peptide_data_filtered;
					}
				}
				itr_protein_analysis.protein_score += itr_homology_data.blastp_parameter_score * itr_homology_data.denovo_replicate_count;
			}
		}
		par_sample_analysis.v_protein_analysis = temp_v_protein_analysis;
		create_v_protein_analysis_map(par_sample_analysis.v_protein_analysis);
	}

	inline bool predicate_homology_data_with_spectralcount(const homology_data& i, const homology_data& j) {
		return ((i.blastp_parameter_score * i.denovo_replicate_count * i.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average) > (j.blastp_parameter_score * j.denovo_replicate_count * i.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average));
	}

	inline void sort_v_homology_data_with_spectralcount(vector<homology_data>& par_v_homology_data) {
		std::sort(par_v_homology_data.begin(), par_v_homology_data.end(), predicate_homology_data_with_spectralcount);
	}

	inline bool predicate_protein_analysis(const protein_analysis& i, const protein_analysis& j) {
		return (i.protein_score > j.protein_score);
	}

	inline void sort_v_protein_analysis(vector<protein_analysis>& par_v_homology_data) {
		std::sort(par_v_homology_data.begin(), par_v_homology_data.end(), predicate_protein_analysis);
	}

	void determine_protein_analysis_score_mean(sample_analysis& par_sample_analysis) {
		double temp_protein_analysis_score_mean{};
		size_t count_IG_proteins{};
		for (auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if ((itr_v_protein_analysis.p_protein_data->protein_type == "IG") && (itr_v_protein_analysis.protein_score >= IgFamily::HOMOLOGY_SCORE_THRESHOLD)) {
				++count_IG_proteins;
				temp_protein_analysis_score_mean += itr_v_protein_analysis.protein_score;
			}
		}
		temp_protein_analysis_score_mean /= count_IG_proteins;
		par_sample_analysis.protein_analysis_score_mean = temp_protein_analysis_score_mean;
	}

	void train_homology_analysis_parameter_score(filesystem& par_filesystem, sample_analysis& par_sample_analysis, const size_t& par_select_N_many_gene_families, const bool& par_refined) {
		std::cout << "\n training protein scores...\n\n";
		size_t count_selected_genefamilies{};
		size_t count_iterations{};
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if ((itr_v_protein_analysis.p_protein_data->protein_type == IgFamily::SELECT_TYPE_GENE_FAMILIES) && (itr_v_protein_analysis.protein_score >= IgFamily::HOMOLOGY_SCORE_THRESHOLD)) {
				++count_selected_genefamilies;
			}
		}
		while (count_selected_genefamilies > par_select_N_many_gene_families) {
			determine_protein_analysis_score_mean(par_sample_analysis);
			for (auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis) {
				for (auto& itr_homology_analysis : par_sample_analysis.v_homology_data) {
					if (itr_homology_analysis.blastp_subject_accession == itr_protein_analysis.p_protein_data->protein_name) {
						double score_conjugate{ IgFamily::LOGISTIC_CONJUGATION_FACTOR + (IgFamily::LOGISTIC_CONJUGATION_MIDPOINT / (1 + std::pow(double(2.718), double(-IgFamily::MULTINOMIAL_CONJUGATION_FACTOR) * (itr_protein_analysis.protein_score - par_sample_analysis.protein_analysis_score_mean)))) };
						itr_homology_analysis.blastp_score_transformed_conjugated = std::pow(itr_homology_analysis.blastp_score_transformed_conjugated, score_conjugate);
					}
				}
			}
			++count_iterations;
			std::cout << " ...iteration ";
			std::cout << count_iterations;
			count_selected_genefamilies = 0;
			for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
				if ((itr_v_protein_analysis.p_protein_data->protein_type == IgFamily::SELECT_TYPE_GENE_FAMILIES) && (itr_v_protein_analysis.protein_score > IgFamily::HOMOLOGY_SCORE_THRESHOLD)) {
					++count_selected_genefamilies;
				}
			}
			std::cout << " with ";
			std::cout << count_selected_genefamilies;
			std::cout << " gene families\n";
			fpf_homology_analysis::determine_homology_parameter_density(par_sample_analysis, true);
			IgFamily::HOMOLOGY_SCORE_THRESHOLD = (par_sample_analysis.protein_analysis_score_mean / IgFamily::HOMOLOGY_SCORE_THRESHOLD_FACTOR);
			if (count_selected_genefamilies <= IgFamily::TRAIN_SELECTED_GENE_FAMILY_PRECISION_CONDITION) {
				IgFamily::LOGISTIC_CONJUGATION_FACTOR = 0.9995;
				IgFamily::LOGISTIC_CONJUGATION_MIDPOINT = 0.001;
			}
			if (count_selected_genefamilies <= IgFamily::TRAIN_SELECTED_GENE_FAMILY_PRECISION_CONDITION_2) {
				IgFamily::LOGISTIC_CONJUGATION_FACTOR = 0.99995;
				IgFamily::LOGISTIC_CONJUGATION_MIDPOINT = 0.0001;
			}
			create_v_protein_analysis(par_sample_analysis, count_iterations, par_refined);
		}
	}

	void create_proteinconstruct_from_denovo(sample_analysis& par_sample_analysis) {
		for (auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis) {
			for (size_t i = 0; i < itr_protein_analysis.p_protein_data->protein_protein.length(); ++i) {
				proteinconstruct_aminoacid temp_proteinconstruct_from_denovo{};
				temp_proteinconstruct_from_denovo.aminoacid = '.';
				temp_proteinconstruct_from_denovo.aminoacid_score = 0;
				temp_proteinconstruct_from_denovo.aminoacid_score_conjugated = 0;
				itr_protein_analysis.proteinconstruct_from_denovo.push_back(temp_proteinconstruct_from_denovo);
			}
			sort_v_homology_data_with_spectralcount(itr_protein_analysis.v_homology_data_combined_by_protein);
			vector<homology_data> v_blastp_query_alignment_selected{};
			vector<homology_data> v_blastp_query_alignment_rejected{};
			for (const auto itr_v_homology_data : itr_protein_analysis.v_homology_data_combined_by_protein) {
				const auto find_blastp_query_alignment_rejected = std::find_if(v_blastp_query_alignment_rejected.begin(), v_blastp_query_alignment_rejected.end(),
					[itr_v_homology_data](const homology_data par_s_blastp) {
					return par_s_blastp.query_alignment == itr_v_homology_data.query_alignment;
				});
				if (find_blastp_query_alignment_rejected == v_blastp_query_alignment_rejected.end()) {
					homology_data temp_blastp_query_alignment{};
					temp_blastp_query_alignment.query_alignment = itr_v_homology_data.query_alignment;
					temp_blastp_query_alignment.blastp_score = itr_v_homology_data.blastp_score_transformed;
					temp_blastp_query_alignment.blastp_score_transformed = (itr_v_homology_data.blastp_score_transformed * itr_v_homology_data.blastp_parameter_density_conjugated);
					temp_blastp_query_alignment.blastp_parameter_score = itr_v_homology_data.blastp_parameter_score;
					temp_blastp_query_alignment.blastp_parameter_density_conjugated = itr_v_homology_data.blastp_parameter_density_conjugated;
					temp_blastp_query_alignment.p_peptide_analysis = itr_v_homology_data.p_peptide_analysis;
					temp_blastp_query_alignment.denovo_replicate_count = itr_v_homology_data.denovo_replicate_count;
					for (const auto itr_v_homology_data_2 : itr_protein_analysis.v_homology_data_combined_by_protein) {
						for (auto i = 0; i < itr_v_homology_data.query_alignment.length(); ++i) {
							if ((itr_v_homology_data.query_alignment.at(i) != '.') && (itr_v_homology_data_2.query_alignment.at(i) != '.')) {
								if (itr_v_homology_data.query_alignment.at(i) != itr_v_homology_data_2.query_alignment.at(i)) {
									v_blastp_query_alignment_rejected.push_back(itr_v_homology_data_2);
								}
							}
						}
					}
					if (temp_blastp_query_alignment.query_alignment != "") {
						v_blastp_query_alignment_selected.push_back(temp_blastp_query_alignment);
					}
				}
			}
			for (auto i = 0; i < v_blastp_query_alignment_selected.size(); ++i) {
				bool skip_blastp_query_alignment{};
				for (auto j = 0; j < itr_protein_analysis.proteinconstruct_from_denovo.size(); ++j) {
					if ((itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid != '.')
						&& (v_blastp_query_alignment_selected[i].query_alignment.at(j) != '.')
						&& (itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid != v_blastp_query_alignment_selected[i].query_alignment.at(j))) {
						if (v_blastp_query_alignment_selected[i].blastp_parameter_score < itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid_parameter_score) {
							skip_blastp_query_alignment = true;
							break;
						}
					}
				}
				if (!skip_blastp_query_alignment) {
					for (auto j = 0; j < itr_protein_analysis.proteinconstruct_from_denovo.size(); j) {
						if (v_blastp_query_alignment_selected[i].query_alignment.at(j) == '.') {
							++j;
						}
						else {
							for (const auto& itr_v_denovo_aminoacid : v_blastp_query_alignment_selected[i].p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->v_denovo_aminoacid) {
								if (j == itr_protein_analysis.proteinconstruct_from_denovo.size()) {
									break;
								}
								else {
									if ((v_blastp_query_alignment_selected[i].blastp_parameter_score >= itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid_parameter_score)
										&& (v_blastp_query_alignment_selected[i].blastp_parameter_score >= IgFamily::REPORT_QUERY_ALIGNMENT_PARAMETER_SCORE_THRESHOLD)
										&& ((v_blastp_query_alignment_selected[i].blastp_parameter_score * v_blastp_query_alignment_selected[i].denovo_replicate_count) >= IgFamily::REPORT_QUERY_ALIGNMENT_TOTALSCORE_OUTPUT_THRESHOLD)) {
										itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid = v_blastp_query_alignment_selected[i].query_alignment.at(j);
										itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid_localconfidence = itr_v_denovo_aminoacid.aminoacid_localconfidence;
										itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid_score = v_blastp_query_alignment_selected[i].blastp_score;
										itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid_score_conjugated = v_blastp_query_alignment_selected[i].blastp_score_transformed;
										itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid_parameter_score = v_blastp_query_alignment_selected[i].blastp_parameter_score;
										itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid_parameter_density_conjugated = v_blastp_query_alignment_selected[i].blastp_parameter_density_conjugated;
									}
									++j;
								}
							}
						}
					}
				}
			}
		}
	}

	void determine_sequence_coverage(sample_analysis& par_sample_analysis) {
		for (auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis) {
			size_t temp_sequencetrue{};
			for (const auto& itr_v_proteinconstruct : itr_protein_analysis.proteinconstruct_from_denovo) {
				if (itr_v_proteinconstruct.aminoacid != '.') {
					++temp_sequencetrue;
				}
			}
			itr_protein_analysis.proteinconstruct_sequencecoverage = double(100 * (double(temp_sequencetrue) / itr_protein_analysis.proteinconstruct_from_denovo.size()));
		}
	}

	void select_protein_analysis_by_score(sample_analysis& par_sample_analysis) {
		vector<protein_analysis> temp_v_protein_analysis_selected_by_polymorphism{};
		for (const auto itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if (itr_v_protein_analysis.p_protein_data->protein_type == "IG") {
				string protein_name_polymorphism_reduced{};
				bool switch_protein_name_polymorphism = bool();
				for (const auto& itr_protein_name : itr_v_protein_analysis.p_protein_data->protein_name) {
					if (itr_protein_name == '*') {
						switch_protein_name_polymorphism = true;
					}
					if (!switch_protein_name_polymorphism) {
						protein_name_polymorphism_reduced += itr_protein_name;
					}
				}
				auto& find_protein_analysis = std::find_if(temp_v_protein_analysis_selected_by_polymorphism.begin(), temp_v_protein_analysis_selected_by_polymorphism.end(),
					[protein_name_polymorphism_reduced](const protein_analysis par_protein_analysis) {
					return par_protein_analysis.p_protein_data->protein_name == protein_name_polymorphism_reduced;
				});
				if (find_protein_analysis != temp_v_protein_analysis_selected_by_polymorphism.end()) {
					if (itr_v_protein_analysis.protein_score > find_protein_analysis->protein_score) {
						*find_protein_analysis = itr_v_protein_analysis;
					}
				}
				else {
					temp_v_protein_analysis_selected_by_polymorphism.push_back(itr_v_protein_analysis);
					temp_v_protein_analysis_selected_by_polymorphism.back().p_protein_data->protein_name = protein_name_polymorphism_reduced;
				}
			}
			else {
				temp_v_protein_analysis_selected_by_polymorphism.push_back(itr_v_protein_analysis);
			}
		}
		par_sample_analysis.v_protein_analysis_selected_by_polymorphism = temp_v_protein_analysis_selected_by_polymorphism;
	}

	void determine_protein_score_density(sample_analysis& par_sample_analysis) {
		double temp_protein_score_sum{};
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if (itr_v_protein_analysis.p_protein_data->protein_type == "IG") {
				temp_protein_score_sum += itr_v_protein_analysis.protein_score;
			}
		}
		for (auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if (itr_v_protein_analysis.p_protein_data->protein_type == "IG") {
				itr_v_protein_analysis.protein_density = itr_v_protein_analysis.protein_score / temp_protein_score_sum;
			}
		}
	}
}

#endif