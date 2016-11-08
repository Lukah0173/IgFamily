// * * fpf_data_analysis.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_DATA_ANALYSIS
#define	FPF_DATA_ANALYSIS

#include <algorithm>
#include <cstdlib>
#include <functional>
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
	using std::multimap;
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

	//map<string, protein_analysis*> create_v_protein_analysis_map(vector<protein_analysis>& par_v_protein_analysis) {
	//	map<string, protein_analysis*> temp_v_protein_analysis_map{};
	//	for (auto& itr_v_peptide_analysis : par_v_protein_analysis) {
	//		temp_v_protein_analysis_map[itr_v_peptide_analysis.p_protein_data->protein_name] = &itr_v_peptide_analysis;
	//	}
	//	return temp_v_protein_analysis_map;
	//}

	void create_v_protein_analysis(sample_analysis& par_sample_analysis, const size_t& par_iteration, const bool& par_refined, const bool& par_filter_by_score) {
		vector<protein_analysis> temp_v_protein_analysis{};
		size_t temp_key_protein_analysis{};
		for (const auto& itr_homology_data : par_sample_analysis.v_homology_data)
		{
			protein_analysis temp_protein_analysis{};
			auto find_protein_analysis = std::find_if(temp_v_protein_analysis.begin(), temp_v_protein_analysis.end(),
				[itr_homology_data](const protein_analysis& par_protein_analysis)
			{
				return par_protein_analysis.p_protein_data->key_protein_data == itr_homology_data.key_blastp_subject_accession;
			});
			if (find_protein_analysis != temp_v_protein_analysis.end())
			{
				find_protein_analysis->v_homology_data_combined_by_protein.push_back(itr_homology_data);
			}
			else
			{
				double temp_protein_score{};
				temp_protein_analysis.p_protein_data = itr_homology_data.p_protein_data;
				temp_protein_analysis.v_homology_data_combined_by_protein.push_back(itr_homology_data);
				temp_protein_analysis.protein_score = temp_protein_score;
				temp_protein_analysis.key_protein_analysis = temp_key_protein_analysis;
				temp_v_protein_analysis.push_back(temp_protein_analysis);
				++temp_key_protein_analysis;
				temp_protein_analysis.v_homology_data_combined_by_protein.clear();
			}
		}
		for (auto& itr_protein_analysis : temp_v_protein_analysis)
		{
			for (auto& itr_homology_data : itr_protein_analysis.v_homology_data_combined_by_protein) 
			{					
				itr_protein_analysis.protein_score += itr_homology_data.score;
				itr_protein_analysis.protein_effective_spectral_count += (itr_homology_data.blastp_homology_density_conjugated * itr_homology_data.denovo_replicate_count);
			}
			for (auto& itr_homology_data : itr_protein_analysis.v_homology_data_combined_by_protein)
			{
				if (itr_protein_analysis.protein_score != 0)
				{
					itr_homology_data.score_density = (itr_homology_data.score / itr_protein_analysis.protein_score);
				}
			}
		}
		par_sample_analysis.v_protein_analysis = temp_v_protein_analysis;
	}

	inline bool predicate_homology_data_combined_by_protein(const homology_data& i, const homology_data& j) {
		return ((i.score) > (j.score));
	}

	inline void sort_v_homology_data_combined_by_protein(vector<homology_data>& par_v_homology_data) {
		std::sort(par_v_homology_data.begin(), par_v_homology_data.end(), predicate_homology_data_combined_by_protein);
	}

	inline bool predicate_protein_analysis(const protein_analysis& i, const protein_analysis& j) {
		return (i.protein_score > j.protein_score);
	}

	inline void sort_v_protein_analysis(vector<protein_analysis>& par_v_homology_data) {
		std::sort(par_v_homology_data.begin(), par_v_homology_data.end(), predicate_protein_analysis);
	}

	void determine_protein_analysis_score_mean(sample_analysis& par_sample_analysis)
	{
		double temp_protein_analysis_score_max{};
		double temp_protein_analysis_score_mean{};
		size_t count_IG_proteins{};
		for (auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis)
		{
			if (itr_protein_analysis.protein_score > temp_protein_analysis_score_max)
			{
				temp_protein_analysis_score_max = itr_protein_analysis.protein_score;
			}
			if (itr_protein_analysis.p_protein_data->protein_type == "IGV")
			{
				++count_IG_proteins;
				temp_protein_analysis_score_mean += itr_protein_analysis.protein_score;
			}
		}
		temp_protein_analysis_score_mean /= count_IG_proteins;
		par_sample_analysis.protein_analysis_score_max = temp_protein_analysis_score_max;
		par_sample_analysis.protein_analysis_score_mean = temp_protein_analysis_score_mean;
	}

	size_t determine_countClusterProportion(sample_analysis& par_sample_analysis, double par_ClusterProportion)
	{
		vector<double> GenefamiliesScores{};
		double sum_GenefamilyScore{};
		for (const auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis)
		{
			if (itr_protein_analysis.p_protein_data->protein_type == "IGV")
			{
				GenefamiliesScores.push_back(itr_protein_analysis.protein_score);
				sum_GenefamilyScore += itr_protein_analysis.protein_score;
			}
		}
		std::sort(GenefamiliesScores.begin(), GenefamiliesScores.end(), std::greater<double>());
		double proportion_Cluster{};
		size_t inefficient_countCluster{}; // ugh
		for (const auto& itr_GenefamiliesScores : GenefamiliesScores)
		{
			proportion_Cluster += itr_GenefamiliesScores;
			++inefficient_countCluster;
			if ((proportion_Cluster / sum_GenefamilyScore) >= par_ClusterProportion)
			{
				return inefficient_countCluster;
			}
		}
		return 0;
	}

	void conjugate_homology(filesystem& par_filesystem, sample_analysis& par_sample_analysis, const size_t& par_select_N_many_gene_families, const bool& par_refined) {
		std::cout << "\n training protein scores...\n\n";
		size_t count_selected_genefamilies{};
		size_t count_iterations{};
		count_selected_genefamilies = determine_countClusterProportion(par_sample_analysis, IgFamily::CLUSTER_PROPORTION_THRESHOLD);
		while ((count_selected_genefamilies > par_select_N_many_gene_families) && (count_iterations <= 5000)) {
			determine_protein_analysis_score_mean(par_sample_analysis);
			for (auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis) {
				double score_conjugate{ std::pow((itr_protein_analysis.protein_score / par_sample_analysis.protein_analysis_score_max), IgFamily::PARAMETER_LOGISTIC_CONJUGATION_FACTOR) };
				for (auto& itr_homology_data : par_sample_analysis.v_homology_data) {
					if (itr_homology_data.blastp_subject_accession == itr_protein_analysis.p_protein_data->protein_name) {
						itr_homology_data.blastp_homology_transformed_conjugated = (itr_homology_data.blastp_homology_transformed * score_conjugate);
					}
				}
			}
			fpf_homology_analysis::determine_HomologyDataParameters(par_sample_analysis, true);
			IgFamily::PARAMETER_LOGISTIC_CONJUGATION_FACTOR = (IgFamily::PARAMETER_LOGISTIC_CONJUGATION_FACTOR + (IgFamily::PARAMETER_LOGISTIC_ITERATION_FACTOR * std::pow(1.01, (count_selected_genefamilies - par_select_N_many_gene_families))));
			create_v_protein_analysis(par_sample_analysis, count_iterations, par_refined, true);
			count_selected_genefamilies = {};
			count_selected_genefamilies = determine_countClusterProportion(par_sample_analysis, IgFamily::CLUSTER_PROPORTION_THRESHOLD);
			++count_iterations;
			if (count_iterations == 1) {
				std::cout << " ...iteration ";
				std::cout << count_iterations;
				std::cout << " with ";
				std::cout << count_selected_genefamilies;
				std::cout << " gene families \n";
			}
			if (count_iterations % 50 == 0) {
				std::cout << " ...iteration ";
				std::cout << count_iterations;
				std::cout << " with ";
				std::cout << count_selected_genefamilies;
				std::cout << " gene families \n";
			}
			if ((count_selected_genefamilies <= par_select_N_many_gene_families) && (count_iterations % 50 != 0)) {
				std::cout << " ...iteration ";
				std::cout << count_iterations;
				std::cout << " with ";
				std::cout << count_selected_genefamilies;
				std::cout << " gene families \n";
			}
		}
	}

	void create_ProteinConstruct(sample_analysis& par_sample_analysis) {
		static homology_data temp_homology_data{};
		for (auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis) {
			for (size_t i = 0; i < itr_protein_analysis.p_protein_data->protein_protein.length(); ++i) {
				proteinconstruct_aminoacid temp_proteinconstruct_from_denovo{};
				temp_proteinconstruct_from_denovo.aminoacid = '.';
				temp_proteinconstruct_from_denovo.proteinconstruct_homology_data = temp_homology_data;
				itr_protein_analysis.proteinconstruct.push_back(temp_proteinconstruct_from_denovo);
			}
			sort_v_homology_data_combined_by_protein(itr_protein_analysis.v_homology_data_combined_by_protein);
			vector<homology_data> v_blastp_query_alignment_selected{};
			vector<homology_data> v_blastp_query_alignment_rejected{};
			for (const auto& itr_homology_data : itr_protein_analysis.v_homology_data_combined_by_protein) {
				const auto find_blastp_query_alignment_rejected = std::find_if(v_blastp_query_alignment_rejected.begin(), v_blastp_query_alignment_rejected.end(),
					[itr_homology_data](const auto& par_homology_data) {
					return par_homology_data.alignment == itr_homology_data.alignment;
				});
				if (find_blastp_query_alignment_rejected == v_blastp_query_alignment_rejected.end()) {
					for (const auto& itr_v_homology_data_2 : itr_protein_analysis.v_homology_data_combined_by_protein) {
						for (auto i = 0; i < itr_homology_data.alignment.length(); ++i) {
							if ((itr_homology_data.alignment.at(i) != '.') && (itr_v_homology_data_2.alignment.at(i) != '.')) {
								if (itr_homology_data.alignment.at(i) != itr_v_homology_data_2.alignment.at(i)) {
									v_blastp_query_alignment_rejected.push_back(itr_v_homology_data_2);
								}
							}
						}
					}
					if (itr_homology_data.alignment != "") {
						v_blastp_query_alignment_selected.push_back(itr_homology_data);
					}
				}
			}
			for (auto i = 0; i < v_blastp_query_alignment_selected.size(); ++i) {
				bool skip_blastp_query_alignment{};
				for (auto j = 0; j < itr_protein_analysis.proteinconstruct.size(); ++j) {
					if ((itr_protein_analysis.proteinconstruct[j].aminoacid != '.')
						&& (v_blastp_query_alignment_selected[i].alignment.at(j) != '.')
						&& (itr_protein_analysis.proteinconstruct[j].aminoacid != v_blastp_query_alignment_selected[i].alignment.at(j))) {
						if ((v_blastp_query_alignment_selected[i].score) 
							< itr_protein_analysis.proteinconstruct[j].proteinconstruct_homology_data.score) {
							skip_blastp_query_alignment = true;
							break;
						}
					}
				}
				if (!skip_blastp_query_alignment) {
					for (auto j = 0; j < itr_protein_analysis.proteinconstruct.size(); j) {
						if (v_blastp_query_alignment_selected[i].alignment.at(j) == '.') {
							++j;
						}
						else {
							for (const auto& itr_v_denovo_aminoacid : v_blastp_query_alignment_selected[i].p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->v_denovo_aminoacid) {
								if (j == itr_protein_analysis.proteinconstruct.size()) {
									break;
								}
								else {
									if ((v_blastp_query_alignment_selected[i].score >= itr_protein_analysis.proteinconstruct[j].proteinconstruct_homology_data.score)
										&& (v_blastp_query_alignment_selected[i].blastp_homology_density_conjugated >= IgFamily::PROTEIN_CONSTRUCT_PARAMETER_DENSITY_CONJUGATED_THRESHOLD)
										&& (v_blastp_query_alignment_selected[i].score >= IgFamily::PROTEIN_CONSTRUCT_SCORE_THRESHOLD)) {
										itr_protein_analysis.proteinconstruct[j].aminoacid = v_blastp_query_alignment_selected[i].alignment.at(j);
										itr_protein_analysis.proteinconstruct[j].aminoacid_localconfidence = itr_v_denovo_aminoacid.aminoacid_localconfidence;
										itr_protein_analysis.proteinconstruct[j].proteinconstruct_homology_data = v_blastp_query_alignment_selected[i];
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

	void determine_SequenceCoverage(sample_analysis& par_sample_analysis) {
		for (auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis) {
			size_t temp_sequencetrue{};
			for (const auto& itr_v_proteinconstruct : itr_protein_analysis.proteinconstruct) {
				if (itr_v_proteinconstruct.aminoacid != '.') {
					++temp_sequencetrue;
				}
			}
			itr_protein_analysis.proteinconstruct_sequencecoverage = double(100 * (double(temp_sequencetrue) / itr_protein_analysis.proteinconstruct.size()));
		}
	}

	void select_protein_analysis_by_score(sample_analysis& par_sample_analysis) {
		vector<protein_analysis> temp_v_protein_analysis_selected_by_polymorphism{};
		for (const auto itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if (itr_v_protein_analysis.p_protein_data->protein_type == "IGV") {
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
					temp_v_protein_analysis_selected_by_polymorphism.back().p_protein_data->protein_nonreduced_name = temp_v_protein_analysis_selected_by_polymorphism.back().p_protein_data->protein_name;
					temp_v_protein_analysis_selected_by_polymorphism.back().p_protein_data->protein_name = protein_name_polymorphism_reduced;
				}
			}
			else {
				temp_v_protein_analysis_selected_by_polymorphism.push_back(itr_v_protein_analysis);
			}
		}
		par_sample_analysis.v_protein_analysis_selected_by_polymorphism = temp_v_protein_analysis_selected_by_polymorphism;
	}

	void determine_ProteinScoreDensity(sample_analysis& par_sample_analysis) {
		double temp_protein_score_sum{};
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if (itr_v_protein_analysis.p_protein_data->protein_type == "IGV") {
				temp_protein_score_sum += itr_v_protein_analysis.protein_score;
			}
		}
		for (auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if (itr_v_protein_analysis.p_protein_data->protein_type == "IGV") {
				itr_v_protein_analysis.protein_density = itr_v_protein_analysis.protein_score / temp_protein_score_sum;
			}
		}
	}
}

#endif