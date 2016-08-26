// * * fpf_data_analysis.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_DATA_ANALYSIS
#define	FPF_DATA_ANALYSIS

#include <cstdlib>						// provides - size_t
#include <string>						// provides - std::string
#include <vector>						// provides - std::vector
#include <iomanip>						// provides - std::setprecision
#include <iostream>						// provides - std::get
#include <utility>						// provides - std::pair

#include "fpf_filesystem.h"
#include "fpf_data.h"


namespace fpf_data_analysis {

	using std::string;
	using std::vector;

	typedef fpf_filesystem::filesystem filesystem;
	typedef fpf_data::peptide_data peptide_data;
	typedef fpf_data::peptide_analysis peptide_analysis;
	typedef fpf_data::denovo_peptide denovo_peptide;
	typedef fpf_data::FASTA_category FASTA_category;
	typedef fpf_data::blastp_data blastp_data;
	typedef fpf_data::v_proteinconstruct_from_denovo v_proteinconstruct_from_denovo;
	typedef fpf_data::category_analysis category_analysis;
	typedef fpf_data::multinomial multinomial;

	void create_category_analysis(filesystem& par_filesystem) {
		category_analysis temp_category_analysis{};
		vector<category_analysis> temp_v_category_analysis{};
		for (const auto itr_blastp_data : par_filesystem.v_blastp_data) {
			if (itr_blastp_data.blastp_evalue_transformed > BLASTP_EVALUETRANSFORMED_THRESHOLD) {
				auto find_category_analysis = std::find_if(temp_v_category_analysis.begin(), temp_v_category_analysis.end(),
					[itr_blastp_data](const category_analysis par_category_analysis) {
					return par_category_analysis.p_FASTA_category->category_name == itr_blastp_data.blastp_subject_accession;
				});
				if (find_category_analysis != temp_v_category_analysis.end()) {
					find_category_analysis->v_blastp_data_combined_by_category.push_back(itr_blastp_data);
				}
				else {
					double temp_category_score{};
					temp_category_analysis.p_FASTA_category = itr_blastp_data.p_FASTA_category;
					temp_category_analysis.v_blastp_data_combined_by_category.push_back(itr_blastp_data);
					temp_category_analysis.category_score = temp_category_score;
					temp_v_category_analysis.push_back(temp_category_analysis);
					temp_category_analysis.v_blastp_data_combined_by_category.clear();
				}
			}
		}
		denovo_peptide default_denovo_peptide{};
		for (auto& itr_category_analysis : temp_v_category_analysis) {
			for (auto& itr_blastp_data : itr_category_analysis.v_blastp_data_combined_by_category) {
				auto& find_peptide_analysis = std::find_if(par_filesystem.v_peptide_analysis.begin(), par_filesystem.v_peptide_analysis.end(),
					[itr_blastp_data](const peptide_analysis par_peptide_analysis) {
					return par_peptide_analysis.peptide_filtered == itr_blastp_data.blastp_query;
				});
				itr_blastp_data.p_peptide_analysis = &(*find_peptide_analysis);
				itr_blastp_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence = &default_denovo_peptide;
				for (auto& itr_peptide_data : find_peptide_analysis->v_peptide_data) {
					++itr_blastp_data.denovo_replicate_count;
					if (itr_blastp_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average < itr_peptide_data->denovo_peptide_data.localconfidence_average) {
						itr_blastp_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence = &itr_peptide_data->denovo_peptide_data;
					}
				}
				itr_category_analysis.category_score += (itr_blastp_data.blastp_parameter_score * itr_blastp_data.denovo_replicate_count);
			}
		}
		par_filesystem.v_category_analysis = temp_v_category_analysis;
	}

	inline bool predicate_blastp_data_with_spectralcount(const blastp_data& i, const blastp_data& j) {
		return ((i.blastp_parameter_score * i.denovo_replicate_count * i.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average) > (j.blastp_parameter_score * j.denovo_replicate_count * i.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average));
	}

	inline void sort_v_blastp_data_with_spectralcount(vector<blastp_data>& par_v_blastp_data) {
		std::sort(par_v_blastp_data.begin(), par_v_blastp_data.end(), predicate_blastp_data_with_spectralcount);
	}

	inline bool predicate_category_analysis(const category_analysis& i, const category_analysis& j) {
		return (i.category_score > j.category_score);
	}

	inline void sort_v_category_analysis(vector<category_analysis>& par_v_blastp_data) {
		std::sort(par_v_blastp_data.begin(), par_v_blastp_data.end(), predicate_category_analysis);
	}

	void create_proteinconstruct_from_denovo(filesystem& par_filesystem) {
		for (auto& itr_category_analysis : par_filesystem.v_category_analysis) {
			for (size_t i = 0; i < itr_category_analysis.p_FASTA_category->category_protein.length(); ++i) {
				v_proteinconstruct_from_denovo temp_proteinconstruct_from_denovo{};
				temp_proteinconstruct_from_denovo.aminoacid = '.';
				temp_proteinconstruct_from_denovo.aminoacid_evalue_transformed = 0;
				itr_category_analysis.v_proteinconstruct_from_denovo.push_back(temp_proteinconstruct_from_denovo);
			}
			sort_v_blastp_data_with_spectralcount(itr_category_analysis.v_blastp_data_combined_by_category);
			vector<blastp_data> v_blastp_query_alignment_selected{};
			vector<blastp_data> v_blastp_query_alignment_rejected{};
			for (const auto itr_v_blastp_data : itr_category_analysis.v_blastp_data_combined_by_category) {
				const auto find_blastp_query_alignment_rejected = std::find_if(v_blastp_query_alignment_rejected.begin(), v_blastp_query_alignment_rejected.end(),
					[itr_v_blastp_data](const blastp_data par_s_blastp) {
					return par_s_blastp.query_alignment == itr_v_blastp_data.query_alignment;
				});
				if (find_blastp_query_alignment_rejected == v_blastp_query_alignment_rejected.end()) {
					blastp_data temp_blastp_query_alignment{};
					temp_blastp_query_alignment.query_alignment = itr_v_blastp_data.query_alignment;
					temp_blastp_query_alignment.blastp_evalue_transformed = itr_v_blastp_data.blastp_evalue_transformed;
					temp_blastp_query_alignment.blastp_parameter_score = itr_v_blastp_data.blastp_parameter_score;
					temp_blastp_query_alignment.p_peptide_analysis = itr_v_blastp_data.p_peptide_analysis;
					temp_blastp_query_alignment.denovo_replicate_count = itr_v_blastp_data.denovo_replicate_count;
					for (const auto itr_v_blastp_data_2 : itr_category_analysis.v_blastp_data_combined_by_category) {
						for (auto i = 0; i < itr_v_blastp_data.query_alignment.length(); ++i) {
							if ((itr_v_blastp_data.query_alignment.at(i) != '.') && (itr_v_blastp_data_2.query_alignment.at(i) != '.')) {
								if (itr_v_blastp_data.query_alignment.at(i) != itr_v_blastp_data_2.query_alignment.at(i)) {
									v_blastp_query_alignment_rejected.push_back(itr_v_blastp_data_2);
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
				for (auto j = 0; j < itr_category_analysis.v_proteinconstruct_from_denovo.size(); ++j) {
					if ((itr_category_analysis.v_proteinconstruct_from_denovo[j].aminoacid != '.')
						&& (v_blastp_query_alignment_selected[i].query_alignment.at(j) != '.')
						&& (itr_category_analysis.v_proteinconstruct_from_denovo[j].aminoacid != v_blastp_query_alignment_selected[i].query_alignment.at(j))) {
						if (v_blastp_query_alignment_selected[i].blastp_parameter_score < itr_category_analysis.v_proteinconstruct_from_denovo[j].aminoacid_parameter_score) {
							skip_blastp_query_alignment = true;
							break;
						}
					}
				}
				if (!skip_blastp_query_alignment) {
					for (auto j = 0; j < itr_category_analysis.v_proteinconstruct_from_denovo.size(); j) {
						if (v_blastp_query_alignment_selected[i].query_alignment.at(j) == '.') {
							++j;
						}
						else {
							for (const auto& itr_v_denovo_aminoacid : v_blastp_query_alignment_selected[i].p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->v_denovo_aminoacid) {
								if (j == itr_category_analysis.v_proteinconstruct_from_denovo.size()) {
									break;
								}
								else {
									if (v_blastp_query_alignment_selected[i].blastp_parameter_score > itr_category_analysis.v_proteinconstruct_from_denovo[j].aminoacid_parameter_score) {
										if (v_blastp_query_alignment_selected[i].blastp_evalue_transformed > BLASTP_EVALUETRANSFORMED_THRESHOLD) {
											itr_category_analysis.v_proteinconstruct_from_denovo[j].aminoacid = v_blastp_query_alignment_selected[i].query_alignment.at(j);
											itr_category_analysis.v_proteinconstruct_from_denovo[j].aminoacid_localconfidence = itr_v_denovo_aminoacid.aminoacid_localconfidence;
											itr_category_analysis.v_proteinconstruct_from_denovo[j].aminoacid_evalue_transformed = v_blastp_query_alignment_selected[i].blastp_evalue_transformed;
											itr_category_analysis.v_proteinconstruct_from_denovo[j].aminoacid_parameter_score = v_blastp_query_alignment_selected[i].blastp_parameter_score;
										}
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

	void determine_sequence_coverage(filesystem& par_filesystem) {
		for (auto& itr_category_analysis : par_filesystem.v_category_analysis) {
			size_t temp_sequencetrue{};
			for (const auto& itr_v_proteinconstruct : itr_category_analysis.v_proteinconstruct_from_denovo) {
				if (itr_v_proteinconstruct.aminoacid != '.') {
					++temp_sequencetrue;
				}
			}
			itr_category_analysis.proteinconstruct_sequencecoverage = double(100 * (double(temp_sequencetrue) / itr_category_analysis.v_proteinconstruct_from_denovo.size()));
		}
	}

	void select_category_analysis_by_score(filesystem& par_filesystem) {
		vector<category_analysis> temp_v_category_analysis_selected_by_polymorphism{};
		for (const auto itr_v_category_analysis : par_filesystem.v_category_analysis) {
			if (itr_v_category_analysis.p_FASTA_category->category_type == "IG") {
				string category_name_polymorphism_reduced{};
				bool switch_category_name_polymorphism = bool();
				for (const auto& itr_category_name : itr_v_category_analysis.p_FASTA_category->category_name) {
					if (itr_category_name == '*') {
						switch_category_name_polymorphism = true;
					}
					if (!switch_category_name_polymorphism) {
						category_name_polymorphism_reduced += itr_category_name;
					}
				}
				auto& find_category_analysis = std::find_if(temp_v_category_analysis_selected_by_polymorphism.begin(), temp_v_category_analysis_selected_by_polymorphism.end(),
					[category_name_polymorphism_reduced](const category_analysis par_category_analysis) {
					return par_category_analysis.p_FASTA_category->category_name == category_name_polymorphism_reduced;
				});
				if (find_category_analysis != temp_v_category_analysis_selected_by_polymorphism.end()) {
					if (itr_v_category_analysis.category_score > find_category_analysis->category_score) {
						*find_category_analysis = itr_v_category_analysis;
					}
				}
				else {
					temp_v_category_analysis_selected_by_polymorphism.push_back(itr_v_category_analysis);
					temp_v_category_analysis_selected_by_polymorphism.back().p_FASTA_category->category_name = category_name_polymorphism_reduced;
				}
			}
			else {
				temp_v_category_analysis_selected_by_polymorphism.push_back(itr_v_category_analysis);
			}
		}
		par_filesystem.v_category_analysis_selected_by_polymorphism = temp_v_category_analysis_selected_by_polymorphism;
	}
}

#endif