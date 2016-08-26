// * * fpf_data_analysis.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_DATA_ANALYSIS
#define	FPF_DATA_ANALYSIS

#include <cstdlib>						// provides - size_t
#include <iomanip>						// provides - std::setprecision
#include <iostream>						// provides - std::get
#include <string>						// provides - std::string
#include <utility>						// provides - std::pair
#include <vector>						// provides - std::vector

#include "fpf_data.h"
#include "fpf_filesystem.h"


namespace fpf_data_analysis {

	using std::string;
	using std::vector;

	typedef fpf_data::blastp_data blastp_data;
	typedef fpf_data::denovo_peptide denovo_peptide;
	typedef fpf_data::protein_data protein_data;
	typedef fpf_data::peptide_analysis peptide_analysis;
	typedef fpf_data::peptide_data peptide_data;
	typedef fpf_data::protein_analysis protein_analysis;
	typedef fpf_data::multinomial multinomial;
	typedef fpf_data::proteinconstruct_aminoacid proteinconstruct_aminoacid;
	typedef fpf_filesystem::filesystem filesystem;
	typedef fpf_filesystem::sample_analysis sample_analysis;
	typedef fpf_parse::csv_data csv_data;
	typedef fpf_parse::FASTA_data FASTA_data;

	void create_protein_analysis(sample_analysis& par_sample_analysis) {
		protein_analysis temp_protein_analysis{};
		vector<protein_analysis> temp_v_protein_analysis{};
		for (const auto itr_blastp_data : par_sample_analysis.v_blastp_data) {
			if (itr_blastp_data.blastp_evalue_transformed > BLASTP_EVALUETRANSFORMED_THRESHOLD) {
				auto find_protein_analysis = std::find_if(temp_v_protein_analysis.begin(), temp_v_protein_analysis.end(),
					[itr_blastp_data](const protein_analysis par_protein_analysis) {
					return par_protein_analysis.p_protein_data->protein_name == itr_blastp_data.blastp_subject_accession;
				});
				if (find_protein_analysis != temp_v_protein_analysis.end()) {
					find_protein_analysis->v_blastp_data_combined_by_protein.push_back(itr_blastp_data);
				}
				else {
					double temp_protein_score{};
					temp_protein_analysis.p_protein_data = itr_blastp_data.p_protein_data;
					temp_protein_analysis.v_blastp_data_combined_by_protein.push_back(itr_blastp_data);
					temp_protein_analysis.protein_score = temp_protein_score;
					temp_v_protein_analysis.push_back(temp_protein_analysis);
					temp_protein_analysis.v_blastp_data_combined_by_protein.clear();
				}
			}
		}
		denovo_peptide default_denovo_peptide{};
		for (auto& itr_protein_analysis : temp_v_protein_analysis) {
			for (auto& itr_blastp_data : itr_protein_analysis.v_blastp_data_combined_by_protein) {
				auto& find_peptide_analysis = std::find_if(par_sample_analysis.v_peptide_analysis.begin(), par_sample_analysis.v_peptide_analysis.end(),
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
				itr_protein_analysis.protein_score += (itr_blastp_data.blastp_parameter_score * itr_blastp_data.denovo_replicate_count);
			}
		}
		par_sample_analysis.v_protein_analysis = temp_v_protein_analysis;
	}

	inline bool predicate_blastp_data_with_spectralcount(const blastp_data& i, const blastp_data& j) {
		return ((i.blastp_parameter_score * i.denovo_replicate_count * i.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average) > (j.blastp_parameter_score * j.denovo_replicate_count * i.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average));
	}

	inline void sort_v_blastp_data_with_spectralcount(vector<blastp_data>& par_v_blastp_data) {
		std::sort(par_v_blastp_data.begin(), par_v_blastp_data.end(), predicate_blastp_data_with_spectralcount);
	}

	inline bool predicate_protein_analysis(const protein_analysis& i, const protein_analysis& j) {
		return (i.protein_score > j.protein_score);
	}

	inline void sort_v_protein_analysis(vector<protein_analysis>& par_v_blastp_data) {
		std::sort(par_v_blastp_data.begin(), par_v_blastp_data.end(), predicate_protein_analysis);
	}

	void create_proteinconstruct_from_denovo(sample_analysis& par_sample_analysis) {
		for (auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis) {
			for (size_t i = 0; i < itr_protein_analysis.p_protein_data->protein_protein.length(); ++i) {
				proteinconstruct_aminoacid temp_proteinconstruct_from_denovo{};
				temp_proteinconstruct_from_denovo.aminoacid = '.';
				temp_proteinconstruct_from_denovo.aminoacid_evalue_transformed = 0;
				itr_protein_analysis.proteinconstruct_from_denovo.push_back(temp_proteinconstruct_from_denovo);
			}
			sort_v_blastp_data_with_spectralcount(itr_protein_analysis.v_blastp_data_combined_by_protein);
			vector<blastp_data> v_blastp_query_alignment_selected{};
			vector<blastp_data> v_blastp_query_alignment_rejected{};
			for (const auto itr_v_blastp_data : itr_protein_analysis.v_blastp_data_combined_by_protein) {
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
					for (const auto itr_v_blastp_data_2 : itr_protein_analysis.v_blastp_data_combined_by_protein) {
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
									if (v_blastp_query_alignment_selected[i].blastp_parameter_score > itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid_parameter_score) {
										if (v_blastp_query_alignment_selected[i].blastp_evalue_transformed > BLASTP_EVALUETRANSFORMED_THRESHOLD) {
											itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid = v_blastp_query_alignment_selected[i].query_alignment.at(j);
											itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid_localconfidence = itr_v_denovo_aminoacid.aminoacid_localconfidence;
											itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid_evalue_transformed = v_blastp_query_alignment_selected[i].blastp_evalue_transformed;
											itr_protein_analysis.proteinconstruct_from_denovo[j].aminoacid_parameter_score = v_blastp_query_alignment_selected[i].blastp_parameter_score;
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
}

#endif