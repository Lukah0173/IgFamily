// * * fpf_multinomial.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_MULTINOMIAL
#define	FPF_MULTINOMIAL

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "fpf_data.h"
#include "fpf_filesystem.h"


namespace fpf_multinomial {

	using std::string;
	using std::vector;

	using fpf_data::multinomial_frequency_type;
	using fpf_data::protein_data;
	using fpf_data::peptide_analysis;
	using fpf_filesystem::filesystem;
	using fpf_filesystem::sample_analysis;

	void create_MultinomialData(sample_analysis& par_sample_analysis) {
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis_with_selected_polymorphism) {
			par_sample_analysis.multinomial_data.v_p_protein_data.push_back(itr_v_protein_analysis.p_protein_data);
		}
		for (const auto& itr_v_peptide_analysis_map : par_sample_analysis.v_peptide_analysis_map) {
			par_sample_analysis.multinomial_data.v_p_peptide_analysis.push_back(itr_v_peptide_analysis_map.second);
		}
		par_sample_analysis.multinomial_data.v2_frequency = vector<vector<double>>(par_sample_analysis.multinomial_data.v_p_peptide_analysis.size(), vector<double>(par_sample_analysis.multinomial_data.v_p_protein_data.size(), 0));
		par_sample_analysis.multinomial_data.v2_density = vector<vector<double>>(par_sample_analysis.multinomial_data.v_p_peptide_analysis.size(), vector<double>(par_sample_analysis.multinomial_data.v_p_protein_data.size(), 0));
		par_sample_analysis.multinomial_data.v_frequency_marginal_sum = vector<double>(par_sample_analysis.multinomial_data.v_p_peptide_analysis.size(), 0);

		for (const auto itr_v_homology_data : par_sample_analysis.v_homology_data) {	
			const auto& find_multinomial_element = std::find_if(par_sample_analysis.multinomial_data.v_p_peptide_analysis.begin(), par_sample_analysis.multinomial_data.v_p_peptide_analysis.end(),
				[itr_v_homology_data](peptide_analysis*& par_peptide_analysis) {
				return par_peptide_analysis->peptide_filtered == itr_v_homology_data.blastp_query;
			});
			const auto& find_multinomial_protein = std::find_if(par_sample_analysis.multinomial_data.v_p_protein_data.begin(), par_sample_analysis.multinomial_data.v_p_protein_data.end(),
				[itr_v_homology_data](protein_data*& par_protein_data) {
				return par_protein_data->protein_name == itr_v_homology_data.blastp_subject_accession;
			});
			if ((find_multinomial_element != par_sample_analysis.multinomial_data.v_p_peptide_analysis.end()) && (find_multinomial_protein) != par_sample_analysis.multinomial_data.v_p_protein_data.end()){
				const size_t i = (find_multinomial_element - par_sample_analysis.multinomial_data.v_p_peptide_analysis.begin());
				const size_t j = (find_multinomial_protein - par_sample_analysis.multinomial_data.v_p_protein_data.begin());
				par_sample_analysis.multinomial_data.v2_frequency[i][j] = itr_v_homology_data.blastp_homology_transformed;
				par_sample_analysis.multinomial_data.v_frequency_marginal_sum[i] += itr_v_homology_data.blastp_homology_transformed;
			}
			else {
				std::cout << "\n\n ~~~ query / peptide or subject / accession mismatch";
				std::cout << "\n\n ~~~ input any key to continue";
				string catch_error;
				std::cin >> catch_error;
			}
		}
	}

	inline bool predicate_multinomial_frequency(const multinomial_frequency_type& i, const multinomial_frequency_type& j) {
		return (i.multinomial_frequency > j.multinomial_frequency);
	}

	inline void sort_v_multinomial_frequency(vector<multinomial_frequency_type>& par_v_multinomial_frequency) {
		std::sort(par_v_multinomial_frequency.begin(), par_v_multinomial_frequency.end(), predicate_multinomial_frequency);
	}
}

#endif