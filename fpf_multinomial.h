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

	typedef fpf_data::multinomial_frequency_type multinomial_frequency_type;
	typedef fpf_filesystem::filesystem filesystem;
	typedef fpf_filesystem::sample_analysis sample_analysis;

	void create_filesystem_multinomial_data(sample_analysis& par_sample_analysis) {
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis_selected_by_polymorphism) {
			par_sample_analysis.multinomial_data.v_protein_name.push_back(itr_v_protein_analysis.p_protein_data->protein_name);
			par_sample_analysis.multinomial_data.v_protein_class.push_back(itr_v_protein_analysis.p_protein_data->protein_class);
		}
		for (const auto& itr_v_peptide_analysis : par_sample_analysis.v_peptide_analysis) {
			par_sample_analysis.multinomial_data.v_element_name.push_back(itr_v_peptide_analysis.peptide_filtered);
		}
		par_sample_analysis.multinomial_data.v2_frequency = vector<vector<double>>(par_sample_analysis.multinomial_data.v_element_name.size(), vector<double>(par_sample_analysis.multinomial_data.v_protein_name.size(), 0));
		par_sample_analysis.multinomial_data.v2_density = vector<vector<double>>(par_sample_analysis.multinomial_data.v_element_name.size(), vector<double>(par_sample_analysis.multinomial_data.v_protein_name.size(), 0));
		par_sample_analysis.multinomial_data.v_frequency_marginal_sum = vector<double>(par_sample_analysis.multinomial_data.v_element_name.size(), 0);

		for (const auto itr_v_homology_data : par_sample_analysis.v_homology_data) {
			const auto find_multinomial_element = std::find(par_sample_analysis.multinomial_data.v_element_name.begin(), par_sample_analysis.multinomial_data.v_element_name.end(), itr_v_homology_data.blastp_query);		
			const auto find_multinomial_protein = std::find(par_sample_analysis.multinomial_data.v_protein_name.begin(), par_sample_analysis.multinomial_data.v_protein_name.end(), itr_v_homology_data.blastp_subject_accession);			
			
			if ((find_multinomial_element != par_sample_analysis.multinomial_data.v_element_name.end()) && (find_multinomial_protein) != par_sample_analysis.multinomial_data.v_protein_name.end()){
				const size_t i = (find_multinomial_element - par_sample_analysis.multinomial_data.v_element_name.begin());
				const size_t j = (find_multinomial_protein - par_sample_analysis.multinomial_data.v_protein_name.begin());
				par_sample_analysis.multinomial_data.v2_frequency[i][j] = itr_v_homology_data.blastp_evalue_transformed_conjugated;
				par_sample_analysis.multinomial_data.v_frequency_marginal_sum[i] += itr_v_homology_data.blastp_evalue_transformed_conjugated;
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