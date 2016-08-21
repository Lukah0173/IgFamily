// * * fpf_multinomial.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_MULTINOMIAL
#define	FPF_MULTINOMIAL

#include <cstdlib>						// provides - size_t
#include <string>						// provides - string
#include <vector>						// provides - vector
#include <iomanip>						// provides - std::setprecision
#include <iostream>						// provides - std::get
#include <utility>						// provides - pair

#include "fpf_data.h"
#include "fpf_filesystem.h"


namespace fpf_multinomial {

	using std::string;
	using std::vector;

	typedef fpf_data::multinomial multinomial;
	typedef fpf_filesystem::filesystem filesystem;

	void create_filesystem_multinomial_data(filesystem& par_filesystem) {
		for (const auto& itr_v_category_analysis : par_filesystem.v_category_analysis_selected_by_polymorphism) {
			par_filesystem.multinomial_data.v_category_name.push_back(itr_v_category_analysis.category_name);
			par_filesystem.multinomial_data.v_category_class.push_back(itr_v_category_analysis.category_class);
		}
		for (const auto& itr_v_peptide_data : par_filesystem.v_peptide_data) {
			par_filesystem.multinomial_data.v_element_name.push_back(itr_v_peptide_data.peptide_filtered);
		}
		par_filesystem.multinomial_data.v2_frequency = vector<vector<double>>(par_filesystem.multinomial_data.v_element_name.size(), vector<double>(par_filesystem.multinomial_data.v_category_name.size(), 0));
		par_filesystem.multinomial_data.v2_density = vector<vector<double>>(par_filesystem.multinomial_data.v_element_name.size(), vector<double>(par_filesystem.multinomial_data.v_category_name.size(), 0));
		par_filesystem.multinomial_data.v_frequency_marginal_sum = vector<double>(par_filesystem.multinomial_data.v_element_name.size(), 0);

		for (const auto itr_v_blastp_data : par_filesystem.v_blastp_data) {
			const auto find_multinomial_element = std::find(par_filesystem.multinomial_data.v_element_name.begin(), par_filesystem.multinomial_data.v_element_name.end(), itr_v_blastp_data.blastp_query);		
			const auto find_multinomial_category = std::find(par_filesystem.multinomial_data.v_category_name.begin(), par_filesystem.multinomial_data.v_category_name.end(), itr_v_blastp_data.blastp_subject_accession);			
			
			if ((find_multinomial_element != par_filesystem.multinomial_data.v_element_name.end()) && (find_multinomial_category) != par_filesystem.multinomial_data.v_category_name.end()){
				const int i = (find_multinomial_element - par_filesystem.multinomial_data.v_element_name.begin());
				const int j = (find_multinomial_category - par_filesystem.multinomial_data.v_category_name.begin());
				par_filesystem.multinomial_data.v2_frequency[i][j] = itr_v_blastp_data.blastp_evalue_transformed;
				par_filesystem.multinomial_data.v_frequency_marginal_sum[i] += itr_v_blastp_data.blastp_evalue_transformed;
			}
			else {
				std::cout << "\n\n ~~~ query / peptide or subject / accession mismatch";
				std::cout << "\n\n ~~~ input any key to continue";
				string catch_error;
				std::cin >> catch_error;
			}
		}
	}

	void fout_multinomial(filesystem& par_filesystem) {
		string output_multinomial = par_filesystem.directory + "multinomial.csv";
		std::ofstream fout_multinomial;
		fout_multinomial.open(output_multinomial);
		std::cout << "\n\n ...outputting multinomial data frame for " << par_filesystem.filename;
		fout_multinomial << ",";
		fout_multinomial << "TOTAL,";
		for (const auto& itr_multinomial_category : par_filesystem.multinomial_data.v_category_name) {
			fout_multinomial << itr_multinomial_category << ",";
		}
		fout_multinomial << "\n";
		fout_multinomial << ",,";
		for (const auto& itr_multinomial_category_class : par_filesystem.multinomial_data.v_category_class) {
			fout_multinomial << itr_multinomial_category_class << ",";
		}
		fout_multinomial << "\n";
		for (auto i = 0; i < par_filesystem.multinomial_data.v_element_name.size(); ++i) {
			fout_multinomial << par_filesystem.multinomial_data.v_element_name[i] << ",";
			fout_multinomial << par_filesystem.multinomial_data.v_frequency_marginal_sum[i] << ",";
			for (auto j = 0; j < par_filesystem.multinomial_data.v_category_name.size(); ++j) {
				fout_multinomial << par_filesystem.multinomial_data.v2_frequency[i][j] << ",";
			}
			fout_multinomial << "\n";
		}
	}

	struct multinomial_frequency_type {
		string multinomial_category;
		double multinomial_frequency;
	};

	inline bool predicate_multinomial_frequency(const multinomial_frequency_type& i, const multinomial_frequency_type& j) {
		return (i.multinomial_frequency > j.multinomial_frequency);
	}

	inline void sort_v_multinomial_frequency(vector<multinomial_frequency_type>& par_v_multinomial_frequency) {
		std::sort(par_v_multinomial_frequency.begin(), par_v_multinomial_frequency.end(), predicate_multinomial_frequency);
	}

	void fout_multinomial_element(filesystem& par_filesystem) {
		string output_multinomial_element = par_filesystem.directory + "multinomial_peptide.txt";
		std::ofstream fout_multinomial_element;
		fout_multinomial_element.open(output_multinomial_element);
		std::cout << "\n\n ...outputting multinomial peptide list for " << par_filesystem.filename;
		fout_multinomial_element << "\n";
		fout_multinomial_element << par_filesystem.filename;
		fout_multinomial_element << "\n\n\n";
		int format_ws_length = int();
		for (auto i = 0; i < par_filesystem.multinomial_data.v_element_name.size(); ++i) {
			if (par_filesystem.multinomial_data.v_element_name[i].length() > format_ws_length) {
				format_ws_length = par_filesystem.multinomial_data.v_element_name[i].length();
			}
		}
		for (auto i = 0; i < par_filesystem.multinomial_data.v_element_name.size(); ++i) {
			fout_multinomial_element << par_filesystem.multinomial_data.v_element_name[i];
			for (auto j = 0; j < (format_ws_length - par_filesystem.multinomial_data.v_element_name[i].length() + 5); ++j) {
				fout_multinomial_element << " ";
			}
			vector<multinomial_frequency_type> temp_v_multinomial_frequency;
			for (auto j = 0; j < par_filesystem.multinomial_data.v_category_name.size(); ++j) {
				multinomial_frequency_type temp_multinomial_frequency;
				temp_multinomial_frequency.multinomial_category = par_filesystem.multinomial_data.v_category_name[j];
				temp_multinomial_frequency.multinomial_frequency = par_filesystem.multinomial_data.v2_frequency[i][j];
				temp_v_multinomial_frequency.push_back(temp_multinomial_frequency);
			}
			sort_v_multinomial_frequency(temp_v_multinomial_frequency);
			for (auto j = 0; j < par_filesystem.multinomial_data.v_category_name.size(); ++j) {
				if (temp_v_multinomial_frequency[j].multinomial_frequency > 0.2) {
					fout_multinomial_element << temp_v_multinomial_frequency[j].multinomial_category << " (";
					fout_multinomial_element << temp_v_multinomial_frequency[j].multinomial_frequency << "), ";
				}
			}
			fout_multinomial_element << "\n";
		}
	}

	void fout_multinomial_element_nomatch(filesystem& par_filesystem) {
		string output_multinomial_element_nomatch = par_filesystem.directory + "multinomial_peptide_filtered.txt";
		std::ofstream fout_multinomial_element_nomatch;
		fout_multinomial_element_nomatch.open(output_multinomial_element_nomatch);
		std::cout << "\n\n ...outputting filtered multinomial peptide list for " << par_filesystem.filename;
		fout_multinomial_element_nomatch << "\n";
		fout_multinomial_element_nomatch << par_filesystem.filename;
		fout_multinomial_element_nomatch << "\n\n\n";
		int format_ws_length = int();
		for (auto i = 0; i < par_filesystem.multinomial_data.v_element_name.size(); ++i) {
			if (par_filesystem.multinomial_data.v_element_name[i].length() > format_ws_length) {
				format_ws_length = par_filesystem.multinomial_data.v_element_name[i].length();
			}
		}
		for (auto i = 0; i < par_filesystem.multinomial_data.v_element_name.size(); ++i) {
			double format_frequency_threshold = double();
			for (auto j = 0; j < par_filesystem.multinomial_data.v_category_name.size(); ++j) {
				if (par_filesystem.multinomial_data.v2_frequency[i][j] > format_frequency_threshold) {
					format_frequency_threshold = par_filesystem.multinomial_data.v2_frequency[i][j];
				}
			}
			if (format_frequency_threshold < 0.1) {
				fout_multinomial_element_nomatch << par_filesystem.multinomial_data.v_element_name[i];
				for (auto j = 0; j < (format_ws_length - par_filesystem.multinomial_data.v_element_name[i].length() + 5); ++j) {
					fout_multinomial_element_nomatch << " ";
				}
				vector<multinomial_frequency_type> temp_v_multinomial_frequency;
				for (auto j = 0; j < par_filesystem.multinomial_data.v_category_name.size(); ++j) {
					multinomial_frequency_type temp_multinomial_frequency;
					temp_multinomial_frequency.multinomial_category = par_filesystem.multinomial_data.v_category_name[j];
					temp_multinomial_frequency.multinomial_frequency = par_filesystem.multinomial_data.v2_frequency[i][j];
					temp_v_multinomial_frequency.push_back(temp_multinomial_frequency);
				}
				sort_v_multinomial_frequency(temp_v_multinomial_frequency);
				for (auto j = 0; j < par_filesystem.multinomial_data.v_category_name.size(); ++j) {
					if (temp_v_multinomial_frequency[j].multinomial_frequency > 0.1) {
						fout_multinomial_element_nomatch << temp_v_multinomial_frequency[j].multinomial_category << " (";
						fout_multinomial_element_nomatch << temp_v_multinomial_frequency[j].multinomial_frequency << "), ";
					}
				}
				fout_multinomial_element_nomatch << "\n";
			}
		}
	}
}

#endif