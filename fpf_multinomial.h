// * * fpf_multinomial.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_MULTINOMIAL
#define	FPF_MULTINOMIAL
#include <cstdlib>						// provides - size_t
#include <string>						// provides - std::string
#include <vector>						// provides - std::vector
#include <iomanip>						// provides - std::setprecision
#include <iostream>						// provides - std::get
#include <utility>						// provides - std::pair
#include "fpf_data.h"
#include "fpf_filesystem.h"



namespace fpf_multinomial {

	typedef size_t size_type;
	typedef std::string string_type;
	typedef fpf_data::multinomial_type multinomial_type;
	typedef fpf_filesystem::filesystem_type filesystem_type;

	void create_s_filesystem_mnom(filesystem_type& par_s_filesystem) {
		for (auto itr_v_c_multinomial_catagory : par_s_filesystem.v_c_multinomial_catagory) {
			par_s_filesystem.s_multinomial.v_str_multinomial_category.push_back(itr_v_c_multinomial_catagory.str_multinomial_catagory_name);
			par_s_filesystem.s_multinomial.v_str_multinomial_category_class.push_back(itr_v_c_multinomial_catagory.str_multinomial_catagory_class);
		}
		for (auto itr_v_s_peptide_data : par_s_filesystem.v_s_peptide_data) {
			par_s_filesystem.s_multinomial.v_str_multinomial_element.push_back(itr_v_s_peptide_data.str_peptide_filtered);
		}
		par_s_filesystem.s_multinomial.v2_d_multinomial_frequency = std::vector<std::vector<double>>(par_s_filesystem.s_multinomial.v_str_multinomial_element.size(), std::vector<double>(par_s_filesystem.s_multinomial.v_str_multinomial_category.size(), 0));
		par_s_filesystem.s_multinomial.v2_d_multinomial_density = std::vector<std::vector<double>>(par_s_filesystem.s_multinomial.v_str_multinomial_element.size(), std::vector<double>(par_s_filesystem.s_multinomial.v_str_multinomial_category.size(), 0));
		par_s_filesystem.s_multinomial.v_d_multinomial_frequency_sum = std::vector<double>(par_s_filesystem.s_multinomial.v_str_multinomial_element.size(), 0);

		for (auto itr_v_s_blastp : par_s_filesystem.v_s_blastp) {
			auto find_v_s_multinomial_element = std::find(par_s_filesystem.s_multinomial.v_str_multinomial_element.begin(), par_s_filesystem.s_multinomial.v_str_multinomial_element.end(), itr_v_s_blastp.str_blastp_query);		
			auto find_v_s_multinomial_category = std::find(par_s_filesystem.s_multinomial.v_str_multinomial_category.begin(), par_s_filesystem.s_multinomial.v_str_multinomial_category.end(), itr_v_s_blastp.str_blastp_subject_accession);			
			
			if ((find_v_s_multinomial_element != par_s_filesystem.s_multinomial.v_str_multinomial_element.end()) && (find_v_s_multinomial_category) != par_s_filesystem.s_multinomial.v_str_multinomial_category.end()){
				int i = (find_v_s_multinomial_element - par_s_filesystem.s_multinomial.v_str_multinomial_element.begin());
				int j = (find_v_s_multinomial_category - par_s_filesystem.s_multinomial.v_str_multinomial_category.begin());
				par_s_filesystem.s_multinomial.v2_d_multinomial_frequency[i][j] = itr_v_s_blastp.d_blastp_par_prop;
				par_s_filesystem.s_multinomial.v_d_multinomial_frequency_sum[i] += itr_v_s_blastp.d_blastp_par_prop;
			}
			else {
				std::cout << "\n\n ~~~ query / peptide or subject / accession mismatch";
				std::cout << "\n\n ~~~ input any key to continue";
				string_type catch_error;
				//std::cin >> catch_error;
			}
		}
	}

	void fout_s_multinomial(filesystem_type& par_s_filesystem) {
		std::string output_s_multinomial = par_s_filesystem.str_directory + "multinomial.csv";
		std::ofstream fout_s_multinomial;
		fout_s_multinomial.open(output_s_multinomial);
		std::cout << "\n\n ...outputting multinomial data frame for " << par_s_filesystem.str_filename;
		fout_s_multinomial << ",";
		fout_s_multinomial << "TOTAL,";
		for (auto itr_v_str_multinomial_category : par_s_filesystem.s_multinomial.v_str_multinomial_category) {
			fout_s_multinomial << itr_v_str_multinomial_category << ",";
		}
		fout_s_multinomial << "\n";
		fout_s_multinomial << ",,";
		for (auto itr_v_str_multinomial_category_class : par_s_filesystem.s_multinomial.v_str_multinomial_category_class) {
			fout_s_multinomial << itr_v_str_multinomial_category_class << ",";
		}
		fout_s_multinomial << "\n";
		for (auto i = 0; i < par_s_filesystem.s_multinomial.v_str_multinomial_element.size(); ++i) {
			fout_s_multinomial << par_s_filesystem.s_multinomial.v_str_multinomial_element[i] << ",";
			fout_s_multinomial << par_s_filesystem.s_multinomial.v_d_multinomial_frequency_sum[i] << ",";
			for (auto j = 0; j < par_s_filesystem.s_multinomial.v_str_multinomial_category.size(); ++j) {
				fout_s_multinomial << par_s_filesystem.s_multinomial.v2_d_multinomial_frequency[i][j] << ",";
			}
			fout_s_multinomial << "\n";
		}
	}

	struct multinomial_frequency_type {
		string_type str_multinomial_category;
		double d_multinomial_frequency;
	};

	inline bool predicate_v_s_multinomial_frequency(const multinomial_frequency_type& i, const multinomial_frequency_type& j) {
		return (i.d_multinomial_frequency > j.d_multinomial_frequency);
	}

	inline void sort_v_s_multinomial_frequency(std::vector<multinomial_frequency_type>& par_v_s_multinomial_frequency) {
		std::sort(par_v_s_multinomial_frequency.begin(), par_v_s_multinomial_frequency.end(), predicate_v_s_multinomial_frequency);
	}

	void fout_s_multinomial_peptide_list(filesystem_type& par_s_filesystem) {
		std::string output_s_multinomial = par_s_filesystem.str_directory + "multinomial_peptide.txt";
		std::ofstream fout_s_multinomial;
		fout_s_multinomial.open(output_s_multinomial);
		std::cout << "\n\n ...outputting multinomial peptide list for " << par_s_filesystem.str_filename;
		fout_s_multinomial << "\n";
		fout_s_multinomial << par_s_filesystem.str_filename;
		fout_s_multinomial << "\n\n\n";
		int i_format_ws_length = int();
		for (auto i = 0; i < par_s_filesystem.s_multinomial.v_str_multinomial_element.size(); ++i) {
			if (par_s_filesystem.s_multinomial.v_str_multinomial_element[i].length() > i_format_ws_length) {
				i_format_ws_length = par_s_filesystem.s_multinomial.v_str_multinomial_element[i].length();
			}
		}
		for (auto i = 0; i < par_s_filesystem.s_multinomial.v_str_multinomial_element.size(); ++i) {
			fout_s_multinomial << par_s_filesystem.s_multinomial.v_str_multinomial_element[i];
			for (auto j = 0; j < (i_format_ws_length - par_s_filesystem.s_multinomial.v_str_multinomial_element[i].length() + 5); ++j) {
				fout_s_multinomial << " ";
			}
			std::vector<multinomial_frequency_type> con_v_s_multinomial_frequency;
			for (auto j = 0; j < par_s_filesystem.s_multinomial.v_str_multinomial_category.size(); ++j) {
				multinomial_frequency_type con_s_multinomial_frequency;
				con_s_multinomial_frequency.str_multinomial_category = par_s_filesystem.s_multinomial.v_str_multinomial_category[j];
				con_s_multinomial_frequency.d_multinomial_frequency = par_s_filesystem.s_multinomial.v2_d_multinomial_frequency[i][j];
				con_v_s_multinomial_frequency.push_back(con_s_multinomial_frequency);
			}
			sort_v_s_multinomial_frequency(con_v_s_multinomial_frequency);
			for (auto j = 0; j < par_s_filesystem.s_multinomial.v_str_multinomial_category.size(); ++j) {
				if (con_v_s_multinomial_frequency[j].d_multinomial_frequency > 0.2) {
					//fout_s_multinomial << par_s_filesystem.s_multinomial.v_str_multinomial_category[j] << " (";
					//fout_s_multinomial << par_s_filesystem.s_multinomial.v2_d_multinomial_frequency[i][j] << "), ";
					fout_s_multinomial << con_v_s_multinomial_frequency[j].str_multinomial_category << " (";
					fout_s_multinomial << con_v_s_multinomial_frequency[j].d_multinomial_frequency << "), ";
				}
			}
			fout_s_multinomial << "\n";
		}
	}

	void fout_s_multinomial_peptide_list_filtered(filesystem_type& par_s_filesystem) {
		std::string output_s_multinomial = par_s_filesystem.str_directory + "multinomial_peptide_filtered.txt";
		std::ofstream fout_s_multinomial;
		fout_s_multinomial.open(output_s_multinomial);
		std::cout << "\n\n ...outputting filtered multinomial peptide list for " << par_s_filesystem.str_filename;
		fout_s_multinomial << "\n";
		fout_s_multinomial << par_s_filesystem.str_filename;
		fout_s_multinomial << "\n\n\n";
		int i_format_ws_length = int();
		for (auto i = 0; i < par_s_filesystem.s_multinomial.v_str_multinomial_element.size(); ++i) {
			if (par_s_filesystem.s_multinomial.v_str_multinomial_element[i].length() > i_format_ws_length) {
				i_format_ws_length = par_s_filesystem.s_multinomial.v_str_multinomial_element[i].length();
			}
		}
		for (auto i = 0; i < par_s_filesystem.s_multinomial.v_str_multinomial_element.size(); ++i) {
			double d_format_frequency_threshold = double();
			for (auto j = 0; j < par_s_filesystem.s_multinomial.v_str_multinomial_category.size(); ++j) {
				if (par_s_filesystem.s_multinomial.v2_d_multinomial_frequency[i][j] > d_format_frequency_threshold) {
					d_format_frequency_threshold = par_s_filesystem.s_multinomial.v2_d_multinomial_frequency[i][j];
				}
			}
			if (d_format_frequency_threshold < 0.1) {
				fout_s_multinomial << par_s_filesystem.s_multinomial.v_str_multinomial_element[i];
				for (auto j = 0; j < (i_format_ws_length - par_s_filesystem.s_multinomial.v_str_multinomial_element[i].length() + 5); ++j) {
					fout_s_multinomial << " ";
				}
				std::vector<multinomial_frequency_type> con_v_s_multinomial_frequency;
				for (auto j = 0; j < par_s_filesystem.s_multinomial.v_str_multinomial_category.size(); ++j) {
					multinomial_frequency_type con_s_multinomial_frequency;
					con_s_multinomial_frequency.str_multinomial_category = par_s_filesystem.s_multinomial.v_str_multinomial_category[j];
					con_s_multinomial_frequency.d_multinomial_frequency = par_s_filesystem.s_multinomial.v2_d_multinomial_frequency[i][j];
					con_v_s_multinomial_frequency.push_back(con_s_multinomial_frequency);
				}
				sort_v_s_multinomial_frequency(con_v_s_multinomial_frequency);
				for (auto j = 0; j < par_s_filesystem.s_multinomial.v_str_multinomial_category.size(); ++j) {
					if (con_v_s_multinomial_frequency[j].d_multinomial_frequency > 0.1) {
						//fout_s_multinomial << par_s_filesystem.s_multinomial.v_str_multinomial_category[j] << " (";
						//fout_s_multinomial << par_s_filesystem.s_multinomial.v2_d_multinomial_frequency[i][j] << "), ";
						fout_s_multinomial << con_v_s_multinomial_frequency[j].str_multinomial_category << " (";
						fout_s_multinomial << con_v_s_multinomial_frequency[j].d_multinomial_frequency << "), ";
					}
				}
				fout_s_multinomial << "\n";
			}
		}
	}
}

#endif