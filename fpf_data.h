// * * fpf_data.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_DATA
#define	FPF_DATA
#include <cstdlib> // provides - size_t
#include <vector> // provides - std::vector
#include <sstream> // provides - std::istringstream
#include <algorithm> // provides - std::find
#include <string> // provides - std::string, std::string::pop_back
#include <complex> // provides - std::log
#include <math.h> // provides - std::floor
#include <assert.h> 
#include "IgFamily.h"
#include "fpf_parse.h"



namespace fpf_data {

	struct multinomial_category_data_type;
	struct peptide_data_type;
	struct denovo_aminoacid_type;
	struct denovo_peptide_type;
	struct blastp_type;
	struct proteinconstruct_from_denovo_type;
	struct report_type;
	struct multinomial_type;

	typedef std::string string_type;
	typedef size_t size_type;
	typedef bool bool_type;

	struct multinomial_category_data_type {
	public:
		multinomial_category_data_type() {
			d_score = { 1 };
		};

		string_type str_multinomial_category_class;
		string_type str_multinomial_category_name;
		string_type str_protein;
		string_type str_species;
		std::vector<peptide_data_type> v_s_peptide_data;
		std::vector<peptide_data_type> v_s_peptide_data_distinct_filtered;
		double d_coverage;
		double d_score;
		std::vector<multinomial_category_data_type*> v_s_multinomial_category_polyassociation;
	};

	struct denovo_aminoacid_type {
		char ch_aminoacid;
		double d_denovo_localconfidence;
	};

	struct denovo_peptide_type {
		std::vector<denovo_aminoacid_type> v_s_denovo_aminoacid;
		double d_denovo_peptide_localconfidence_average;
	};

	struct peptide_data_type {
	public:
		string_type str_peptide;
		string_type str_peptide_filtered;
		size_type st_spectralcount;
		size_type st_IgP;
		size_type st_count_denovo_replicate;
		std::vector<denovo_peptide_type> v_s_denovo_peptide;
		bool b_replicate_merged = bool();
		size_type st_filesystem_replicate = size_type{ 1 };
		std::vector<std::tuple<string_type, size_type, size_type>> v_p_replicate_data;
		std::vector<string_type> v_str_peptideassociation;
		std::vector<string_type> v_str_peptideassociation_distinct;
		std::vector<std::pair<multinomial_category_data_type*, double>> v_p_peptideassociation;
		std::vector<std::pair<multinomial_category_data_type*, double>> v_p_peptideassociation_distinct;
	};

	struct blastp_type {
	public:
		string_type str_blastp_query;
		string_type str_blastp_query_aligned;
		string_type str_blastp_subject;
		string_type str_blastp_subject_accession;
		string_type str_blastp_subject_accession_class;
		size_type st_blastp_query_alignment_index;
		size_type st_blastp_subject_alignment_index;
		string_type str_protein;
		string_type str_blastp_query_alignment;
		double d_blastp_evalue;
		double d_blastp_par_prop;
		denovo_peptide_type s_denovo_peptide_best;
		size_type st_count_denovo_replicates;
	};

	struct proteinconstruct_from_denovo_type {
		char ch_aminoacid;
		double d_score;
	};

	struct report_type {
	public:
		std::vector<blastp_type> v_s_blastp_genefamily_combined;
		string_type str_report_multinomial_category_class;
		string_type str_protein_accession;
		string_type str_protein;
		std::vector<proteinconstruct_from_denovo_type> proteinconstruct_from_denovo_type;
		double d_score;
	};

	struct multinomial_type {
	public:
		std::vector<string_type> v_str_multinomial_element;
		std::vector<string_type> v_str_multinomial_category;
		std::vector<string_type> v_str_multinomial_category_class;
		std::vector<std::vector<double>> v2_d_multinomial_frequency;
		std::vector<double> v_d_multinomial_frequency_sum;
		std::vector<std::vector<double>> v2_d_multinomial_density;
	};

	std::vector<multinomial_category_data_type> create_v_s_multinomial_element_data(std::vector<fpf_parse::parse_FASTA_type> par_v_c_parse_FASTA) {

		// Vector reallocation will invalidate all pointers to multinomial_category_data_type references.
		// Do not reallocate after return!

		std::vector<multinomial_category_data_type> con_v_s_multinomial_element_data;
		multinomial_category_data_type con_c_analysis;
		size_type con_st_ID = size_type();
		for (auto itr_par_v_c_parse_FASTA : par_v_c_parse_FASTA) {
			++con_st_ID;
			auto find_v_s_multinomial_element_data = std::find_if(con_v_s_multinomial_element_data.begin(), con_v_s_multinomial_element_data.end(),
				[itr_par_v_c_parse_FASTA](multinomial_category_data_type par_s_multinomial_element_data) {
				return par_s_multinomial_element_data.str_multinomial_category_name == itr_par_v_c_parse_FASTA.return_str_parse_FASTA_genefamily(); });
			if (find_v_s_multinomial_element_data == con_v_s_multinomial_element_data.end()) {
				con_c_analysis.str_multinomial_category_name = itr_par_v_c_parse_FASTA.return_str_parse_FASTA_genefamily();
				con_c_analysis.str_multinomial_category_class = itr_par_v_c_parse_FASTA.return_str_parse_FASTA_genefamily_class();
				con_c_analysis.str_species = itr_par_v_c_parse_FASTA.return_str_parse_FASTA_species();
				con_c_analysis.str_protein = itr_par_v_c_parse_FASTA.return_str_parse_FASTA_protein();
				con_v_s_multinomial_element_data.push_back(con_c_analysis);
			}
			else {
				find_v_s_multinomial_element_data->str_protein = find_v_s_multinomial_element_data->str_multinomial_category_name + itr_par_v_c_parse_FASTA.return_str_parse_FASTA_protein();
			}
		}
		return con_v_s_multinomial_element_data;
	}

	std::vector<multinomial_category_data_type> create_v_s_multinomial_element_data_distinct(std::vector<multinomial_category_data_type>& par_v_s_multinomial_element_data) {
		std::vector<multinomial_category_data_type> con_v_s_multinomial_element_data_distinct;
		std::vector<string_type> con_v_str_multinomial_category_name_distinct;
		multinomial_category_data_type con_c_analysis_distinct;
		for (auto itr_v_s_multinomial_element_data_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data_data) {
			string_type con_str_multinomial_category_name = itr_v_s_multinomial_element_data_data->str_multinomial_category_name;
			size_type sw_peptideassociation_distinct = size_type();
			string_type con_str_multinomial_category_name_distinct = string_type();
			string_type con_str_multinomial_category_name_distinct_class = string_type();
			for (auto itr_str_multinomial_category_name = size_type(); itr_str_multinomial_category_name < con_str_multinomial_category_name.length(); ++itr_str_multinomial_category_name) {
				if (con_str_multinomial_category_name.at(itr_str_multinomial_category_name) == '*') {
					sw_peptideassociation_distinct = 1;
				}
				if (sw_peptideassociation_distinct == 0) {
					con_str_multinomial_category_name_distinct += con_str_multinomial_category_name.at(itr_str_multinomial_category_name);
				}
			}
			con_str_multinomial_category_name_distinct_class = itr_v_s_multinomial_element_data_data->str_multinomial_category_class;
			if (std::find(con_v_str_multinomial_category_name_distinct.begin(), con_v_str_multinomial_category_name_distinct.end(), con_str_multinomial_category_name_distinct) == con_v_str_multinomial_category_name_distinct.end()) {
				con_v_str_multinomial_category_name_distinct.push_back(con_str_multinomial_category_name_distinct);
				con_c_analysis_distinct.str_multinomial_category_name = con_str_multinomial_category_name_distinct;
				con_c_analysis_distinct.str_multinomial_category_class = con_str_multinomial_category_name_distinct_class;
				con_c_analysis_distinct.str_species = itr_v_s_multinomial_element_data_data->str_species;
				con_v_s_multinomial_element_data_distinct.push_back(con_c_analysis_distinct);
			}
			for (auto itr_v_s_multinomial_element_data_distinct = con_v_s_multinomial_element_data_distinct.begin(); itr_v_s_multinomial_element_data_distinct != con_v_s_multinomial_element_data_distinct.end(); ++itr_v_s_multinomial_element_data_distinct) {
				if (itr_v_s_multinomial_element_data_distinct->str_multinomial_category_name == con_str_multinomial_category_name_distinct) {
					itr_v_s_multinomial_element_data_distinct->v_s_multinomial_category_polyassociation.push_back(&(*itr_v_s_multinomial_element_data_data));
				}
			}
			sw_peptideassociation_distinct = 0;
			con_str_multinomial_category_name_distinct.clear();
		}
		return con_v_s_multinomial_element_data_distinct;
	}

	void update_v_s_multinomial_element_distinctpolymorphism_data(std::vector<multinomial_category_data_type>& par_v_s_multinomial_element_data_distinct) {
		multinomial_category_data_type con_c_analysis_distinct_update;
		for (auto itr_v_s_multinomial_element_data_distinct = par_v_s_multinomial_element_data_distinct.begin(); itr_v_s_multinomial_element_data_distinct != par_v_s_multinomial_element_data_distinct.end(); ++itr_v_s_multinomial_element_data_distinct) {
			std::vector<multinomial_category_data_type*> con_v_s_multinomial_element_data_distinct_update = itr_v_s_multinomial_element_data_distinct->v_s_multinomial_category_polyassociation;
			for (auto itr_v_s_multinomial_element_data_distinct_polyassociation = con_v_s_multinomial_element_data_distinct_update.begin(); itr_v_s_multinomial_element_data_distinct_polyassociation != con_v_s_multinomial_element_data_distinct_update.end(); ++itr_v_s_multinomial_element_data_distinct_polyassociation) {
				if (itr_v_s_multinomial_element_data_distinct_polyassociation == con_v_s_multinomial_element_data_distinct_update.begin()) {
					con_c_analysis_distinct_update = **itr_v_s_multinomial_element_data_distinct_polyassociation;
				}
				else {
					if (con_c_analysis_distinct_update.d_score < (*itr_v_s_multinomial_element_data_distinct_polyassociation)->d_score) {
						con_c_analysis_distinct_update = **itr_v_s_multinomial_element_data_distinct_polyassociation;
					}
				}
			}
			con_c_analysis_distinct_update.v_s_multinomial_category_polyassociation = itr_v_s_multinomial_element_data_distinct->v_s_multinomial_category_polyassociation;
			con_c_analysis_distinct_update.str_multinomial_category_name = itr_v_s_multinomial_element_data_distinct->str_multinomial_category_name;
			con_c_analysis_distinct_update.str_multinomial_category_class = itr_v_s_multinomial_element_data_distinct->str_multinomial_category_class;
			con_c_analysis_distinct_update.str_species = itr_v_s_multinomial_element_data_distinct->str_species;
			*itr_v_s_multinomial_element_data_distinct = con_c_analysis_distinct_update;
		}
	}

	std::vector<peptide_data_type> create_v_s_peptide_data(std::vector<fpf_parse::parse_peptides_csv_type> par_c_parse_csv_peptide_data) {
		std::vector<peptide_data_type> con_v_s_peptide_data;
		for (auto itr_c_parse_csv_peptide_data : par_c_parse_csv_peptide_data) {
			peptide_data_type con_s_peptide_data = peptide_data_type();
			denovo_peptide_type con_s_denovo_peptide = denovo_peptide_type();
			denovo_aminoacid_type con_s_denovo_aminoacid = denovo_aminoacid_type();
			string_type con_str_peptide_filtered = string_type();
			size_type sw_peptide_filtered = size_type();
			for (auto j = size_type(); j < itr_c_parse_csv_peptide_data.str_parse_peptides_csv_peptide.length(); ++j) {
				if (itr_c_parse_csv_peptide_data.str_parse_peptides_csv_peptide.at(j) == '(' || ((itr_c_parse_csv_peptide_data.str_parse_peptides_csv_peptide.at(j) == '.') && (itr_c_parse_csv_peptide_data.str_parse_peptides_csv_peptide.length() > 2))) {
					sw_peptide_filtered = 1;
					if (itr_c_parse_csv_peptide_data.str_parse_peptides_csv_peptide.at(j + 1) == 's') {
						con_str_peptide_filtered.pop_back();
						con_str_peptide_filtered += itr_c_parse_csv_peptide_data.str_parse_peptides_csv_peptide.at(j + 5);
					}
				}
				if (sw_peptide_filtered == 0) {
					con_str_peptide_filtered += itr_c_parse_csv_peptide_data.str_parse_peptides_csv_peptide.at(j);
				}
				if (itr_c_parse_csv_peptide_data.str_parse_peptides_csv_peptide.at(j) == ')') {
					sw_peptide_filtered = 0;
				}
				if ((itr_c_parse_csv_peptide_data.str_parse_peptides_csv_peptide.at(j) == '.') && (itr_c_parse_csv_peptide_data.str_parse_peptides_csv_peptide.length() <= 2) && (sw_peptide_filtered == 0)) {
					con_str_peptide_filtered.clear();
				}
			}
			con_s_peptide_data.str_peptide_filtered = con_str_peptide_filtered;
			auto& find_s_peptide_data = std::find_if(con_v_s_peptide_data.begin(), con_v_s_peptide_data.end(),
				[con_str_peptide_filtered](peptide_data_type par_s_peptide_data) {
				return par_s_peptide_data.str_peptide_filtered == con_str_peptide_filtered;
			});
			if (find_s_peptide_data == con_v_s_peptide_data.end()) {
				size_type ss_st_csv_IgP = size_type();
				std::istringstream(itr_c_parse_csv_peptide_data.str_parse_csv_IgP) >> ss_st_csv_IgP;
				con_s_peptide_data.str_peptide = itr_c_parse_csv_peptide_data.str_parse_peptides_csv_peptide;
				std::istringstream ss_spectralcount(itr_c_parse_csv_peptide_data.str_parse_peptides_csv_spectralcount);
				size_type ss_st_spectralcount;
				ss_spectralcount >> ss_st_spectralcount;
				con_s_peptide_data.st_spectralcount = ss_st_spectralcount;
				for (size_type i = 0; i < con_str_peptide_filtered.size(); ++i) {
					con_s_denovo_aminoacid.ch_aminoacid = con_str_peptide_filtered[i];
					con_s_denovo_aminoacid.d_denovo_localconfidence = itr_c_parse_csv_peptide_data.v_d_denovo_localconfidence[i];
					con_s_denovo_peptide.v_s_denovo_aminoacid.push_back(con_s_denovo_aminoacid);
				}
				con_s_denovo_peptide.d_denovo_peptide_localconfidence_average = double();
				for (auto itr_s_denovo_aminoacid : con_s_denovo_peptide.v_s_denovo_aminoacid) {
					con_s_denovo_peptide.d_denovo_peptide_localconfidence_average += itr_s_denovo_aminoacid.d_denovo_localconfidence;
				}
				con_s_denovo_peptide.d_denovo_peptide_localconfidence_average /= con_s_denovo_peptide.v_s_denovo_aminoacid.size();
				con_s_peptide_data.v_s_denovo_peptide.push_back(con_s_denovo_peptide);
				++con_s_peptide_data.st_count_denovo_replicate;
				con_v_s_peptide_data.push_back(con_s_peptide_data);
			}
			else {
				for (size_type i = 0; i < con_str_peptide_filtered.size(); ++i) {
					con_s_denovo_aminoacid.ch_aminoacid = con_str_peptide_filtered[i];
					con_s_denovo_aminoacid.d_denovo_localconfidence = itr_c_parse_csv_peptide_data.v_d_denovo_localconfidence[i];
					con_s_denovo_peptide.v_s_denovo_aminoacid.push_back(con_s_denovo_aminoacid);
				}
				con_s_denovo_peptide.d_denovo_peptide_localconfidence_average = double();
				for (auto itr_s_denovo_aminoacid : con_s_denovo_peptide.v_s_denovo_aminoacid) {
					con_s_denovo_peptide.d_denovo_peptide_localconfidence_average += itr_s_denovo_aminoacid.d_denovo_localconfidence;
				}
				con_s_denovo_peptide.d_denovo_peptide_localconfidence_average /= con_s_denovo_peptide.v_s_denovo_aminoacid.size();
				find_s_peptide_data->v_s_denovo_peptide.push_back(con_s_denovo_peptide);
				++find_s_peptide_data->st_count_denovo_replicate;
			}
		}
		return con_v_s_peptide_data;
	}

	void create_v_s_multinomial_element_data_v_s_peptide_data(std::vector<multinomial_category_data_type>& par_v_s_multinomial_element_data, std::vector<peptide_data_type> par_v_s_peptide_data) {
		string_type con_str_peptide = string_type();
		std::vector<peptide_data_type> con_v_s_peptide_data;
		size_type sw_peptide_filtered = size_type();
		for (auto& itr_v_s_multinomial_element_data : par_v_s_multinomial_element_data) {
			for (auto itr_par_v_data : par_v_s_peptide_data) {
				for (auto j = size_type(); j < itr_par_v_data.str_peptide.length(); ++j) {
					if (itr_par_v_data.str_peptide.at(j) == '(' || ((itr_par_v_data.str_peptide.at(j) == '.') && (itr_par_v_data.str_peptide.length() > 2))) {
						sw_peptide_filtered = 1;
						if (itr_par_v_data.str_peptide.at(j + 1) == 's') {
							con_str_peptide.pop_back();
							con_str_peptide += itr_par_v_data.str_peptide.at(j + 5);
						}
					}
					if (sw_peptide_filtered == 0) {
						con_str_peptide += itr_par_v_data.str_peptide.at(j);
					}
					if (itr_par_v_data.str_peptide.at(j) == ')') {
						sw_peptide_filtered = 0;
					}
					if ((itr_par_v_data.str_peptide.at(j) == '.') && (itr_par_v_data.str_peptide.length() <= 2) && (sw_peptide_filtered == 0)) {
						con_str_peptide.clear();
					}
				}
				if (itr_v_s_multinomial_element_data.str_protein.find(con_str_peptide) != std::string::npos) {
					con_v_s_peptide_data.push_back(itr_par_v_data);
				}
				con_str_peptide.clear();
				sw_peptide_filtered = 0;
			}
			itr_v_s_multinomial_element_data.v_s_peptide_data = con_v_s_peptide_data;
			con_v_s_peptide_data.clear();
			if ((IgFamily::DEBUG_MODE == 2) && (itr_v_s_multinomial_element_data.v_s_peptide_data.size() != 0)) {
				std::cout << "\n\n" << itr_v_s_multinomial_element_data.str_multinomial_category_name;
				std::cout << "   " << itr_v_s_multinomial_element_data.str_protein;
				std::vector<peptide_data_type> con_itr_v_s_peptide_data = itr_v_s_multinomial_element_data.v_s_peptide_data;
				for (auto itr_v_s_peptide_data : con_itr_v_s_peptide_data) {
					std::cout << "\n * " << itr_v_s_peptide_data.str_peptide;
					std::cout << "   " << itr_v_s_peptide_data.st_spectralcount;
				}
			}
		}
	}

	void create_v_s_multinomial_element_data_v_s_peptide_data_distinct_filtered(std::vector<multinomial_category_data_type>& par_v_s_multinomial_element_data, std::vector<peptide_data_type> par_v_data_filtered_distinct) {
		string_type con_str_peptide_distinct_filtered = string_type();
		std::vector<peptide_data_type> con_v_s_peptide_data_distinct_filtered;
		for (auto& itr_v_s_multinomial_element_data : par_v_s_multinomial_element_data) {
			for (auto itr_v_data_filtered_distinct : par_v_data_filtered_distinct) {
				if (itr_v_s_multinomial_element_data.str_protein.find(itr_v_data_filtered_distinct.str_peptide) != std::string::npos) {
					peptide_data_type con_s_peptide_data_distinct_filtered;
					con_s_peptide_data_distinct_filtered = itr_v_data_filtered_distinct;
					con_v_s_peptide_data_distinct_filtered.push_back(con_s_peptide_data_distinct_filtered);
				}
			}
			itr_v_s_multinomial_element_data.v_s_peptide_data_distinct_filtered = con_v_s_peptide_data_distinct_filtered;
			con_v_s_peptide_data_distinct_filtered.clear();
		}
	}
}

#endif