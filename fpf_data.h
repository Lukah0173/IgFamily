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
#include <string> // provides - std::string, string::pop_back
#include <math.h> // provides - std::floor

#include "IgFamily.h"
#include "fpf_parse.h"


namespace fpf_data {

	using std::string;
	using std::vector;
	using std::tuple;

	struct peptide_data;
	struct denovo_peptide;
	struct denovo_aminoacid;
	struct FASTA_category;
	struct blastp_data;
	struct proteinconstruct_from_denovo;
	struct category_analysis;
	struct multinomial;

	struct denovo_peptide {
		vector<denovo_aminoacid> v_denovo_aminoacid;
		double localconfidence_average;
	};

	struct denovo_aminoacid {
		char aminoacid;
		double aminoacid_score;
	};

	struct peptide_data {
	public:
		string peptide_withmod;
		string peptide_filtered;
		size_t spectralcount;
		size_t denovo_replicate_count;
		vector<denovo_peptide> v_denovo_peptide_data;
		size_t filesystem_sample_replicate_count = size_t{ 1 };
		bool filesystem_sample_replicate_merged = bool();
		vector<tuple<string, size_t, size_t>> v_filesystem_sample_replicate_data;
	};

	struct FASTA_category {
	public:
		string category_name;
		string category_class;
		string category_type;
		string category_protein;
		string category_species;
	};

	struct blastp_data {
	public:
		string blastp_query;
		string blastp_query_aligned;
		string blastp_subject;
		string blastp_subject_accession;
		string blastp_subject_accession_class;
		size_t blastp_query_alignment_index;
		size_t blastp_subject_alignment_index;
		double blastp_evalue;
		double blastp_evalue_transformed;
		double blastp_parameter_density;
		FASTA_category* p_FASTA_category;
		string query_alignment;
		denovo_peptide denovo_peptide_best_averagelocalconfidence;
		size_t denovo_replicate_count;
	};

	struct proteinconstruct_from_denovo {
		char aminoacid;
		double aminoacid_score;
	};

	struct category_analysis {
	public:
		FASTA_category* p_FASTA_category;
		vector<blastp_data> v_blastp_data_combined_by_category;
		double category_score;
		vector<proteinconstruct_from_denovo> proteinconstruct_from_denovo;
	};

	struct multinomial {
	public:
		vector<string> v_category_name;
		vector<string> v_category_class;
		vector<string> v_element_name;
		vector<vector<double>> v2_frequency;
		vector<double> v_frequency_marginal_sum;
		vector<vector<double>> v2_density;
	};

	vector<FASTA_category> create_FASTA_category(vector<fpf_parse::FASTA_data> par_parse_FASTA) {
		vector<FASTA_category> temp_v_FASTA_category;
		for (const auto itr_parse_FASTA : par_parse_FASTA) {
			const auto find_v_FASTA_element = std::find_if(temp_v_FASTA_category.begin(), temp_v_FASTA_category.end(),
				[itr_parse_FASTA](const FASTA_category& par_FASTA_element) {
				return par_FASTA_element.category_name == itr_parse_FASTA.return_FASTA_name(); });
			if (find_v_FASTA_element == temp_v_FASTA_category.end()) {
				FASTA_category temp_FASTA_category;
				temp_FASTA_category.category_name = itr_parse_FASTA.return_FASTA_name();
				temp_FASTA_category.category_class = itr_parse_FASTA.return_FASTA_class();
				temp_FASTA_category.category_type = itr_parse_FASTA.return_FASTA_type();
				temp_FASTA_category.category_species = itr_parse_FASTA.return_FASTA_species();
				temp_FASTA_category.category_protein = itr_parse_FASTA.return_FASTA_protein();
				temp_v_FASTA_category.push_back(temp_FASTA_category);
			}
			else {
				find_v_FASTA_element->category_protein = find_v_FASTA_element->category_name + itr_parse_FASTA.return_FASTA_protein();
			}
		}
		return temp_v_FASTA_category;
	}

	vector<peptide_data> create_peptide_data(vector<fpf_parse::csv_data> par_parse_csv_peptide_data) {
		vector<peptide_data> temp_v_peptide_data;
		for (const auto itr_parse_csv_peptide_data : par_parse_csv_peptide_data) {
			peptide_data temp_peptide_data = peptide_data();
			denovo_peptide temp_denovo_peptide = denovo_peptide();
			denovo_aminoacid temp_denovo_aminoacid = denovo_aminoacid();
			string temp_peptide_filtered = string();
			size_t sw_peptide_filtered = size_t();
			for (auto j = size_t(); j < itr_parse_csv_peptide_data.csv_peptide.length(); ++j) {
				if (itr_parse_csv_peptide_data.csv_peptide.at(j) == '(' || ((itr_parse_csv_peptide_data.csv_peptide.at(j) == '.') && (itr_parse_csv_peptide_data.csv_peptide.length() > 2))) {
					sw_peptide_filtered = 1;
					if (itr_parse_csv_peptide_data.csv_peptide.at(j + 1) == 's') {
						temp_peptide_filtered.pop_back();
						temp_peptide_filtered += itr_parse_csv_peptide_data.csv_peptide.at(j + 5);
					}
				}
				if (sw_peptide_filtered == 0) {
					temp_peptide_filtered += itr_parse_csv_peptide_data.csv_peptide.at(j);
				}
				if (itr_parse_csv_peptide_data.csv_peptide.at(j) == ')') {
					sw_peptide_filtered = 0;
				}
				if ((itr_parse_csv_peptide_data.csv_peptide.at(j) == '.') && (itr_parse_csv_peptide_data.csv_peptide.length() <= 2) && (sw_peptide_filtered == 0)) {
					temp_peptide_filtered.clear();
				}
			}
			temp_peptide_data.peptide_filtered = temp_peptide_filtered;
			const auto find_peptide_data = std::find_if(temp_v_peptide_data.begin(), temp_v_peptide_data.end(),
				[temp_peptide_filtered](peptide_data par_peptide_data) {
				return par_peptide_data.peptide_filtered == temp_peptide_filtered;
			});
			if (find_peptide_data == temp_v_peptide_data.end()) {
				temp_peptide_data.peptide_withmod = itr_parse_csv_peptide_data.csv_peptide;
				for (size_t i = 0; i < temp_peptide_filtered.size(); ++i) {
					temp_denovo_aminoacid.aminoacid = temp_peptide_filtered[i];
					temp_denovo_aminoacid.aminoacid_score = itr_parse_csv_peptide_data.v_csv_denovo_localconfidence[i];
					temp_denovo_peptide.v_denovo_aminoacid.push_back(temp_denovo_aminoacid);
				}
				temp_denovo_peptide.localconfidence_average = double();
				for (const auto itr_s_denovo_aminoacid : temp_denovo_peptide.v_denovo_aminoacid) {
					temp_denovo_peptide.localconfidence_average += itr_s_denovo_aminoacid.aminoacid_score;
				}
				temp_denovo_peptide.localconfidence_average /= temp_denovo_peptide.v_denovo_aminoacid.size();
				temp_peptide_data.v_denovo_peptide_data.push_back(temp_denovo_peptide);
				++temp_peptide_data.denovo_replicate_count;
				temp_v_peptide_data.push_back(temp_peptide_data);
			}
			else {
				for (size_t i = 0; i < temp_peptide_filtered.size(); ++i) {
					temp_denovo_aminoacid.aminoacid = temp_peptide_filtered[i];
					temp_denovo_aminoacid.aminoacid_score = itr_parse_csv_peptide_data.v_csv_denovo_localconfidence[i];
					temp_denovo_peptide.v_denovo_aminoacid.push_back(temp_denovo_aminoacid);
				}
				temp_denovo_peptide.localconfidence_average = double();
				for (const auto itr_denovo_aminoacid : temp_denovo_peptide.v_denovo_aminoacid) {
					temp_denovo_peptide.localconfidence_average += itr_denovo_aminoacid.aminoacid_score;
				}
				temp_denovo_peptide.localconfidence_average /= temp_denovo_peptide.v_denovo_aminoacid.size();
				find_peptide_data->v_denovo_peptide_data.push_back(temp_denovo_peptide);
				++find_peptide_data->denovo_replicate_count;
			}
		}
		return temp_v_peptide_data;
	}
}

#endif