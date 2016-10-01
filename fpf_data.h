// * * fpf_data.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_DATA
#define	FPF_DATA

#include <algorithm>
#include <cstdlib>
#include <math.h>
#include <sstream>
#include <string>
#include <map>
#include <utility>
#include <vector>

#include "IgFamily.h"
#include "fpf_parse.h"


namespace fpf_data {

	using std::string;
	using std::tuple;
	using std::vector;
	using std::map;
	using std::multimap;
	
	struct denovo_peptide;
	struct denovo_aminoacid;
	struct protein_data;
	struct homology_data;
	struct peptide_data;
	struct protein_analysis;
	struct proteinconstruct_aminoacid;
	struct multinomial;

	struct denovo_peptide {
		vector<denovo_aminoacid> v_denovo_aminoacid;
		double localconfidence_average;
	};

	struct denovo_aminoacid {
		char aminoacid;
		double aminoacid_localconfidence;
	};

	struct peptide_data {
	public:
		size_t key_peptide_data;
		size_t scan_ID;
		string peptide_withmod;
		string peptide_withoutmod;
		string peptide_filtered;
		denovo_peptide denovo_peptide_data;
		denovo_peptide denovo_peptide_data_filtered;
	};
	
	struct peptide_analysis {
		size_t key_peptide_analysis;
		string peptide_filtered;
		size_t replicate_count;
		vector<peptide_data*> v_peptide_data;
		denovo_peptide* p_denovo_peptide_best_by_averagelocalconfidence;
		double v_denovo_peptide_averagescore;
	};

	struct protein_data {
	public:
		size_t key_protein_data;
		string protein_name;
		string protein_class;
		string protein_type;
		string protein_species;
		string protein_protein;
	};

	struct protein_analysis {
	public:
		size_t key_protein_analysis;
		protein_data* p_protein_data;
		vector<homology_data> v_homology_data_combined_by_protein;
		double protein_score;
		double protein_density;
		vector<proteinconstruct_aminoacid> proteinconstruct_from_denovo;
		double proteinconstruct_sequencecoverage;
	};

	struct homology_data {
	public:
		peptide_analysis* p_peptide_analysis;
		protein_data* p_protein_data;
		size_t key_blastp_query;
		string blastp_query;
		string blastp_query_aligned;
		string blastp_subject;
		size_t key_blastp_subject_accession;
		string blastp_subject_accession;
		string blastp_subject_accession_class;
		size_t blastp_query_alignment_index;
		size_t blastp_subject_alignment_index;
		double blastp_score;
		double blastp_score_transformed;
		double blastp_score_transformed_conjugated;
		double blastp_parameter_density;
		double blastp_parameter_density_conjugated;
		double blastp_parameter_score;
		string query_alignment;
		size_t denovo_replicate_count;
		bool is_IG;
	};

	struct proteinconstruct_aminoacid {
		char aminoacid;
		double aminoacid_localconfidence;
		double aminoacid_score;
		double aminoacid_score_conjugated;
		double aminoacid_parameter_score;
	};

	struct multinomial_frequency_type {
		string protein_data;
		double multinomial_frequency;
	};

	struct multinomial {
	public:
		vector<string> v_protein_name;
		vector<string> v_protein_class;
		vector<string> v_element_name;
		vector<vector<double>> v2_frequency;
		vector<double> v_frequency_marginal_sum;
		vector<vector<double>> v2_density;
	};

	map<string, protein_data*> create_v_protein_data_map(vector<protein_data>& par_v_protein_data) {
		map<string, protein_data*> temp_v_protein_data_map{};
		for (auto& itr_v_protein_data : par_v_protein_data) {
			temp_v_protein_data_map[itr_v_protein_data.protein_name] = &itr_v_protein_data;
		}
		return temp_v_protein_data_map;
	}

	vector<protein_data> create_v_protein_data(vector<fpf_parse::FASTA_data> par_parse_FASTA) {
		vector<protein_data> temp_v_protein_data{};
		size_t temp_key_protein_data{};
		for (const auto itr_parse_FASTA : par_parse_FASTA) {
			const auto find_v_FASTA_element = std::find_if(temp_v_protein_data.begin(), temp_v_protein_data.end(),
				[itr_parse_FASTA](const protein_data& par_FASTA_element) {
					return par_FASTA_element.protein_name == itr_parse_FASTA.return_FASTA_name();
				});
			if (find_v_FASTA_element == temp_v_protein_data.end()) {
				protein_data temp_protein_data{
					temp_key_protein_data,
					itr_parse_FASTA.return_FASTA_name(),
					itr_parse_FASTA.return_FASTA_class(),
					itr_parse_FASTA.return_FASTA_type(),
					itr_parse_FASTA.return_FASTA_species(),
					itr_parse_FASTA.return_FASTA_protein()};
				temp_v_protein_data.push_back(temp_protein_data);
				++temp_key_protein_data;
			}
		}
		return temp_v_protein_data;
	}

	multimap<string, peptide_data*> create_v_peptide_data_map(vector<peptide_data>& par_v_peptide_data) {
		multimap<string, peptide_data*> temp_v_peptide_data_map{};
		for (auto& itr_v_peptide_data : par_v_peptide_data) {
			temp_v_peptide_data_map.insert(std::make_pair(itr_v_peptide_data.peptide_withmod, &itr_v_peptide_data));
		}
		return temp_v_peptide_data_map;
	}

	vector<peptide_data> create_v_peptide_data(vector<fpf_parse::csv_data> par_parse_csv_peptide_data) {
		vector<peptide_data> temp_v_peptide_data{};
		size_t temp_key_peptide_data{};
		for (const auto itr_parse_csv_peptide_data : par_parse_csv_peptide_data) {
			peptide_data temp_peptide_data{};
			denovo_peptide temp_denovo_peptide{};
			denovo_aminoacid temp_denovo_aminoacid{};
			string temp_peptide_withoutmod{};
			size_t sw_peptide_filtered{};
			for (auto i = size_t(); i < itr_parse_csv_peptide_data.csv_peptide.length(); ++i) {
				if (itr_parse_csv_peptide_data.csv_peptide.at(i) == '(' || ((itr_parse_csv_peptide_data.csv_peptide.at(i) == '.') && (itr_parse_csv_peptide_data.csv_peptide.length() > 2))) {
					sw_peptide_filtered = 1;
					if (itr_parse_csv_peptide_data.csv_peptide.at(i + 1) == 's') {
						temp_peptide_withoutmod.pop_back();
						temp_peptide_withoutmod += itr_parse_csv_peptide_data.csv_peptide.at(i + 5);
					}
				}
				if (sw_peptide_filtered == 0) {
					temp_peptide_withoutmod += itr_parse_csv_peptide_data.csv_peptide.at(i);
				}
				if (itr_parse_csv_peptide_data.csv_peptide.at(i) == ')') {
					sw_peptide_filtered = 0;
				}
				if ((itr_parse_csv_peptide_data.csv_peptide.at(i) == '.') && (itr_parse_csv_peptide_data.csv_peptide.length() <= 2) && (sw_peptide_filtered == 0)) {
					temp_peptide_withoutmod.clear();
				}
			}
			temp_peptide_data.peptide_withmod = itr_parse_csv_peptide_data.csv_peptide;
			temp_peptide_data.peptide_withoutmod = temp_peptide_withoutmod;
			temp_peptide_data.peptide_filtered = temp_peptide_withoutmod;
			for (auto j = 0; j < temp_peptide_withoutmod.size(); ++j) {
				temp_denovo_aminoacid.aminoacid = temp_peptide_withoutmod[j];
				temp_denovo_aminoacid.aminoacid_localconfidence = itr_parse_csv_peptide_data.v_csv_denovo_localconfidence[j];
				temp_denovo_peptide.v_denovo_aminoacid.push_back(temp_denovo_aminoacid);
			}
			temp_denovo_peptide.localconfidence_average = double();
			for (const auto itr_denovo_aminoacid : temp_denovo_peptide.v_denovo_aminoacid) {
				temp_denovo_peptide.localconfidence_average += itr_denovo_aminoacid.aminoacid_localconfidence;
			}
			temp_denovo_peptide.localconfidence_average /= temp_denovo_peptide.v_denovo_aminoacid.size();
			temp_peptide_data.denovo_peptide_data = temp_denovo_peptide;
			temp_peptide_data.denovo_peptide_data_filtered = temp_denovo_peptide;
			if (IgFamily::BLASTP_BY_SELECTED_PEPTIDE) {
				string temp_peptide_selected{};
				vector<string> temp_v_peptide_selected{};
				denovo_peptide temp_denovo_peptide_2{};
				vector<denovo_peptide> temp_v_denovo_peptide{};
				vector<double> v_moving_value{};
				for (auto j = 0; j < temp_peptide_data.denovo_peptide_data_filtered.v_denovo_aminoacid.size(); ++j) {
					v_moving_value.push_back(temp_peptide_data.denovo_peptide_data_filtered.v_denovo_aminoacid[j].aminoacid_localconfidence);
					if (v_moving_value.size() == 4) {
						v_moving_value.erase(v_moving_value.begin());
					}
					double temp_moving_average{};
					for (const auto& itr_v_moving_value : v_moving_value) {
						temp_moving_average += itr_v_moving_value;
					}
					temp_moving_average /= v_moving_value.size();
					if ((temp_moving_average < IgFamily::DENOVO_LOCAL_CONFIDENCE_MOVING_AVERAGE_THRESHOLD)
						&& (temp_peptide_data.denovo_peptide_data_filtered.v_denovo_aminoacid[j].aminoacid_localconfidence < IgFamily::DENOVO_LOCAL_CONFIDENCE_MOVING_AVERAGE_THRESHOLD)) {
						if (temp_peptide_selected.size() > IgFamily::DENOVO_PEPTIDE_SIZE_THRESHOLD) {
							temp_peptide_selected.pop_back();
							temp_denovo_peptide_2.v_denovo_aminoacid.pop_back();
							temp_v_peptide_selected.push_back(temp_peptide_selected);
							temp_v_denovo_peptide.push_back(temp_denovo_peptide_2);
						}
						temp_peptide_selected.clear();
						temp_denovo_peptide_2.v_denovo_aminoacid.clear();
					}
					else {
						temp_peptide_selected += temp_peptide_data.denovo_peptide_data_filtered.v_denovo_aminoacid[j].aminoacid;
						temp_denovo_peptide_2.v_denovo_aminoacid.push_back(temp_peptide_data.denovo_peptide_data_filtered.v_denovo_aminoacid[j]);
					}
					if (((j + 1) == temp_peptide_data.denovo_peptide_data_filtered.v_denovo_aminoacid.size())
						&& (temp_peptide_selected.size() > IgFamily::DENOVO_PEPTIDE_SIZE_THRESHOLD)) {
						temp_v_peptide_selected.push_back(temp_peptide_selected);
						temp_v_denovo_peptide.push_back(temp_denovo_peptide_2);
					}
				}
				temp_peptide_selected = string();
				for (auto j = 0; j < temp_v_peptide_selected.size(); ++j) {
					if (temp_v_peptide_selected[j].size() > temp_peptide_selected.size()) {
						temp_peptide_selected = temp_v_peptide_selected[j];
						temp_v_denovo_peptide[j].localconfidence_average = double();
						for (const auto itr_denovo_aminoacid : temp_v_denovo_peptide[j].v_denovo_aminoacid) {
							temp_v_denovo_peptide[j].localconfidence_average += itr_denovo_aminoacid.aminoacid_localconfidence;
						}
						temp_v_denovo_peptide[j].localconfidence_average /= temp_v_denovo_peptide[j].v_denovo_aminoacid.size();
						temp_peptide_data.denovo_peptide_data_filtered = temp_v_denovo_peptide[j];
					}
				}
				temp_peptide_data.peptide_filtered = temp_peptide_selected;
			}
			if ((temp_denovo_peptide.localconfidence_average > IgFamily::DENOVO_PEPTIDE_CONFIDENCE_THRESHOLD) && (temp_peptide_data.peptide_filtered != "")) {
				temp_peptide_data.key_peptide_data = temp_key_peptide_data;
				temp_peptide_data.scan_ID = std::stoi(itr_parse_csv_peptide_data.csv_scan_ID);
				temp_v_peptide_data.push_back(temp_peptide_data);
				++temp_key_peptide_data;
			}
		}
		return temp_v_peptide_data;
	}
}

#endif