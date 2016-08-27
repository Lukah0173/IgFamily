// * * fpf_data.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_DATA
#define	FPF_DATA

#include <cstdlib>
#include <vector>
#include <sstream>
#include <algorithm>
#include <string>
#include <utility>
#include <math.h>

#include "IgFamily.h"
#include "fpf_parse.h"


namespace fpf_data {

	using std::string;
	using std::tuple;
	using std::vector;

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
		size_t scan_ID;
		string peptide_withmod;
		string peptide_filtered;
		denovo_peptide denovo_peptide_data;
	};
	
	struct peptide_analysis {
		string peptide_filtered;
		size_t replicate_count;
		vector<peptide_data*> v_peptide_data; 
		denovo_peptide* p_denovo_peptide_best_by_averagelocalconfidence;
		double v_denovo_peptide_averagescore;
	};

	struct protein_data {
	public:
		string protein_name;
		string protein_class;
		string protein_type;
		string protein_species;
		string protein_protein;
	};

	struct protein_analysis {
	public:
		protein_data* p_protein_data;
		vector<homology_data> v_homology_data_combined_by_protein;
		double protein_score;
		vector<proteinconstruct_aminoacid> proteinconstruct_from_denovo;
		double proteinconstruct_sequencecoverage;
	};

	struct homology_data {
	public:
		peptide_analysis* p_peptide_analysis;
		protein_data* p_protein_data;
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
		double blastp_parameter_score;
		string query_alignment;
		size_t denovo_replicate_count;
	};

	struct proteinconstruct_aminoacid {
		char aminoacid;
		double aminoacid_localconfidence;
		double aminoacid_evalue_transformed;
		double aminoacid_parameter_score;
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

	vector<protein_data> create_protein_data(vector<fpf_parse::FASTA_data> par_parse_FASTA) {
		vector<protein_data> temp_v_protein_data{};
		for (const auto itr_parse_FASTA : par_parse_FASTA) {
			const auto find_v_FASTA_element = std::find_if(temp_v_protein_data.begin(), temp_v_protein_data.end(),
				[itr_parse_FASTA](const protein_data& par_FASTA_element) {
				return par_FASTA_element.protein_name == itr_parse_FASTA.return_FASTA_name(); });
			if (find_v_FASTA_element == temp_v_protein_data.end()) {
				protein_data temp_protein_data{
					itr_parse_FASTA.return_FASTA_name(),
					itr_parse_FASTA.return_FASTA_class(),
					itr_parse_FASTA.return_FASTA_type(),
					itr_parse_FASTA.return_FASTA_species(),
					itr_parse_FASTA.return_protein_data()};
				temp_v_protein_data.push_back(temp_protein_data);
			}
		}
		return temp_v_protein_data;
	}

	vector<peptide_data> create_peptide_data(vector<fpf_parse::csv_data> par_parse_csv_peptide_data) {
		vector<peptide_data> temp_v_peptide_data{};
		for (const auto itr_parse_csv_peptide_data : par_parse_csv_peptide_data) {
			peptide_data temp_peptide_data{};
			denovo_peptide temp_denovo_peptide{};
			denovo_aminoacid temp_denovo_aminoacid{};
			string temp_peptide_filtered{};
			size_t sw_peptide_filtered{};
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
			temp_peptide_data.peptide_withmod = itr_parse_csv_peptide_data.csv_peptide;
			temp_peptide_data.peptide_filtered = temp_peptide_filtered;
			for (size_t i = 0; i < temp_peptide_filtered.size(); ++i) {
				temp_denovo_aminoacid.aminoacid = temp_peptide_filtered[i];
				temp_denovo_aminoacid.aminoacid_localconfidence = itr_parse_csv_peptide_data.v_csv_denovo_localconfidence[i];
				temp_denovo_peptide.v_denovo_aminoacid.push_back(temp_denovo_aminoacid);
			}
			temp_denovo_peptide.localconfidence_average = double();
			for (const auto itr_s_denovo_aminoacid : temp_denovo_peptide.v_denovo_aminoacid) {
				temp_denovo_peptide.localconfidence_average += itr_s_denovo_aminoacid.aminoacid_localconfidence;
			}
			temp_denovo_peptide.localconfidence_average /= temp_denovo_peptide.v_denovo_aminoacid.size();
			temp_peptide_data.denovo_peptide_data = temp_denovo_peptide;
			if (temp_denovo_peptide.localconfidence_average > DENOVO_PEPTIDE_CONFIDENCE_THRESHOLD) {
				temp_v_peptide_data.push_back(temp_peptide_data);
			}
		}
		return temp_v_peptide_data;
	}

	vector<peptide_analysis> create_peptide_analysis(vector<peptide_data>& par_v_peptide_data) {
		vector<peptide_analysis> temp_v_peptide_analysis{};
		for (auto& itr_v_peptide_data : par_v_peptide_data) {
			peptide_analysis temp_peptide_analysis{};
			auto& find_peptide_analysis = std::find_if(temp_v_peptide_analysis.begin(), temp_v_peptide_analysis.end(),
				[itr_v_peptide_data](peptide_analysis par_peptide_analysis) {
				return par_peptide_analysis.peptide_filtered == itr_v_peptide_data.peptide_filtered;
			});
			if (find_peptide_analysis == temp_v_peptide_analysis.end()) {
				temp_peptide_analysis.peptide_filtered = itr_v_peptide_data.peptide_filtered;
				++temp_peptide_analysis.replicate_count;
				temp_peptide_analysis.v_peptide_data.push_back(&itr_v_peptide_data);
				temp_peptide_analysis.v_denovo_peptide_averagescore = itr_v_peptide_data.denovo_peptide_data.localconfidence_average;
				temp_v_peptide_analysis.push_back(temp_peptide_analysis);
			}
			else {
				++find_peptide_analysis->replicate_count;
				find_peptide_analysis->v_peptide_data.push_back(&itr_v_peptide_data);
				find_peptide_analysis->v_denovo_peptide_averagescore
					= ((find_peptide_analysis->v_denovo_peptide_averagescore * (find_peptide_analysis->replicate_count - 1)) + itr_v_peptide_data.denovo_peptide_data.localconfidence_average) / find_peptide_analysis->replicate_count;
			}
		}
		return temp_v_peptide_analysis;
	}
}

#endif