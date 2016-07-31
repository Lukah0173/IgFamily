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

	struct s_multinomial_element_data;
	struct s_peptide_data;
	struct s_denovo_peptide;
	struct s_blastp;
	struct s_mnom;
	struct s_proteinconstruct_from_denovo;
	struct s_report;

	typedef std::string string_type;
	typedef size_t size_type;

	struct s_multinomial_element_data {
	public:
		s_multinomial_element_data() {
			d_score = { 1 };
		};

		~s_multinomial_element_data() {
		};

		string_type str_multinomial_element_name;
		string_type str_protein;
		string_type str_species;
		string_type str_multinomial_element_class;	
		std::vector<s_peptide_data*> v_s_peptide_data;
		std::vector<s_peptide_data*> v_s_peptide_data_distinct_filtered;
		string_type str_alignment;
		size_type st_totalspectralcount;
		double d_coverage;
		double d_score;
		std::vector<s_multinomial_element_data*> v_s_multinomial_element_polyassociation;
	};

	struct s_denovo_peptide {
		char ch_aminoacid;
		double d_denovo_localconfidence;
	};

	struct s_peptide_data {
	public:
		s_peptide_data() {
		};

		~s_peptide_data() {
		};

		string_type str_peptide;
		size_type st_spectralcount;
		size_type st_IgP;
		std::vector<s_denovo_peptide> v_s_denovo_peptide;
		double d_denovo_peptide_localconfidence_average;
		bool b_replicate_merged = bool();
		size_type st_filesystem_replicate = size_type{ 1 };
		std::vector<std::tuple<string_type, size_type, size_type>> v_p_replicate_data;
		std::vector<string_type> v_str_peptideassociation;
		std::vector<string_type> v_str_peptideassociation_distinct;
		std::vector<std::pair<s_multinomial_element_data*, double>> v_p_peptideassociation;
		std::vector<std::pair<s_multinomial_element_data*, double>> v_p_peptideassociation_distinct;
		s_denovo_peptide s_blastp_denovo_peptide;
	};

	struct s_blastp {
	public:
		string_type str_blastp_query;
		string_type str_blastp_query_aligned;
		string_type str_blastp_subject;
		string_type str_blastp_subject_accession;
		size_type st_blastp_query_alignment_index;
		size_type st_blastp_subject_alignment_index;
		string_type str_protein;
		string_type str_blastp_query_alignment;
		double d_blastp_evalue;
		double d_blastp_par_prop;
		s_denovo_peptide s_blastp_denovo_peptide;
	};

	struct s_proteinconstruct_from_denovo {
		char ch_aminoacid;
		double d_score;
	};

	struct s_mnom {
	public:
		string_type str_mnom_comp;
		double d_mnom_value;
	};

	struct s_report {
	public:
		std::vector<s_blastp> v_s_blastp;
		string_type str_protein_accession;
		string_type str_protein;
		std::vector<s_proteinconstruct_from_denovo> s_proteinconstruct_from_denovo;
		double d_score;
		std::vector<double> v_d_aminoacid_scores;
	};

	std::vector<s_multinomial_element_data> create_v_s_multinomial_element_data(std::vector<fpf_parse::c_parse_FASTA> par_v_c_parse_FASTA) {

		// Vector reallocation will invalidate all pointers to s_multinomial_element_data references.
		// Do not reallocate after return!

		std::vector<s_multinomial_element_data> con_v_s_multinomial_element_data;
		s_multinomial_element_data con_c_analysis;
		size_type con_st_ID = size_type();
		for (auto itr_par_v_c_parse_FASTA : par_v_c_parse_FASTA) {
			++con_st_ID;
			auto find_v_s_multinomial_element_data = std::find_if(con_v_s_multinomial_element_data.begin(), con_v_s_multinomial_element_data.end(),
				[itr_par_v_c_parse_FASTA](s_multinomial_element_data par_s_multinomial_element_data) {
				return par_s_multinomial_element_data.str_multinomial_element_name == itr_par_v_c_parse_FASTA.return_str_parse_FASTA_genefamily(); });
			if (find_v_s_multinomial_element_data == con_v_s_multinomial_element_data.end()) {
				con_c_analysis.str_multinomial_element_name = itr_par_v_c_parse_FASTA.return_str_parse_FASTA_genefamily();
				con_c_analysis.str_multinomial_element_class = itr_par_v_c_parse_FASTA.return_str_parse_FASTA_genefamily_class();
				con_c_analysis.str_species = itr_par_v_c_parse_FASTA.return_str_parse_FASTA_species();
				con_c_analysis.str_protein = itr_par_v_c_parse_FASTA.return_str_parse_FASTA_protein();
				con_v_s_multinomial_element_data.push_back(con_c_analysis);
			}
			else {
				find_v_s_multinomial_element_data->str_protein = find_v_s_multinomial_element_data->str_multinomial_element_name + itr_par_v_c_parse_FASTA.return_str_parse_FASTA_protein();
			}
		}
		return con_v_s_multinomial_element_data;
	}

	std::vector<s_multinomial_element_data> create_v_s_multinomial_element_data_distinct(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data) {
		std::vector<s_multinomial_element_data> con_v_s_multinomial_element_data_distinct;
		std::vector<string_type> con_v_str_multinomial_element_name_distinct;
		s_multinomial_element_data con_c_analysis_distinct;
		for (auto itr_v_s_multinomial_element_data_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data_data) {
			string_type con_str_multinomial_element_name = itr_v_s_multinomial_element_data_data->str_multinomial_element_name;
			size_type sw_peptideassociation_distinct = size_type();
			string_type con_str_multinomial_element_name_distinct = string_type();
			string_type con_str_multinomial_element_name_distinct_class = string_type();
			for (auto itr_str_multinomial_element_name = size_type(); itr_str_multinomial_element_name < con_str_multinomial_element_name.length(); ++itr_str_multinomial_element_name) {
				if (con_str_multinomial_element_name.at(itr_str_multinomial_element_name) == '*') {
					sw_peptideassociation_distinct = 1;
				}
				if (sw_peptideassociation_distinct == 0) {
					con_str_multinomial_element_name_distinct += con_str_multinomial_element_name.at(itr_str_multinomial_element_name);
				}
			}
			con_str_multinomial_element_name_distinct_class = itr_v_s_multinomial_element_data_data->str_multinomial_element_class;
			if (std::find(con_v_str_multinomial_element_name_distinct.begin(), con_v_str_multinomial_element_name_distinct.end(), con_str_multinomial_element_name_distinct) == con_v_str_multinomial_element_name_distinct.end()) {
				con_v_str_multinomial_element_name_distinct.push_back(con_str_multinomial_element_name_distinct);
				con_c_analysis_distinct.str_multinomial_element_name = con_str_multinomial_element_name_distinct;
				con_c_analysis_distinct.str_multinomial_element_class = con_str_multinomial_element_name_distinct_class;
				con_c_analysis_distinct.str_species = itr_v_s_multinomial_element_data_data->str_species;
				con_v_s_multinomial_element_data_distinct.push_back(con_c_analysis_distinct);
			}
			for (auto itr_v_s_multinomial_element_data_distinct = con_v_s_multinomial_element_data_distinct.begin(); itr_v_s_multinomial_element_data_distinct != con_v_s_multinomial_element_data_distinct.end(); ++itr_v_s_multinomial_element_data_distinct) {
				if (itr_v_s_multinomial_element_data_distinct->str_multinomial_element_name == con_str_multinomial_element_name_distinct) {
					itr_v_s_multinomial_element_data_distinct->v_s_multinomial_element_polyassociation.push_back(&(*itr_v_s_multinomial_element_data_data));
				}
			}
			sw_peptideassociation_distinct = 0;
			con_str_multinomial_element_name_distinct.clear();
		}
		return con_v_s_multinomial_element_data_distinct;
	}

	void update_v_s_multinomial_element_distinctpolymorphism_data(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data_distinct) {
		s_multinomial_element_data con_c_analysis_distinct_update;
		for (auto itr_v_s_multinomial_element_data_distinct = par_v_s_multinomial_element_data_distinct.begin(); itr_v_s_multinomial_element_data_distinct != par_v_s_multinomial_element_data_distinct.end(); ++itr_v_s_multinomial_element_data_distinct) {
			std::vector<s_multinomial_element_data*> con_v_s_multinomial_element_data_distinct_update = itr_v_s_multinomial_element_data_distinct->v_s_multinomial_element_polyassociation;
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
			con_c_analysis_distinct_update.v_s_multinomial_element_polyassociation = itr_v_s_multinomial_element_data_distinct->v_s_multinomial_element_polyassociation;
			con_c_analysis_distinct_update.str_multinomial_element_name = itr_v_s_multinomial_element_data_distinct->str_multinomial_element_name;
			con_c_analysis_distinct_update.str_multinomial_element_class = itr_v_s_multinomial_element_data_distinct->str_multinomial_element_class;
			con_c_analysis_distinct_update.str_species = itr_v_s_multinomial_element_data_distinct->str_species;
			*itr_v_s_multinomial_element_data_distinct = con_c_analysis_distinct_update;
		}
	}

	std::vector<s_peptide_data> create_v_s_peptide_data(std::vector<fpf_parse::s_parse_peptides_csv> par_c_parse_csv_peptide_data) {
		std::vector<s_peptide_data> con_v_s_peptide_data;
		for (auto itr_c_parse_csv_peptide_data = par_c_parse_csv_peptide_data.begin(); itr_c_parse_csv_peptide_data != par_c_parse_csv_peptide_data.end(); ++itr_c_parse_csv_peptide_data) {
			s_peptide_data con_s_peptide_data = s_peptide_data();
			s_denovo_peptide con_s_denovo_peptide = s_denovo_peptide();
			size_type ss_st_csv_IgP = size_type();
			std::istringstream(itr_c_parse_csv_peptide_data->str_parse_csv_IgP) >> ss_st_csv_IgP;
			con_s_peptide_data.str_peptide = itr_c_parse_csv_peptide_data->str_parse_peptides_csv_peptide;
			std::istringstream ss_spectralcount(itr_c_parse_csv_peptide_data->str_parse_peptides_csv_spectralcount);
			size_type ss_st_spectralcount;
			ss_spectralcount >> ss_st_spectralcount;
			con_s_peptide_data.st_spectralcount = ss_st_spectralcount;
			for (size_type i = 0; i < itr_c_parse_csv_peptide_data->str_parse_peptides_csv_peptide.size(); ++i) {
				con_s_denovo_peptide.ch_aminoacid = itr_c_parse_csv_peptide_data->str_parse_peptides_csv_peptide[i];
				con_s_denovo_peptide.d_denovo_localconfidence = itr_c_parse_csv_peptide_data->v_d_denovo_localconfidence[i];
				con_s_peptide_data.v_s_denovo_peptide.push_back(con_s_denovo_peptide);
			}
			con_s_peptide_data.d_denovo_peptide_localconfidence_average = double();
			for (auto itr_s_denovo_peptide : con_s_peptide_data.v_s_denovo_peptide) {
				con_s_peptide_data.d_denovo_peptide_localconfidence_average += itr_s_denovo_peptide.d_denovo_localconfidence;
			}
			con_s_peptide_data.d_denovo_peptide_localconfidence_average /= con_s_peptide_data.v_s_denovo_peptide.size();
			con_v_s_peptide_data.push_back(con_s_peptide_data);
		}
		return con_v_s_peptide_data;
	}

	std::vector<s_peptide_data> create_v_s_peptide_data_filtered(std::vector<s_peptide_data> par_v_s_peptide_data) {
		std::vector<s_peptide_data> v_s_peptide_data_filtered;
		size_type sw_peptide_filtered = size_type();
		string_type con_str_peptide_filtered = string_type();
		for (auto itr_par_v_data = par_v_s_peptide_data.begin(); itr_par_v_data != par_v_s_peptide_data.end(); ++itr_par_v_data) {
			for (auto j = size_type(); j < itr_par_v_data->str_peptide.length(); ++j) {
				if (itr_par_v_data->str_peptide.at(j) == '(' || ((itr_par_v_data->str_peptide.at(j) == '.') && (itr_par_v_data->str_peptide.length() > 2))) {
					sw_peptide_filtered = 1;
					if (itr_par_v_data->str_peptide.at(j + 1) == 's') {
						con_str_peptide_filtered.pop_back();
						con_str_peptide_filtered += itr_par_v_data->str_peptide.at(j + 5);
					}
				}
				if (sw_peptide_filtered == 0) {
					con_str_peptide_filtered += itr_par_v_data->str_peptide.at(j);
				}
				if (itr_par_v_data->str_peptide.at(j) == ')') {
					sw_peptide_filtered = 0;
				}
				if ((itr_par_v_data->str_peptide.at(j) == '.') && (itr_par_v_data->str_peptide.length() <= 2) && (sw_peptide_filtered == 0)) {
					con_str_peptide_filtered.clear();
				}
			}
			s_peptide_data con_s_peptide_data_filtered;
			con_s_peptide_data_filtered.str_peptide = con_str_peptide_filtered;
			con_s_peptide_data_filtered.st_spectralcount = itr_par_v_data->st_spectralcount;
			con_s_peptide_data_filtered.st_IgP = itr_par_v_data->st_IgP;
			con_s_peptide_data_filtered.v_s_denovo_peptide = itr_par_v_data->v_s_denovo_peptide;
			con_s_peptide_data_filtered.d_denovo_peptide_localconfidence_average = itr_par_v_data->d_denovo_peptide_localconfidence_average;
			con_s_peptide_data_filtered.st_filesystem_replicate = itr_par_v_data->st_filesystem_replicate;
			con_s_peptide_data_filtered.v_p_replicate_data = itr_par_v_data->v_p_replicate_data;
			v_s_peptide_data_filtered.push_back(con_s_peptide_data_filtered);
			con_str_peptide_filtered.clear();
			sw_peptide_filtered = 0;
		}

		if (IgFamily::DEBUG_MODE == 1) {
			string_type output_v_s_peptide_data_filtered = "z_" + IgFamily::output + "_v_s_peptide_data_filtered.txt";
			std::ofstream fout_v_s_peptide_data_filtered;
			fout_v_s_peptide_data_filtered.open(output_v_s_peptide_data_filtered);
			fout_v_s_peptide_data_filtered << "-- IgFamily " << IgFamily::version << " --\n\n\n";
			fout_v_s_peptide_data_filtered << "Input file: " << IgFamily::INPUT_CSV << "\n\n\n";
			for (auto itr_v_s_peptide_data_filtered = v_s_peptide_data_filtered.begin(); itr_v_s_peptide_data_filtered != v_s_peptide_data_filtered.end(); ++itr_v_s_peptide_data_filtered) {
				fout_v_s_peptide_data_filtered << "\n" << itr_v_s_peptide_data_filtered->str_peptide;
				fout_v_s_peptide_data_filtered << " " << itr_v_s_peptide_data_filtered->st_spectralcount;
			}
		}
		return v_s_peptide_data_filtered;
	}

	std::vector<s_peptide_data> create_v_s_peptide_data_distinct(std::vector<s_peptide_data> par_v_s_peptide_data) {
		std::vector<s_peptide_data> con_v_s_peptide_data_distinct;
		for (auto itr_par_v_data = par_v_s_peptide_data.begin(); itr_par_v_data != par_v_s_peptide_data.end(); ++itr_par_v_data) {
			s_peptide_data con_s_peptide_data_distinct;
			if (con_v_s_peptide_data_distinct.size() == 0) {
				con_s_peptide_data_distinct.str_peptide = itr_par_v_data->str_peptide;
				con_s_peptide_data_distinct.st_spectralcount = itr_par_v_data->st_spectralcount;
				con_s_peptide_data_distinct.st_IgP = itr_par_v_data->st_IgP;
				con_s_peptide_data_distinct.v_s_denovo_peptide = itr_par_v_data->v_s_denovo_peptide;
				con_s_peptide_data_distinct.d_denovo_peptide_localconfidence_average = itr_par_v_data->d_denovo_peptide_localconfidence_average;
				con_s_peptide_data_distinct.st_filesystem_replicate = itr_par_v_data->st_filesystem_replicate;
				con_s_peptide_data_distinct.v_p_replicate_data = itr_par_v_data->v_p_replicate_data;
				con_v_s_peptide_data_distinct.push_back(con_s_peptide_data_distinct);
			}
			else {
				size_type st_data_distinct = size_type();
				for (auto itr_v_data_distinct = con_v_s_peptide_data_distinct.begin(); itr_v_data_distinct != con_v_s_peptide_data_distinct.end(); ++itr_v_data_distinct) {
					++st_data_distinct;
					if (itr_v_data_distinct->str_peptide == itr_par_v_data->str_peptide) {
						break;
					}
					if (st_data_distinct == con_v_s_peptide_data_distinct.size()) {
						con_s_peptide_data_distinct.str_peptide = (itr_par_v_data->str_peptide);
						con_s_peptide_data_distinct.st_spectralcount = (itr_par_v_data->st_spectralcount);
						con_s_peptide_data_distinct.st_IgP = itr_par_v_data->st_IgP;
						con_s_peptide_data_distinct.v_s_denovo_peptide = itr_par_v_data->v_s_denovo_peptide;
						con_s_peptide_data_distinct.d_denovo_peptide_localconfidence_average = itr_par_v_data->d_denovo_peptide_localconfidence_average;
						con_s_peptide_data_distinct.st_filesystem_replicate = itr_par_v_data->st_filesystem_replicate;
						con_s_peptide_data_distinct.v_p_replicate_data = itr_par_v_data->v_p_replicate_data;
						con_v_s_peptide_data_distinct.push_back(con_s_peptide_data_distinct);
						break;
					}
				}
			}
		}

		if (IgFamily::DEBUG_MODE == 1) {
			string_type output_v_s_peptide_data_distinct = "z_" + IgFamily::output + "_v_s_peptide_data_distinct.txt";
			std::ofstream fout_v_c_data_distinct;
			fout_v_c_data_distinct.open(output_v_s_peptide_data_distinct);
			fout_v_c_data_distinct << "-- IgFamily " << IgFamily::version << " --\n\n\n";
			fout_v_c_data_distinct << "Input file: " << IgFamily::INPUT_CSV << "\n\n\n";
			for (auto itr_par_v_data_distinct = con_v_s_peptide_data_distinct.begin(); itr_par_v_data_distinct != con_v_s_peptide_data_distinct.end(); ++itr_par_v_data_distinct) {
				fout_v_c_data_distinct << "\n" << itr_par_v_data_distinct->str_peptide;
				fout_v_c_data_distinct << " " << itr_par_v_data_distinct->st_spectralcount;
			}
		}

		return con_v_s_peptide_data_distinct;
	}

	std::vector<s_peptide_data> create_v_s_peptide_data_filtered_distinct(std::vector<s_peptide_data> par_v_s_peptide_data_filtered) {
		std::vector<s_peptide_data> v_s_peptide_data_filtered_distinct;
		for (auto itr_par_v_data_filtered = par_v_s_peptide_data_filtered.begin(); itr_par_v_data_filtered != par_v_s_peptide_data_filtered.end(); ++itr_par_v_data_filtered) {
			s_peptide_data con_s_peptide_data_filtered_distinct;
			if (v_s_peptide_data_filtered_distinct.size() == 0) {
				con_s_peptide_data_filtered_distinct.str_peptide = itr_par_v_data_filtered->str_peptide;
				con_s_peptide_data_filtered_distinct.st_spectralcount = itr_par_v_data_filtered->st_spectralcount;
				con_s_peptide_data_filtered_distinct.st_IgP = itr_par_v_data_filtered->st_IgP;
				con_s_peptide_data_filtered_distinct.v_s_denovo_peptide = itr_par_v_data_filtered->v_s_denovo_peptide;
				con_s_peptide_data_filtered_distinct.d_denovo_peptide_localconfidence_average = itr_par_v_data_filtered->d_denovo_peptide_localconfidence_average;
				con_s_peptide_data_filtered_distinct.st_filesystem_replicate = itr_par_v_data_filtered->st_filesystem_replicate;
				con_s_peptide_data_filtered_distinct.v_p_replicate_data = itr_par_v_data_filtered->v_p_replicate_data;
				v_s_peptide_data_filtered_distinct.push_back(con_s_peptide_data_filtered_distinct);
			}
			else {
				size_type st_data_filtered_distinct = size_type();
				for (auto itr_v_data_filtered_distinct = v_s_peptide_data_filtered_distinct.begin(); itr_v_data_filtered_distinct != v_s_peptide_data_filtered_distinct.end(); ++itr_v_data_filtered_distinct) {
					++st_data_filtered_distinct;
					if (itr_v_data_filtered_distinct->str_peptide == itr_par_v_data_filtered->str_peptide) {
						break;
					}
					if (st_data_filtered_distinct == v_s_peptide_data_filtered_distinct.size()) {
						con_s_peptide_data_filtered_distinct.str_peptide = itr_par_v_data_filtered->str_peptide;
						con_s_peptide_data_filtered_distinct.st_spectralcount = itr_par_v_data_filtered->st_spectralcount;
						con_s_peptide_data_filtered_distinct.st_IgP = itr_par_v_data_filtered->st_IgP;
						con_s_peptide_data_filtered_distinct.v_s_denovo_peptide = itr_par_v_data_filtered->v_s_denovo_peptide;
						con_s_peptide_data_filtered_distinct.d_denovo_peptide_localconfidence_average = itr_par_v_data_filtered->d_denovo_peptide_localconfidence_average;
						con_s_peptide_data_filtered_distinct.st_filesystem_replicate = itr_par_v_data_filtered->st_filesystem_replicate;
						con_s_peptide_data_filtered_distinct.v_p_replicate_data = itr_par_v_data_filtered->v_p_replicate_data;
						v_s_peptide_data_filtered_distinct.push_back(con_s_peptide_data_filtered_distinct);
						break;
					}
				}
			}
		}

		if (IgFamily::DEBUG_MODE == 1) {
			string_type output_v_s_peptide_data_filtered_distinct = "z_" + IgFamily::output + "_v_c_data_filtered_distinct.txt";
			std::ofstream fout_v_s_peptide_data_filtered_distinct;
			fout_v_s_peptide_data_filtered_distinct.open(output_v_s_peptide_data_filtered_distinct);
			fout_v_s_peptide_data_filtered_distinct << "-- IgFamily " << IgFamily::version << " --\n\n\n";
			fout_v_s_peptide_data_filtered_distinct << "Input file: " << IgFamily::INPUT_CSV << "\n\n\n";
			for (auto itr_v_s_peptide_data_filtered_distinct = v_s_peptide_data_filtered_distinct.begin(); itr_v_s_peptide_data_filtered_distinct != v_s_peptide_data_filtered_distinct.end(); ++itr_v_s_peptide_data_filtered_distinct) {
				fout_v_s_peptide_data_filtered_distinct << "\n" << itr_v_s_peptide_data_filtered_distinct->str_peptide;
				fout_v_s_peptide_data_filtered_distinct << " " << itr_v_s_peptide_data_filtered_distinct->st_spectralcount;
			}
		}

		return v_s_peptide_data_filtered_distinct;
	}

	void create_v_s_multinomial_element_data_v_s_peptide_data(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data, std::vector<s_peptide_data> par_v_s_peptide_data) {
		string_type con_str_peptide = string_type();
		std::vector<s_peptide_data*> con_v_s_peptide_data;
		size_type sw_peptide_filtered = size_type();
		for (auto itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data) {
			for (auto itr_par_v_data = par_v_s_peptide_data.begin(); itr_par_v_data != par_v_s_peptide_data.end(); ++itr_par_v_data) {
				for (auto j = size_type(); j < itr_par_v_data->str_peptide.length(); ++j) {
					if (itr_par_v_data->str_peptide.at(j) == '(' || ((itr_par_v_data->str_peptide.at(j) == '.') && (itr_par_v_data->str_peptide.length() > 2))) {
						sw_peptide_filtered = 1;
						if (itr_par_v_data->str_peptide.at(j + 1) == 's') {
							con_str_peptide.pop_back();
							con_str_peptide += itr_par_v_data->str_peptide.at(j + 5);
						}
					}
					if (sw_peptide_filtered == 0) {
						con_str_peptide += itr_par_v_data->str_peptide.at(j);
					}
					if (itr_par_v_data->str_peptide.at(j) == ')') {
						sw_peptide_filtered = 0;
					}
					if ((itr_par_v_data->str_peptide.at(j) == '.') && (itr_par_v_data->str_peptide.length() <= 2) && (sw_peptide_filtered == 0)) {
						con_str_peptide.clear();
					}
				}
				if (itr_v_s_multinomial_element_data->str_protein.find(con_str_peptide) != std::string::npos) {
					s_peptide_data* con_s_peptide_data = new s_peptide_data;
					con_s_peptide_data->str_peptide = itr_par_v_data->str_peptide;
					con_s_peptide_data->st_spectralcount = itr_par_v_data->st_spectralcount;
					con_s_peptide_data->st_IgP = itr_par_v_data->st_IgP;
					con_s_peptide_data->v_s_denovo_peptide = itr_par_v_data->v_s_denovo_peptide;
					con_s_peptide_data->v_str_peptideassociation = itr_par_v_data->v_str_peptideassociation;
					con_s_peptide_data->v_str_peptideassociation_distinct = itr_par_v_data->v_str_peptideassociation_distinct;
					con_s_peptide_data->v_p_peptideassociation = itr_par_v_data->v_p_peptideassociation;
					con_s_peptide_data->v_str_peptideassociation_distinct = itr_par_v_data->v_str_peptideassociation_distinct;
					con_s_peptide_data->st_filesystem_replicate = itr_par_v_data->st_filesystem_replicate;
					con_s_peptide_data->v_p_replicate_data = itr_par_v_data->v_p_replicate_data;
					con_v_s_peptide_data.push_back(con_s_peptide_data);
				}
				con_str_peptide.clear();
				sw_peptide_filtered = 0;
			}
			itr_v_s_multinomial_element_data->v_s_peptide_data = con_v_s_peptide_data;
			con_v_s_peptide_data.clear();
			if ((IgFamily::DEBUG_MODE == 2) && (itr_v_s_multinomial_element_data->v_s_peptide_data.size() != 0)) {
				std::cout << "\n\n" << itr_v_s_multinomial_element_data->str_multinomial_element_name;
				std::cout << "   " << itr_v_s_multinomial_element_data->str_protein;
				std::vector<s_peptide_data*> con_itr_v_s_peptide_data = itr_v_s_multinomial_element_data->v_s_peptide_data;
				for (auto itr_v_s_peptide_data = con_itr_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_itr_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
					std::cout << "\n * " << (*itr_v_s_peptide_data)->str_peptide;
					std::cout << "   " << (*itr_v_s_peptide_data)->st_spectralcount;
				}
			}
		}
	}

	void create_v_s_multinomial_element_data_v_s_peptide_data_distinct_filtered(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data, std::vector<s_peptide_data> par_v_data_filtered_distinct) {
		string_type con_str_peptide_distinct_filtered = string_type();
		std::vector<s_peptide_data*> con_v_s_peptide_data_distinct_filtered;
		for (auto itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data) {
			for (auto itr_v_data_filtered_distinct = par_v_data_filtered_distinct.begin(); itr_v_data_filtered_distinct != par_v_data_filtered_distinct.end(); ++itr_v_data_filtered_distinct) {
				if (itr_v_s_multinomial_element_data->str_protein.find(itr_v_data_filtered_distinct->str_peptide) != std::string::npos) {
					s_peptide_data* con_s_peptide_data_distinct_filtered = new s_peptide_data;
					con_s_peptide_data_distinct_filtered->str_peptide = itr_v_data_filtered_distinct->str_peptide;
					con_s_peptide_data_distinct_filtered->st_spectralcount = itr_v_data_filtered_distinct->st_spectralcount;
					con_s_peptide_data_distinct_filtered->st_IgP = itr_v_data_filtered_distinct->st_IgP;
					con_s_peptide_data_distinct_filtered->v_s_denovo_peptide = itr_v_data_filtered_distinct->v_s_denovo_peptide;
					con_s_peptide_data_distinct_filtered->st_filesystem_replicate = itr_v_data_filtered_distinct->st_filesystem_replicate;
					con_s_peptide_data_distinct_filtered->v_p_replicate_data = itr_v_data_filtered_distinct->v_p_replicate_data;
					con_v_s_peptide_data_distinct_filtered.push_back(con_s_peptide_data_distinct_filtered);
				}
			}
			itr_v_s_multinomial_element_data->v_s_peptide_data_distinct_filtered = con_v_s_peptide_data_distinct_filtered;
			con_v_s_peptide_data_distinct_filtered.clear();
			if ((IgFamily::DEBUG_MODE == 2) && (itr_v_s_multinomial_element_data->v_s_peptide_data_distinct_filtered.size() != 0)) {
				std::cout << "\n\n" << itr_v_s_multinomial_element_data->str_multinomial_element_name;
				std::cout << "   " << itr_v_s_multinomial_element_data->str_protein;
				std::vector<s_peptide_data*> con_itr_v_s_peptide_data_distinct_filtered = itr_v_s_multinomial_element_data->v_s_peptide_data_distinct_filtered;
				for (auto itr_v_s_peptide_data_distinct_filtered = con_itr_v_s_peptide_data_distinct_filtered.begin(); itr_v_s_peptide_data_distinct_filtered != con_itr_v_s_peptide_data_distinct_filtered.end(); ++itr_v_s_peptide_data_distinct_filtered) {
					std::cout << "\n * " << (*itr_v_s_peptide_data_distinct_filtered)->str_peptide;
					std::cout << "   " << (*itr_v_s_peptide_data_distinct_filtered)->st_spectralcount;
				}
			}
		}
	}

	void create_v_s_multinomial_element_data_str_alignment(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data) {
		for (auto itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data) {
			string_type con_str_alignment = string_type();
			con_str_alignment.clear();
			for (auto i = size_type(); i < itr_v_s_multinomial_element_data->str_multinomial_element_name.length(); ++i) {
				con_str_alignment += '.';
			}
			std::vector<s_peptide_data*> con_itr_v_s_peptide_data_distinct_filtered = itr_v_s_multinomial_element_data->v_s_peptide_data_distinct_filtered;
			for (auto itr_v_s_peptide_data_distinct_filtered = con_itr_v_s_peptide_data_distinct_filtered.begin(); itr_v_s_peptide_data_distinct_filtered != con_itr_v_s_peptide_data_distinct_filtered.end(); ++itr_v_s_peptide_data_distinct_filtered) {
				char con_str_alignment_match;
				char con_str_alignment_test;
				size_type sw_alignment = size_type();
				size_type count_str_alignment = size_type();
				string_type con_itr_str_protein = itr_v_s_multinomial_element_data->str_multinomial_element_name;
				for (auto j = size_type(); j < con_itr_str_protein.length(); ++j) {
					con_str_alignment_match = con_itr_str_protein.at(j);
					con_str_alignment_test = (*itr_v_s_peptide_data_distinct_filtered)->str_peptide.at(0);
					if (sw_alignment == 1) {
						con_str_alignment_test = (*itr_v_s_peptide_data_distinct_filtered)->str_peptide.at(count_str_alignment);
						if (con_str_alignment_match == con_str_alignment_test) {
							++count_str_alignment;
						}
						else {
							count_str_alignment = 0;
							sw_alignment = 0;
						}
					}
					if ((con_str_alignment_match == con_str_alignment_test) && (count_str_alignment == 0)) {
						sw_alignment = 1;
						++count_str_alignment;
					}
					if (count_str_alignment == (*itr_v_s_peptide_data_distinct_filtered)->str_peptide.length()) {
						for (auto itr_modify_s_alignment = size_type(); itr_modify_s_alignment < (*itr_v_s_peptide_data_distinct_filtered)->str_peptide.length(); ++itr_modify_s_alignment) {
							con_str_alignment.at(j - itr_modify_s_alignment) = con_itr_str_protein.at(j - itr_modify_s_alignment);
						}
						count_str_alignment = 0;
						sw_alignment = 0;
					}
				}
			}
			itr_v_s_multinomial_element_data->str_alignment = con_str_alignment;
		}
	}

	void create_v_s_multinomial_element_data_st_totalspectralcount(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data) {
		for (auto itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data) {
			size_type con_st_totalspectralcount = size_type();
			std::vector<s_peptide_data*> con_v_s_peptide_data = itr_v_s_multinomial_element_data->v_s_peptide_data;
			for (auto itr_v_s_peptide_data = con_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
				con_st_totalspectralcount += (*itr_v_s_peptide_data)->st_spectralcount;
			}
			itr_v_s_multinomial_element_data->st_totalspectralcount = con_st_totalspectralcount;
		}
	}

	void create_v_s_multinomial_element_data_d_coverage(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data) {
		for (auto itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data) {
			double con_d_coverage;
			double count_str_alignment_true = double();
			for (auto itr_str_alignment = size_type(); itr_str_alignment < itr_v_s_multinomial_element_data->str_alignment.length(); ++itr_str_alignment) {
				if (itr_v_s_multinomial_element_data->str_alignment.at(itr_str_alignment) != '.') {
					++count_str_alignment_true;
				}
			}
			if ((itr_v_s_multinomial_element_data->str_alignment.length()) == 0) {
				std::cout << "error: division by zero";
				string_type str_catch_error;
				std::cin >> str_catch_error;
			}
			con_d_coverage = 100 * (count_str_alignment_true / (itr_v_s_multinomial_element_data->str_alignment.length()));
			itr_v_s_multinomial_element_data->d_coverage = con_d_coverage;
		}
	}

	void create_v_c_peptide_v_p_peptideassociation(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data, std::vector<s_peptide_data>& par_v_s_peptide_data) {
		for (auto itr_v_s_peptide_data = par_v_s_peptide_data.begin(); itr_v_s_peptide_data != par_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
			size_type sw_select_data = size_type();
			string_type con_st_select_data = string_type();
			for (auto j = size_type(); j < itr_v_s_peptide_data->str_peptide.length(); ++j) {
				if (itr_v_s_peptide_data->str_peptide.at(j) == '(' || ((itr_v_s_peptide_data->str_peptide.at(j) == '.') && (itr_v_s_peptide_data->str_peptide.length() > 2))) {
					sw_select_data = 1;
					if (itr_v_s_peptide_data->str_peptide.at(j + 1) == 's') {
						con_st_select_data.pop_back();
						con_st_select_data += itr_v_s_peptide_data->str_peptide.at(j + 5);
					}
				}
				if (sw_select_data == 0) {
					con_st_select_data += itr_v_s_peptide_data->str_peptide.at(j);
				}
				if (itr_v_s_peptide_data->str_peptide.at(j) == ')') {
					sw_select_data = 0;
				}
				if ((itr_v_s_peptide_data->str_peptide.at(j) == '.') && (itr_v_s_peptide_data->str_peptide.length() <= 2) && (sw_select_data == 0)) {
					con_st_select_data.clear();
				}
			}
			std::vector<string_type> con_v_str_peptideassociation;
			for (auto itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data) {
				if (itr_v_s_multinomial_element_data->str_multinomial_element_name.find(con_st_select_data) != std::string::npos) {
					con_v_str_peptideassociation.push_back(itr_v_s_multinomial_element_data->str_multinomial_element_name);
					itr_v_s_peptide_data->v_p_peptideassociation.push_back(std::make_pair(&(*itr_v_s_multinomial_element_data), double{ 0 }));
				}
			}
			itr_v_s_peptide_data->v_str_peptideassociation = con_v_str_peptideassociation;
			con_v_str_peptideassociation.clear();
			con_st_select_data.clear();
			sw_select_data = 0;
		}
	}

	void create_v_c_peptide_v_str_peptideassociation_distinct(std::vector<s_peptide_data>& par_v_s_peptide_data) {
		std::vector<string_type> con_v_str_peptideassociation_distinct;
		for (auto itr_v_s_peptide_data = par_v_s_peptide_data.begin(); itr_v_s_peptide_data != par_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
			std::vector<string_type> con_v_str_peptideassociation = itr_v_s_peptide_data->v_str_peptideassociation;
			for (auto itr_v_str_peptideassociation = con_v_str_peptideassociation.begin(); itr_v_str_peptideassociation != con_v_str_peptideassociation.end(); ++itr_v_str_peptideassociation) {
				size_type sw_peptideassociation_distinct = size_type();
				string_type con_str_peptideassociation_distinct = string_type();
				for (auto i = size_type(); i < (*itr_v_str_peptideassociation).length(); ++i) {
					if ((*itr_v_str_peptideassociation).at(i) == '*') {
						sw_peptideassociation_distinct = 1;
					}
					if (sw_peptideassociation_distinct == 0) {
						con_str_peptideassociation_distinct += (*itr_v_str_peptideassociation).at(i);
					}
				}
				if (std::find(con_v_str_peptideassociation_distinct.begin(), con_v_str_peptideassociation_distinct.end(), con_str_peptideassociation_distinct) == con_v_str_peptideassociation_distinct.end()) {
					con_v_str_peptideassociation_distinct.push_back(con_str_peptideassociation_distinct);
				}
				sw_peptideassociation_distinct = 0;
				con_str_peptideassociation_distinct.clear();
			}
			itr_v_s_peptide_data->v_str_peptideassociation_distinct = con_v_str_peptideassociation_distinct;
			con_v_str_peptideassociation_distinct.clear();
		}
	}

	void create_v_c_peptide_v_p_peptideassociation_distinct(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data_distinct, std::vector<s_peptide_data>& par_v_s_peptide_data) {
		std::vector<s_multinomial_element_data*> con_v_p_peptideassociation_distinct;
		std::vector<string_type> con_v_str_peptideassociation_distinct;
		for (auto itr_v_s_peptide_data = par_v_s_peptide_data.begin(); itr_v_s_peptide_data != par_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
			std::vector<std::pair<s_multinomial_element_data*, double>> con_v_p_peptideassociation = itr_v_s_peptide_data->v_p_peptideassociation;
			for (auto itr_v_p_peptideassociation = con_v_p_peptideassociation.begin(); itr_v_p_peptideassociation != con_v_p_peptideassociation.end(); ++itr_v_p_peptideassociation) {
				size_type sw_peptideassociation_distinct = size_type();
				string_type con_str_peptideassociation_distinct = string_type();
				for (auto i = size_type(); i < std::get<0>(*itr_v_p_peptideassociation)->str_multinomial_element_name.length(); ++i) {
					if (std::get<0>(*itr_v_p_peptideassociation)->str_multinomial_element_name.at(i) == '*') {
						sw_peptideassociation_distinct = 1;
					}
					if (sw_peptideassociation_distinct == 0) {
						con_str_peptideassociation_distinct += std::get<0>(*itr_v_p_peptideassociation)->str_multinomial_element_name.at(i);
					}
				}
				if (std::find(con_v_str_peptideassociation_distinct.begin(), con_v_str_peptideassociation_distinct.end(), con_str_peptideassociation_distinct) == con_v_str_peptideassociation_distinct.end()) {
					con_v_str_peptideassociation_distinct.push_back(con_str_peptideassociation_distinct);
					for (auto itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data_distinct.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data_distinct.end(); ++itr_v_s_multinomial_element_data) {
						if (itr_v_s_multinomial_element_data->str_multinomial_element_name == con_str_peptideassociation_distinct) {
							itr_v_s_peptide_data->v_p_peptideassociation_distinct.push_back(std::make_pair(&(*itr_v_s_multinomial_element_data), double{ 1 }));
						}
					}
				}
				sw_peptideassociation_distinct = 0;
				con_str_peptideassociation_distinct.clear();
			}
			itr_v_s_peptide_data->v_str_peptideassociation_distinct = con_v_str_peptideassociation_distinct;
			con_v_str_peptideassociation_distinct.clear();
		}
	}

	inline void create_global_score_mean(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data) {
		double con_d_score = double();
		size_type con_st_nonzero_score = size_type();
		for (auto itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data) {
			if (itr_v_s_multinomial_element_data->d_score >= double{ 1.5 }) {
				con_d_score += itr_v_s_multinomial_element_data->d_score;
				++con_st_nonzero_score;
			}
		}
		if (con_st_nonzero_score == 0) {
			con_st_nonzero_score = 1;
			con_d_score = 1;
		}
		con_d_score = con_d_score / con_st_nonzero_score;
		IgFamily::SCORE_MEAN = con_d_score;
	}

	inline double log_basechange(double d, double base) {
		return (log(d) / log(base));
	}

	void create_v_s_multinomial_element_data_d_score(std::vector<s_peptide_data>& par_v_s_peptide_data, std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data, std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data_distinct) {
		for (std::vector<s_multinomial_element_data>::iterator itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data) {
			string_type con_str_multinomial_element_name = itr_v_s_multinomial_element_data->str_multinomial_element_name;
			size_type sw_peptideassociation_distinct = size_type();
			string_type con_str_multinomial_element_name_distinct = string_type();
			for (auto i = size_type(); i < con_str_multinomial_element_name.length(); ++i) {
				if (con_str_multinomial_element_name.at(i) == '*') {
					sw_peptideassociation_distinct = 1;
				}
				if (sw_peptideassociation_distinct == 0) {
					con_str_multinomial_element_name_distinct += con_str_multinomial_element_name.at(i);
				}
			}	
			std::vector<s_peptide_data*> con_v_s_peptide_data = itr_v_s_multinomial_element_data->v_s_peptide_data;
			double con_d_score = double();
			double con_d_parameter_gene_family = double();
			double con_d_parameter_gene_family_train = double();
			double con_d_parameter_gene_family_weight = double();
			double con_d_parameter_IgP = double();
			for (auto itr_v_s_peptide_data = con_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
				con_d_parameter_gene_family_train = pow(((itr_v_s_multinomial_element_data->d_score + (double{ 100 } *IgFamily::SCORE_MEAN)) / (double{ 100 } *IgFamily::SCORE_MEAN)), double{ 1.6 });
				con_d_parameter_gene_family_train = 1;
				for (auto itr_v_p_peptideassociation_distinct = (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
					if (std::get<0>(*itr_v_p_peptideassociation_distinct)->str_multinomial_element_name == con_str_multinomial_element_name_distinct) {
						con_d_parameter_gene_family_weight = std::get<1>(*itr_v_p_peptideassociation_distinct);
					}
				}
				con_d_parameter_gene_family = pow((con_d_parameter_gene_family_weight / (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.size()), double{ double{ 1 } / con_d_parameter_gene_family_train });
				con_d_parameter_IgP = (log_basechange(((double((*itr_v_s_peptide_data)->st_IgP) + double(IgFamily::PARSE_THRESHOLD_IgP)) / (double{ 2 } *double(IgFamily::PARSE_THRESHOLD_IgP))), double{ 9 }) + 1);
				con_d_score += (con_d_parameter_IgP * (*itr_v_s_peptide_data)->st_spectralcount * con_d_parameter_gene_family);
			}
			itr_v_s_multinomial_element_data->d_score = con_d_score;
			con_str_multinomial_element_name_distinct.clear();
		}
		update_v_s_multinomial_element_distinctpolymorphism_data(par_v_s_multinomial_element_data_distinct);

		for (auto itr_v_s_peptide_data = par_v_s_peptide_data.begin(); itr_v_s_peptide_data != par_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {		
			double d_v_p_peptideassociation_d_score_sum = double();
			double d_v_p_peptideassociation_d_score_mean = double();
			for (auto itr_v_p_peptideassociation_distinct = itr_v_s_peptide_data->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != itr_v_s_peptide_data->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
				d_v_p_peptideassociation_d_score_sum += std::get<0>(*itr_v_p_peptideassociation_distinct)->d_score;
			}
			d_v_p_peptideassociation_d_score_mean = (d_v_p_peptideassociation_d_score_sum / itr_v_s_peptide_data->v_p_peptideassociation_distinct.size());
			for (auto itr_v_p_peptideassociation_distinct_2 = itr_v_s_peptide_data->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct_2 != itr_v_s_peptide_data->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct_2) {
				std::get<1>(*itr_v_p_peptideassociation_distinct_2) = (((std::get<0>(*itr_v_p_peptideassociation_distinct_2)->d_score * double(1)) + (d_v_p_peptideassociation_d_score_mean)) / (double{ 2 } *d_v_p_peptideassociation_d_score_mean));
			}
		}
	}

	inline bool predicate_v_s_multinomial_element_data_d_score(const s_multinomial_element_data& i, const s_multinomial_element_data& j) {
		return (i.d_score > j.d_score);
	}

	inline void sort_v_s_multinomial_element_data_d_score(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data) {
		std::sort(par_v_s_multinomial_element_data.begin(), par_v_s_multinomial_element_data.end(), predicate_v_s_multinomial_element_data_d_score);
	}

	inline bool predicate_v_s_peptide_data_str_peptide(const s_peptide_data& i, const s_peptide_data& j) {
		return (i.str_peptide < j.str_peptide);
	}

	inline void sort_v_s_peptide_data_str_peptide(std::vector<s_peptide_data>& par_s_peptide_data) {
		std::sort(par_s_peptide_data.begin(), par_s_peptide_data.end(), predicate_v_s_peptide_data_str_peptide);
	}

	inline bool predicate_v_s_peptide_data_st_spectralcount(const s_peptide_data& i, const s_peptide_data& j) {
		return (i.st_spectralcount > j.st_spectralcount);
	}

	inline void sort_v_s_peptide_data_st_spectralcount(std::vector<s_peptide_data>& par_s_peptide_data) {
		std::sort(par_s_peptide_data.begin(), par_s_peptide_data.end(), predicate_v_s_peptide_data_st_spectralcount);
	}

	inline bool test2(const s_peptide_data* i, const s_peptide_data* j) {
		return (i->str_peptide < j->str_peptide);
	}

	inline void sort_c_analysis_v_s_peptide_data_s_peptide(std::vector<s_peptide_data*>& par_v_s_peptide_data) {
		std::sort(par_v_s_peptide_data.begin(), par_v_s_peptide_data.end(), test2);
	}

	std::vector<s_multinomial_element_data*> map_v_s_multinomial_element_data_by_score(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data) {
		std::vector<s_multinomial_element_data*> con_v_s_multinomial_element_data;
		s_multinomial_element_data* con_c_analysis;
		size_type st_count_map = size_type();
		while (st_count_map < par_v_s_multinomial_element_data.size()) {
			double d_compare_map = double();
			for (std::vector<s_multinomial_element_data>::iterator itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data) {
				bool b_compare_map = ((std::find(con_v_s_multinomial_element_data.begin(), con_v_s_multinomial_element_data.end(), &(*itr_v_s_multinomial_element_data))) == con_v_s_multinomial_element_data.end());
				if ((d_compare_map <= itr_v_s_multinomial_element_data->d_score) && b_compare_map) {
					d_compare_map = itr_v_s_multinomial_element_data->d_score;
					con_c_analysis = &(*itr_v_s_multinomial_element_data);
				}
			}
			con_v_s_multinomial_element_data.push_back(con_c_analysis);
			++st_count_map;
		}
		return con_v_s_multinomial_element_data;
	}

	void train_v_s_multinomial_element_data_d_score(std::vector<s_peptide_data>& par_v_s_peptide_data, std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data, std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data_distinct) {
		double delta_score_mean = double();
		for (IgFamily::GLOBAL_ITERATOR; IgFamily::GLOBAL_ITERATOR < IgFamily::ITERATE_TRAIN_SCORE; ++IgFamily::GLOBAL_ITERATOR) {
			delta_score_mean = IgFamily::SCORE_MEAN;
			create_global_score_mean(par_v_s_multinomial_element_data_distinct);
			delta_score_mean -= IgFamily::SCORE_MEAN;
			std::cout << "\n\n iteration - " << IgFamily::GLOBAL_ITERATOR;
			std::cout << "   gene family score mean - " << std::fixed << std::setprecision(2) << IgFamily::SCORE_MEAN;
			if (IgFamily::GLOBAL_ITERATOR != 0) {
				std::cout << "   mean sample difference - " << std::fixed << std::setprecision(2) << abs(delta_score_mean);
			}
			fpf_data::create_v_s_multinomial_element_data_d_score(par_v_s_peptide_data, par_v_s_multinomial_element_data, par_v_s_multinomial_element_data_distinct);
			fpf_data::create_v_s_multinomial_element_data_v_s_peptide_data(par_v_s_multinomial_element_data, par_v_s_peptide_data);
			//if (abs(delta_score_mean) < (IgFamily::SCORE_MEAN / double{ 10000 })) {
			//	break;
			//}
		}
		delta_score_mean = IgFamily::SCORE_MEAN;
		create_global_score_mean(par_v_s_multinomial_element_data_distinct);
		delta_score_mean -= IgFamily::SCORE_MEAN;
		std::cout << "\n\n iteration - " << IgFamily::GLOBAL_ITERATOR;
		std::cout << "   gene family score mean - " << std::fixed << std::setprecision(2) << IgFamily::SCORE_MEAN;
		std::cout << "   mean sample difference - " << std::fixed << std::setprecision(2) << abs(delta_score_mean);
		fpf_data::create_v_s_multinomial_element_data_v_s_peptide_data(par_v_s_multinomial_element_data, par_v_s_peptide_data);
	}

	std::vector<s_multinomial_element_data*> map_v_s_multinomial_element_data_by_genefamily(std::vector<s_multinomial_element_data>& par_v_s_multinomial_element_data) {
		std::vector<s_multinomial_element_data*> con_v_s_multinomial_element_data;
		s_multinomial_element_data* con_c_analysis;
		size_type st_count_map = size_type();
		while (st_count_map < par_v_s_multinomial_element_data.size()) {
			string_type str_compare_map = string_type();
			for (auto itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data) {
				bool b_compare_map = ((std::find(con_v_s_multinomial_element_data.begin(), con_v_s_multinomial_element_data.end(), &(*itr_v_s_multinomial_element_data))) == con_v_s_multinomial_element_data.end());
				if ((str_compare_map < itr_v_s_multinomial_element_data->str_multinomial_element_name) && b_compare_map) {
					str_compare_map = itr_v_s_multinomial_element_data->str_multinomial_element_name;
				}
			}
			con_v_s_multinomial_element_data.push_back(con_c_analysis);
			str_compare_map.clear();
			++st_count_map;
		}
		return con_v_s_multinomial_element_data;
	}

	void output(string_type par_str_fin, std::vector<s_multinomial_element_data*> par_v_s_multinomial_element_data_distinct, std::vector<s_multinomial_element_data*> par_v_s_multinomial_element_data) {
		
		std::cout << "\n\n ...streaming output";
		std::string output_v_s_multinomial_element_data = par_str_fin + "_output.txt";
		std::ofstream fout_v_s_multinomial_element_data;
		fout_v_s_multinomial_element_data.open(output_v_s_multinomial_element_data);
		fout_v_s_multinomial_element_data << "-- IgFamily " << IgFamily::version << " --\n\n\n";
		fout_v_s_multinomial_element_data << "Input file : " << IgFamily::INPUT_CSV;
		for (auto itr_v_s_multinomial_element_data_distinct = par_v_s_multinomial_element_data_distinct.begin(); itr_v_s_multinomial_element_data_distinct != par_v_s_multinomial_element_data_distinct.end(); ++itr_v_s_multinomial_element_data_distinct) {
			if (itr_v_s_multinomial_element_data_distinct == par_v_s_multinomial_element_data_distinct.begin()) {
				fout_v_s_multinomial_element_data << "\n\n\n\nTop gene families - ";
				fout_v_s_multinomial_element_data << "     Score mean - " << std::fixed << std::setprecision(2) << IgFamily::SCORE_MEAN;
				fout_v_s_multinomial_element_data << "\n\n";
			}
			if ((itr_v_s_multinomial_element_data_distinct == par_v_s_multinomial_element_data_distinct.end()) || ((*itr_v_s_multinomial_element_data_distinct)->d_score < double{ 0.01 })) {
				break;
			}
			fout_v_s_multinomial_element_data << "\n  *  " << (*itr_v_s_multinomial_element_data_distinct)->str_multinomial_element_name;
			fout_v_s_multinomial_element_data << "   " << (*itr_v_s_multinomial_element_data_distinct)->str_species;
			for (size_type itr_str_multinomial_element_name_output = (*itr_v_s_multinomial_element_data_distinct)->str_multinomial_element_name.length() + (*itr_v_s_multinomial_element_data_distinct)->str_species.length(); itr_str_multinomial_element_name_output < 35; ++itr_str_multinomial_element_name_output) {
				fout_v_s_multinomial_element_data << " ";
			}
			fout_v_s_multinomial_element_data << "   Score - " << std::fixed << std::setprecision(2) << (*itr_v_s_multinomial_element_data_distinct)->d_score;
			std::vector<fpf_data::s_peptide_data*> con_itr_v_s_peptide_data = (*itr_v_s_multinomial_element_data_distinct)->v_s_peptide_data;
			double d_peptide_1genefamily = double();
			double d_peptide_2genefamily = double();
			double d_peptide_3genefamily = double();
			for (auto itr_v_s_peptide_data = con_itr_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_itr_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
				string_type con_str_multinomial_element_name = (*itr_v_s_multinomial_element_data_distinct)->str_multinomial_element_name;
				size_type sw_peptideassociation_distinct = size_type();
				string_type con_str_multinomial_element_name_distinct = string_type();
				for (auto i = size_type(); i < con_str_multinomial_element_name.length(); ++i) {
					if (con_str_multinomial_element_name.at(i) == '*') {
						sw_peptideassociation_distinct = 1;
					}
					if (sw_peptideassociation_distinct == 0) {
						con_str_multinomial_element_name_distinct += con_str_multinomial_element_name.at(i);
					}
				}
				double con_d_score = double();
				double con_d_parameter_gene_family = double();
				double con_d_parameter_gene_family_train = double();
				double con_d_parameter_gene_family_weight = double();
				double con_d_parameter_IgP = double();
				con_d_parameter_gene_family_train = pow((((*itr_v_s_multinomial_element_data_distinct)->d_score + (double{ 100 } *IgFamily::SCORE_MEAN)) / (double{ 100 } *IgFamily::SCORE_MEAN)), double{ 1.6 });
				con_d_parameter_gene_family_train = 1;
				for (auto itr_v_p_peptideassociation_distinct = (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
					if (std::get<0>(*itr_v_p_peptideassociation_distinct)->str_multinomial_element_name == con_str_multinomial_element_name_distinct) {
						con_d_parameter_gene_family_weight = std::get<1>(*itr_v_p_peptideassociation_distinct);
					}
				}
				con_d_parameter_gene_family = pow((con_d_parameter_gene_family_weight / (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.size()), double{ double{ 1 } / con_d_parameter_gene_family_train });
				con_d_parameter_IgP = (log_basechange(((double((*itr_v_s_peptide_data)->st_IgP) + double(IgFamily::PARSE_THRESHOLD_IgP)) / (double{ 2 } *double(IgFamily::PARSE_THRESHOLD_IgP))), double{ 9 }) + 1);
				con_d_score += (con_d_parameter_IgP * (*itr_v_s_peptide_data)->st_spectralcount * con_d_parameter_gene_family);
				if ((*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.size() == 1) {
					d_peptide_1genefamily += con_d_score;
				}
				if ((*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.size() == 2) {
					d_peptide_2genefamily += con_d_score;
				}
				if ((*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.size() == 3) {
					d_peptide_3genefamily += con_d_score;
				}
				con_str_multinomial_element_name_distinct.clear();
			}
			d_peptide_1genefamily /= (*itr_v_s_multinomial_element_data_distinct)->d_score;
			d_peptide_2genefamily /= (*itr_v_s_multinomial_element_data_distinct)->d_score;
			d_peptide_3genefamily /= (*itr_v_s_multinomial_element_data_distinct)->d_score;
			for (size_type itr_v_s_peptide_data_output = ((log10((*itr_v_s_multinomial_element_data_distinct)->d_score) >= double{ 1 }) ? (size_type(log10((*itr_v_s_multinomial_element_data_distinct)->d_score)) + size_type(1)) : size_type(1)); itr_v_s_peptide_data_output < 6; ++itr_v_s_peptide_data_output) {
				fout_v_s_multinomial_element_data << " ";
			}
			fout_v_s_multinomial_element_data << "( 1GF - " << std::fixed << std::setprecision(2) << d_peptide_1genefamily;
			fout_v_s_multinomial_element_data << ", 2GF - " << std::fixed << std::setprecision(2) << d_peptide_2genefamily;
			fout_v_s_multinomial_element_data << ", 3GF - " << std::fixed << std::setprecision(2) << d_peptide_3genefamily;
			fout_v_s_multinomial_element_data << " )";
			d_peptide_1genefamily = 0;
			d_peptide_2genefamily = 0;
			d_peptide_3genefamily = 0;
			//for (std::vector<fpf_data::s_multinomial_element_data*>::iterator itr_v_s_multinomial_element_data_distinct_v_c_polyassociation = (itr_v_s_multinomial_element_data_distinct + i)->ref_v_s_multinomial_element_polyassociation().begin(); itr_v_s_multinomial_element_data_distinct_v_c_polyassociation != (itr_v_s_multinomial_element_data_distinct + i)->ref_v_s_multinomial_element_polyassociation().end(); ++itr_v_s_multinomial_element_data_distinct_v_c_polyassociation) {
			//	fout_v_s_multinomial_element_data << "\n - - * " << (*itr_v_s_multinomial_element_data_distinct_v_c_polyassociation)->str_multinomial_element_name;
			//	for (size_type itr_str_multinomial_element_name_output = (itr_v_s_multinomial_element_data_distinct + i)->str_multinomial_element_name.length(); itr_str_multinomial_element_name_output < 20; ++itr_str_multinomial_element_name_output) {
			//		fout_v_s_multinomial_element_data << " ";
			//	}
			//	fout_v_s_multinomial_element_data << " - " << (*itr_v_s_multinomial_element_data_distinct_v_c_polyassociation)->d_score;
			//}
		}

		//fout_v_s_multinomial_element_data << "\n\n";
		//for (std::vector<fpf_data::s_multinomial_element_data*>::iterator itr_v_s_multinomial_element_data_distinct = par_v_s_multinomial_element_data_distinct.begin(); itr_v_s_multinomial_element_data_distinct != par_v_s_multinomial_element_data_distinct.end(); ++itr_v_s_multinomial_element_data_distinct) {
		//	if ((*itr_v_s_multinomial_element_data_distinct)->d_score >= 5) {
		//		fout_v_s_multinomial_element_data << "\n" << (*itr_v_s_multinomial_element_data_distinct)->str_multinomial_element_name;
		//	}
		//}
		//fout_v_s_multinomial_element_data << "\n\n";
		//for (std::vector<fpf_data::s_multinomial_element_data*>::iterator itr_v_s_multinomial_element_data_distinct = par_v_s_multinomial_element_data_distinct.begin(); itr_v_s_multinomial_element_data_distinct != par_v_s_multinomial_element_data_distinct.end(); ++itr_v_s_multinomial_element_data_distinct) {
		//	if ((*itr_v_s_multinomial_element_data_distinct)->d_score >= 5) {
		//		fout_v_s_multinomial_element_data << "\n" << (*itr_v_s_multinomial_element_data_distinct)->d_score;
		//	}
		//}

		for (auto itr_v_s_multinomial_element_data = par_v_s_multinomial_element_data.begin(); itr_v_s_multinomial_element_data != par_v_s_multinomial_element_data.end(); ++itr_v_s_multinomial_element_data) {
			sort_c_analysis_v_s_peptide_data_s_peptide((*itr_v_s_multinomial_element_data)->v_s_peptide_data);
			IgFamily::SCORE_THRESHOLD = (IgFamily::SCORE_MEAN / double{ 2 });
			if (((*itr_v_s_multinomial_element_data)->v_s_peptide_data.size() != 0) && ((*itr_v_s_multinomial_element_data)->d_score >= 3)) {
				fout_v_s_multinomial_element_data << "\n\n\n\n\n" << (*itr_v_s_multinomial_element_data)->str_multinomial_element_name;
				fout_v_s_multinomial_element_data << "   " << (*itr_v_s_multinomial_element_data)->str_species;
				fout_v_s_multinomial_element_data << "   (top: ";
				double test2 = double();
				string_type test3 = string_type();
				for (auto test = (*itr_v_s_multinomial_element_data)->v_s_multinomial_element_polyassociation.begin(); test != (*itr_v_s_multinomial_element_data)->v_s_multinomial_element_polyassociation.end(); ++test) {
					if ((*test)->d_score > test2) {
						test2 = (*test)->d_score;
						test3 = (*test)->str_multinomial_element_name;
					}
				}
				fout_v_s_multinomial_element_data << test3 << ")";
				fout_v_s_multinomial_element_data << "   #SC - " << (*itr_v_s_multinomial_element_data)->st_totalspectralcount;
				fout_v_s_multinomial_element_data << "   % - " << std::fixed << std::setprecision(2) << (*itr_v_s_multinomial_element_data)->d_coverage;
				fout_v_s_multinomial_element_data << "     Score - " << std::fixed << std::setprecision(2) << (*itr_v_s_multinomial_element_data)->d_score;
				fout_v_s_multinomial_element_data << "\n\n" << (*itr_v_s_multinomial_element_data)->str_multinomial_element_name;
				fout_v_s_multinomial_element_data << "\n" << (*itr_v_s_multinomial_element_data)->str_alignment;
				fout_v_s_multinomial_element_data << "\n";
				std::vector<fpf_data::s_peptide_data*> con_itr_v_s_peptide_data = (*itr_v_s_multinomial_element_data)->v_s_peptide_data;
				for (auto itr_v_s_peptide_data = con_itr_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_itr_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
					//if ((*itr_v_s_peptide_data)->return_v_str_peptideassociation_distinct().size() <= 5) {
					fout_v_s_multinomial_element_data << "\n * " << (*itr_v_s_peptide_data)->str_peptide;
					for (size_type itr_v_s_peptide_data_output = (*itr_v_s_peptide_data)->str_peptide.length(); itr_v_s_peptide_data_output < 50; ++itr_v_s_peptide_data_output) {
						fout_v_s_multinomial_element_data << " ";
					}

					string_type con_s_genefamily = (*itr_v_s_multinomial_element_data)->str_multinomial_element_name;
					size_type sw_s_peptideassociation_distinct = size_type();
					string_type con_s_genefamily_distinct = string_type();
					for (unsigned itr_s_genefamily = 0; itr_s_genefamily < con_s_genefamily.length(); ++itr_s_genefamily) {
						if (con_s_genefamily.at(itr_s_genefamily) == '*') {
							sw_s_peptideassociation_distinct = 1;
						}
						if (sw_s_peptideassociation_distinct == 0) {
							con_s_genefamily_distinct += con_s_genefamily.at(itr_s_genefamily);
						}
					}
					double con_d_score = double();
					double con_d_parameter_gene_family = double();
					double con_d_parameter_gene_family_train = double();
					double con_d_parameter_gene_family_weight = double();
					double con_d_parameter_IgP = double();
					//con_d_parameter_gene_family_train = pow((((*itr_v_s_multinomial_element_data)->d_score + (double{ 100 } *IgFamily::SCORE_MEAN)) / (double{ 100 } *IgFamily::SCORE_MEAN)), double{ 1.6 });
					con_d_parameter_gene_family_train = 1;
					for (auto itr_v_p_peptideassociation_distinct = (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
						if (std::get<0>(*itr_v_p_peptideassociation_distinct)->str_multinomial_element_name == con_s_genefamily_distinct) {
							con_d_parameter_gene_family_weight = std::get<1>(*itr_v_p_peptideassociation_distinct);
						}
					}
					con_d_parameter_gene_family = pow((con_d_parameter_gene_family_weight / (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.size()), double{ double{ 1 } / con_d_parameter_gene_family_train });
					con_d_parameter_IgP = (log_basechange(((double((*itr_v_s_peptide_data)->st_IgP) + double(IgFamily::PARSE_THRESHOLD_IgP)) / (double{ 2 } *double(IgFamily::PARSE_THRESHOLD_IgP))), double{ 9 }) + 1);
					con_d_score += (con_d_parameter_IgP * (*itr_v_s_peptide_data)->st_spectralcount * con_d_parameter_gene_family);
					con_s_genefamily_distinct.clear();

					fout_v_s_multinomial_element_data << "   Score - " << std::fixed << std::setprecision(2) << con_d_score;
					for (size_type itr_v_s_peptide_data_output = ((log10(con_d_score) >= double{ 1 }) ? (size_type(log10(con_d_score)) + size_type(1)) : size_type(1)); itr_v_s_peptide_data_output < 6; ++itr_v_s_peptide_data_output) {
						fout_v_s_multinomial_element_data << " ";
					}
					//fout_v_s_multinomial_element_data << "   n - " << std::fixed << std::setprecision(2) << con_d_score;
					//for (size_type itr_v_s_peptide_data_output = ((log10(con_d_score) >= double{ 1 }) ? (size_type(log10(con_d_score)) + size_type(1)) : size_type(1)); itr_v_s_peptide_data_output < 6; ++itr_v_s_peptide_data_output) {
					//	fout_v_s_multinomial_element_data << " ";
					//}
					fout_v_s_multinomial_element_data << "   SC - " << (*itr_v_s_peptide_data)->st_spectralcount;
					for (size_type itr_v_s_peptide_data_output = size_type(log10(double((*itr_v_s_peptide_data)->st_spectralcount))); itr_v_s_peptide_data_output < 2; ++itr_v_s_peptide_data_output) {
						fout_v_s_multinomial_element_data << " ";
					}
					//fout_v_s_multinomial_element_data << "   GF par - " << std::fixed << std::setprecision(2) << con_d_parameter_gene_family;
					//for (size_type itr_v_s_peptide_data_output = ((log10(con_d_score) >= double{ 1 }) ? size_type(log10(con_d_score)) : size_type(1)); itr_v_s_peptide_data_output < 4; ++itr_v_s_peptide_data_output) {
					//	fout_v_s_multinomial_element_data << " ";
					//}
					//fout_v_s_multinomial_element_data << "   Hom score - " << (*itr_v_s_peptide_data)->return_st_IgP() << "  ";
					//for (size_type itr_v_s_peptide_data_output = size_type(log10(double((*itr_v_s_peptide_data)->return_st_IgP()))); itr_v_s_peptide_data_output < 2; ++itr_v_s_peptide_data_output) {
					//	fout_v_s_multinomial_element_data << " ";
					//}
					//fout_v_s_multinomial_element_data << "   Train score - " << con_d_parameter_gene_family_train << "  ";
					//for (size_type itr_v_s_peptide_data_output = size_type(log10(con_d_parameter_gene_family_train)); itr_v_s_peptide_data_output < 3; ++itr_v_s_peptide_data_output) {
					//	fout_v_s_multinomial_element_data << " ";
					//}
					//fout_v_s_multinomial_element_data << "   Hom par - " << std::fixed << std::setprecision(2) << con_d_parameter_IgP;
					//for (size_type itr_v_s_peptide_data_output = size_type(log10(con_d_parameter_IgP)); itr_v_s_peptide_data_output < 3; ++itr_v_s_peptide_data_output) {
					//	fout_v_s_multinomial_element_data << " ";
					//}
					fout_v_s_multinomial_element_data << "   " << "Member of " << (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.size() << " gene families - ";
					std::vector<std::pair<fpf_data::s_multinomial_element_data*, double>> con_itr_v_p_peptideassociation_distinct = (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct;
					for (std::vector<std::pair<fpf_data::s_multinomial_element_data*, double>>::iterator itr_v_p_peptideassociation_distinct = con_itr_v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != con_itr_v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
						fout_v_s_multinomial_element_data << std::get<0>(*itr_v_p_peptideassociation_distinct)->str_multinomial_element_name;
						fout_v_s_multinomial_element_data << "(" << std::get<1>(*itr_v_p_peptideassociation_distinct) << ")";
						if ((itr_v_p_peptideassociation_distinct + 1) != con_itr_v_p_peptideassociation_distinct.end()) {
							fout_v_s_multinomial_element_data << ", ";
						}
					}
				}
			}
		}
	}
}

#endif