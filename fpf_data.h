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



namespace fpf_data {

	class c_genefamily_data;
	struct s_peptide_data;
	struct s_denovo_peptide;

	typedef std::string string_type;
	typedef size_t size_type;

	class c_genefamily_data {
	public:
		c_genefamily_data() {
			d_score = { 1 };
		};

		~c_genefamily_data() {
		};

		inline void set_str_genefamily(string_type par_str_genefamily) {
			str_genefamily = par_str_genefamily;
		};

		inline void set_str_genefamily_class(string_type par_str_genefamily_class) {
			str_genefamily_class = par_str_genefamily_class;
		};

		inline void set_str_protein(string_type par_str_protein) {
			str_protein = par_str_protein;
		};

		inline void set_v_s_peptide_data(std::vector<s_peptide_data*> par_v_s_peptide_data) {
			v_s_peptide_data = par_v_s_peptide_data;
		};

		inline void set_v_s_peptide_data_distinct_filtered(std::vector<s_peptide_data*> par_v_s_peptide_data_distinct_filtered) {
			v_s_peptide_data_distinct_filtered = par_v_s_peptide_data_distinct_filtered;
		};

		inline void set_str_alignment(string_type par_str_alignment) {
			str_alignment = par_str_alignment;
		};

		inline void set_st_totalspectralcount(size_type par_st_totalspectralcount) {
			st_totalspectralcount = par_st_totalspectralcount;
		};

		inline void set_d_coverage(double par_d_coverage) {
			d_coverage = par_d_coverage;
		};

		inline void set_d_score(double par_d_score) {
			d_score = par_d_score;
		};

		inline std::vector<c_genefamily_data*>& ref_v_c_analysis_polyassociation() {
			return v_c_analysis_polyassociation;
		};

		inline const string_type return_str_genefamily() const {
			return str_genefamily;
		};

		inline const string_type return_str_genefamily_class() const {
			return str_genefamily_class;
		};

		inline const string_type return_str_protein() const {
			return str_protein;
		};

		inline std::vector<s_peptide_data*>& return_v_s_peptide_data() {
			return v_s_peptide_data;
		};

		inline std::vector<s_peptide_data*>& return_v_s_peptide_data_distinct_filtered() {
			return v_s_peptide_data_distinct_filtered;
		};

		inline const string_type return_str_alignment() const {
			return str_alignment;
		};

		inline const size_type return_st_totalspectralcount() const {
			return st_totalspectralcount;
		};

		inline const double return_d_coverage() const {
			return d_coverage;
		};

		inline const double return_d_score() const {
			return d_score;
		};

	public:
		string_type str_species;
		string_type str_genefamily;
		string_type str_protein;

	private:
		std::vector<c_genefamily_data*> v_c_analysis_polyassociation;		
		string_type str_genefamily_class;	
		std::vector<s_peptide_data*> v_s_peptide_data;
		std::vector<s_peptide_data*> v_s_peptide_data_distinct_filtered;
		string_type str_alignment;
		size_type st_totalspectralcount;
		double d_coverage;
		double d_score;
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
		bool b_replicate_merged = bool();
		size_type st_filesystem_replicate = size_type{ 1 };
		std::vector<std::tuple<string_type, size_type, size_type>> v_p_replicate_data;
		std::vector<string_type> v_str_peptideassociation;
		std::vector<string_type> v_str_peptideassociation_distinct;
		std::vector<std::pair<c_genefamily_data*, double>> v_p_peptideassociation;
		std::vector<std::pair<c_genefamily_data*, double>> v_p_peptideassociation_distinct;		
	};

	struct s_denovo_peptide {
		char ch_aminoacid;
		double d_denovo_localconfidence;
	};

	std::vector<c_genefamily_data> create_v_c_analysis(std::vector<fpf_parse::c_parse_FASTA> par_v_c_parse_FASTA) {

		// NOTE: Vector reallocation will invalidate all pointers to c_genefamily_data references.
		//		 Do not reallocate after return!

		std::vector<c_genefamily_data> con_v_c_analysis;
		c_genefamily_data con_c_analysis;
		size_type con_st_ID = size_type();
		for (auto itr_par_v_c_parse_FASTA : par_v_c_parse_FASTA) {
			++con_st_ID;
			auto find_v_c_analysis = std::find_if(con_v_c_analysis.begin(), con_v_c_analysis.end(),
				[itr_par_v_c_parse_FASTA](c_genefamily_data par_c_genefamily_data) {
				return par_c_genefamily_data.return_str_genefamily() == itr_par_v_c_parse_FASTA.return_str_parse_FASTA_genefamily(); });
			if (find_v_c_analysis == con_v_c_analysis.end()) {
				con_c_analysis.set_str_genefamily(itr_par_v_c_parse_FASTA.return_str_parse_FASTA_genefamily());
				con_c_analysis.set_str_genefamily_class(itr_par_v_c_parse_FASTA.return_str_parse_FASTA_genefamily_class());
				con_c_analysis.str_species = itr_par_v_c_parse_FASTA.return_str_parse_FASTA_species();
				con_c_analysis.set_str_protein(itr_par_v_c_parse_FASTA.return_str_parse_FASTA_protein());
				con_v_c_analysis.push_back(con_c_analysis);
			}
			else {
				find_v_c_analysis->set_str_protein(find_v_c_analysis->return_str_protein() + itr_par_v_c_parse_FASTA.return_str_parse_FASTA_protein());
			}
		}
		return con_v_c_analysis;
	}

	std::vector<c_genefamily_data> create_v_c_analysis_distinct(std::vector<c_genefamily_data>& par_v_c_analysis) {
		std::vector<c_genefamily_data> con_v_c_analysis_distinct;
		std::vector<string_type> con_v_str_genefamily_distinct;
		c_genefamily_data con_c_analysis_distinct;
		for (auto itr_v_c_analysis_data = par_v_c_analysis.begin(); itr_v_c_analysis_data != par_v_c_analysis.end(); ++itr_v_c_analysis_data) {
			string_type con_str_genefamily = itr_v_c_analysis_data->return_str_genefamily();
			size_type sw_peptideassociation_distinct = size_type();
			string_type con_str_genefamily_distinct = string_type();
			string_type con_str_genefamily_distinct_class = string_type();
			for (auto itr_str_genefamily = size_type(); itr_str_genefamily < con_str_genefamily.length(); ++itr_str_genefamily) {
				if (con_str_genefamily.at(itr_str_genefamily) == '*') {
					sw_peptideassociation_distinct = 1;
				}
				if (sw_peptideassociation_distinct == 0) {
					con_str_genefamily_distinct += con_str_genefamily.at(itr_str_genefamily);
				}
			}
			con_str_genefamily_distinct_class = itr_v_c_analysis_data->return_str_genefamily_class();
			if (std::find(con_v_str_genefamily_distinct.begin(), con_v_str_genefamily_distinct.end(), con_str_genefamily_distinct) == con_v_str_genefamily_distinct.end()) {
				con_v_str_genefamily_distinct.push_back(con_str_genefamily_distinct);
				con_c_analysis_distinct.set_str_genefamily(con_str_genefamily_distinct);
				con_c_analysis_distinct.set_str_genefamily_class(con_str_genefamily_distinct_class);
				con_c_analysis_distinct.str_species = itr_v_c_analysis_data->str_species;
				con_v_c_analysis_distinct.push_back(con_c_analysis_distinct);
			}
			for (auto itr_v_c_analysis_distinct = con_v_c_analysis_distinct.begin(); itr_v_c_analysis_distinct != con_v_c_analysis_distinct.end(); ++itr_v_c_analysis_distinct) {
				if (itr_v_c_analysis_distinct->return_str_genefamily() == con_str_genefamily_distinct) {
					itr_v_c_analysis_distinct->ref_v_c_analysis_polyassociation().push_back(&(*itr_v_c_analysis_data));
				}
			}
			sw_peptideassociation_distinct = 0;
			con_str_genefamily_distinct.clear();
		}
		return con_v_c_analysis_distinct;
	}

	void update_v_c_analysis_distinct(std::vector<c_genefamily_data>& par_v_c_analysis_distinct) {
		c_genefamily_data con_c_analysis_distinct_update;
		for (auto itr_v_c_analysis_distinct = par_v_c_analysis_distinct.begin(); itr_v_c_analysis_distinct != par_v_c_analysis_distinct.end(); ++itr_v_c_analysis_distinct) {
			std::vector<c_genefamily_data*> con_v_c_analysis_distinct_update = itr_v_c_analysis_distinct->ref_v_c_analysis_polyassociation();
			for (auto itr_v_c_analysis_distinct_polyassociation = con_v_c_analysis_distinct_update.begin(); itr_v_c_analysis_distinct_polyassociation != con_v_c_analysis_distinct_update.end(); ++itr_v_c_analysis_distinct_polyassociation) {
				if (itr_v_c_analysis_distinct_polyassociation == con_v_c_analysis_distinct_update.begin()) {
					con_c_analysis_distinct_update = **itr_v_c_analysis_distinct_polyassociation;
				}
				else {
					if (con_c_analysis_distinct_update.return_d_score() < (*itr_v_c_analysis_distinct_polyassociation)->return_d_score()) {
						con_c_analysis_distinct_update = **itr_v_c_analysis_distinct_polyassociation;
					}
				}
			}
			con_c_analysis_distinct_update.ref_v_c_analysis_polyassociation() = itr_v_c_analysis_distinct->ref_v_c_analysis_polyassociation();
			con_c_analysis_distinct_update.set_str_genefamily(itr_v_c_analysis_distinct->return_str_genefamily());
			con_c_analysis_distinct_update.set_str_genefamily_class(itr_v_c_analysis_distinct->return_str_genefamily_class());
			con_c_analysis_distinct_update.str_species = itr_v_c_analysis_distinct->str_species;
			*itr_v_c_analysis_distinct = con_c_analysis_distinct_update;
		}
	}

	std::vector<s_peptide_data> create_v_s_peptide_data(std::vector<fpf_parse::s_parse_peptides_csv> par_c_parse_csv_peptide_data) {
		std::vector<s_peptide_data> con_v_s_peptide_data;
		s_peptide_data con_s_peptide_data;
		for (auto itr_c_parse_csv_peptide_data = par_c_parse_csv_peptide_data.begin(); itr_c_parse_csv_peptide_data != par_c_parse_csv_peptide_data.end(); ++itr_c_parse_csv_peptide_data) {
			size_type ss_st_csv_IgP = size_type();
			std::istringstream(itr_c_parse_csv_peptide_data->str_parse_csv_IgP) >> ss_st_csv_IgP;
			//if (ss_st_csv_IgP >= PARSE_THRESHOLD_IgP) {
				con_s_peptide_data.str_peptide = itr_c_parse_csv_peptide_data->str_parse_peptides_csv_peptide;
				std::istringstream ss_spectralcount(itr_c_parse_csv_peptide_data->str_parse_peptides_csv_spectralcount);
				size_type ss_st_spectralcount;
				ss_spectralcount >> ss_st_spectralcount;
				con_s_peptide_data.st_spectralcount = ss_st_spectralcount;
				s_denovo_peptide con_s_denovo_peptide;
				for (auto i = 0; i < itr_c_parse_csv_peptide_data->str_parse_peptides_csv_peptide.size(); ++i) {
					con_s_denovo_peptide.ch_aminoacid = itr_c_parse_csv_peptide_data->str_parse_peptides_csv_peptide[i];
					con_s_denovo_peptide.d_denovo_localconfidence = itr_c_parse_csv_peptide_data->v_d_denovo_localconfidence[i];
					con_s_peptide_data.v_s_denovo_peptide.push_back(con_s_denovo_peptide);
				}
				con_v_s_peptide_data.push_back(con_s_peptide_data);
			//}
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

	void create_v_c_analysis_v_s_peptide_data(std::vector<c_genefamily_data>& par_v_c_analysis, std::vector<s_peptide_data> par_v_s_peptide_data) {
		string_type con_str_peptide = string_type();
		std::vector<s_peptide_data*> con_v_s_peptide_data;
		size_type sw_peptide_filtered = size_type();
		for (auto itr_v_c_analysis = par_v_c_analysis.begin(); itr_v_c_analysis != par_v_c_analysis.end(); ++itr_v_c_analysis) {
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
				if (itr_v_c_analysis->return_str_protein().find(con_str_peptide) != std::string::npos) {
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
			itr_v_c_analysis->set_v_s_peptide_data(con_v_s_peptide_data);
			con_v_s_peptide_data.clear();
			if ((IgFamily::DEBUG_MODE == 2) && (itr_v_c_analysis->return_v_s_peptide_data().size() != 0)) {
				std::cout << "\n\n" << itr_v_c_analysis->return_str_genefamily();
				std::cout << "   " << itr_v_c_analysis->return_str_protein();
				std::vector<s_peptide_data*> con_itr_v_s_peptide_data = itr_v_c_analysis->return_v_s_peptide_data();
				for (auto itr_v_s_peptide_data = con_itr_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_itr_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
					std::cout << "\n * " << (*itr_v_s_peptide_data)->str_peptide;
					std::cout << "   " << (*itr_v_s_peptide_data)->st_spectralcount;
				}
			}
		}
	}

	void create_v_c_analysis_v_s_peptide_data_distinct_filtered(std::vector<c_genefamily_data>& par_v_c_analysis, std::vector<s_peptide_data> par_v_data_filtered_distinct) {
		string_type con_str_peptide_distinct_filtered = string_type();
		std::vector<s_peptide_data*> con_v_s_peptide_data_distinct_filtered;
		for (auto itr_v_c_analysis = par_v_c_analysis.begin(); itr_v_c_analysis != par_v_c_analysis.end(); ++itr_v_c_analysis) {
			for (auto itr_v_data_filtered_distinct = par_v_data_filtered_distinct.begin(); itr_v_data_filtered_distinct != par_v_data_filtered_distinct.end(); ++itr_v_data_filtered_distinct) {
				if (itr_v_c_analysis->return_str_protein().find(itr_v_data_filtered_distinct->str_peptide) != std::string::npos) {
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
			itr_v_c_analysis->set_v_s_peptide_data_distinct_filtered(con_v_s_peptide_data_distinct_filtered);
			con_v_s_peptide_data_distinct_filtered.clear();
			if ((IgFamily::DEBUG_MODE == 2) && (itr_v_c_analysis->return_v_s_peptide_data_distinct_filtered().size() != 0)) {
				std::cout << "\n\n" << itr_v_c_analysis->return_str_genefamily();
				std::cout << "   " << itr_v_c_analysis->return_str_protein();
				std::vector<s_peptide_data*> con_itr_v_s_peptide_data_distinct_filtered = itr_v_c_analysis->return_v_s_peptide_data_distinct_filtered();
				for (auto itr_v_s_peptide_data_distinct_filtered = con_itr_v_s_peptide_data_distinct_filtered.begin(); itr_v_s_peptide_data_distinct_filtered != con_itr_v_s_peptide_data_distinct_filtered.end(); ++itr_v_s_peptide_data_distinct_filtered) {
					std::cout << "\n * " << (*itr_v_s_peptide_data_distinct_filtered)->str_peptide;
					std::cout << "   " << (*itr_v_s_peptide_data_distinct_filtered)->st_spectralcount;
				}
			}
		}
	}

	void create_v_c_analysis_str_alignment(std::vector<c_genefamily_data>& par_v_c_analysis) {
		for (auto itr_v_c_analysis = par_v_c_analysis.begin(); itr_v_c_analysis != par_v_c_analysis.end(); ++itr_v_c_analysis) {
			string_type con_str_alignment = string_type();
			con_str_alignment.clear();
			for (auto i = size_type(); i < itr_v_c_analysis->return_str_protein().length(); ++i) {
				con_str_alignment += '.';
			}
			std::vector<s_peptide_data*> con_itr_v_s_peptide_data_distinct_filtered = itr_v_c_analysis->return_v_s_peptide_data_distinct_filtered();
			for (auto itr_v_s_peptide_data_distinct_filtered = con_itr_v_s_peptide_data_distinct_filtered.begin(); itr_v_s_peptide_data_distinct_filtered != con_itr_v_s_peptide_data_distinct_filtered.end(); ++itr_v_s_peptide_data_distinct_filtered) {
				char con_str_alignment_match;
				char con_str_alignment_test;
				size_type sw_alignment = size_type();
				size_type count_str_alignment = size_type();
				string_type con_itr_str_protein = itr_v_c_analysis->return_str_protein();
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
			itr_v_c_analysis->set_str_alignment(con_str_alignment);
		}
	}

	void create_v_c_analysis_st_totalspectralcount(std::vector<c_genefamily_data>& par_v_c_analysis) {
		for (auto itr_v_c_analysis = par_v_c_analysis.begin(); itr_v_c_analysis != par_v_c_analysis.end(); ++itr_v_c_analysis) {
			size_type con_st_totalspectralcount = size_type();
			std::vector<s_peptide_data*> con_v_s_peptide_data = itr_v_c_analysis->return_v_s_peptide_data();
			for (auto itr_v_s_peptide_data = con_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
				con_st_totalspectralcount += (*itr_v_s_peptide_data)->st_spectralcount;
			}
			itr_v_c_analysis->set_st_totalspectralcount(con_st_totalspectralcount);
		}
	}

	void create_v_c_analysis_d_coverage(std::vector<c_genefamily_data>& par_v_c_analysis) {
		for (auto itr_v_c_analysis = par_v_c_analysis.begin(); itr_v_c_analysis != par_v_c_analysis.end(); ++itr_v_c_analysis) {
			double con_d_coverage;
			double count_str_alignment_true = double();
			for (auto itr_str_alignment = size_type(); itr_str_alignment < itr_v_c_analysis->return_str_alignment().length(); ++itr_str_alignment) {
				if (itr_v_c_analysis->return_str_alignment().at(itr_str_alignment) != '.') {
					++count_str_alignment_true;
				}
			}
			if ((itr_v_c_analysis->return_str_alignment().length()) == 0) {
				std::cout << "error: division by zero";
				string_type str_catch_error;
				std::cin >> str_catch_error;
			}
			con_d_coverage = 100 * (count_str_alignment_true / (itr_v_c_analysis->return_str_alignment().length()));
			itr_v_c_analysis->set_d_coverage(con_d_coverage);
		}
	}

	void create_v_c_peptide_v_p_peptideassociation(std::vector<c_genefamily_data>& par_v_c_analysis, std::vector<s_peptide_data>& par_v_s_peptide_data) {
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
			std::vector<fpf_parse::string_type> con_v_str_peptideassociation;
			for (auto itr_v_c_analysis = par_v_c_analysis.begin(); itr_v_c_analysis != par_v_c_analysis.end(); ++itr_v_c_analysis) {
				if (itr_v_c_analysis->return_str_protein().find(con_st_select_data) != std::string::npos) {
					con_v_str_peptideassociation.push_back(itr_v_c_analysis->return_str_genefamily());
					itr_v_s_peptide_data->v_p_peptideassociation.push_back(std::make_pair(&(*itr_v_c_analysis), double{ 0 }));
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

	void create_v_c_peptide_v_p_peptideassociation_distinct(std::vector<c_genefamily_data>& par_v_c_analysis_distinct, std::vector<s_peptide_data>& par_v_s_peptide_data) {
		std::vector<c_genefamily_data*> con_v_p_peptideassociation_distinct;
		std::vector<string_type> con_v_str_peptideassociation_distinct;
		for (auto itr_v_s_peptide_data = par_v_s_peptide_data.begin(); itr_v_s_peptide_data != par_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
			std::vector<std::pair<c_genefamily_data*, double>> con_v_p_peptideassociation = itr_v_s_peptide_data->v_p_peptideassociation;
			for (auto itr_v_p_peptideassociation = con_v_p_peptideassociation.begin(); itr_v_p_peptideassociation != con_v_p_peptideassociation.end(); ++itr_v_p_peptideassociation) {
				size_type sw_peptideassociation_distinct = size_type();
				string_type con_str_peptideassociation_distinct = string_type();
				for (auto i = size_type(); i < std::get<0>(*itr_v_p_peptideassociation)->return_str_genefamily().length(); ++i) {
					if (std::get<0>(*itr_v_p_peptideassociation)->return_str_genefamily().at(i) == '*') {
						sw_peptideassociation_distinct = 1;
					}
					if (sw_peptideassociation_distinct == 0) {
						con_str_peptideassociation_distinct += std::get<0>(*itr_v_p_peptideassociation)->return_str_genefamily().at(i);
					}
				}
				if (std::find(con_v_str_peptideassociation_distinct.begin(), con_v_str_peptideassociation_distinct.end(), con_str_peptideassociation_distinct) == con_v_str_peptideassociation_distinct.end()) {
					con_v_str_peptideassociation_distinct.push_back(con_str_peptideassociation_distinct);
					for (auto itr_v_c_analysis = par_v_c_analysis_distinct.begin(); itr_v_c_analysis != par_v_c_analysis_distinct.end(); ++itr_v_c_analysis) {
						if (itr_v_c_analysis->return_str_genefamily() == con_str_peptideassociation_distinct) {
							itr_v_s_peptide_data->v_p_peptideassociation_distinct.push_back(std::make_pair(&(*itr_v_c_analysis), double{ 1 }));
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

	inline void create_global_score_mean(std::vector<c_genefamily_data>& par_v_c_analysis) {
		double con_d_score = double();
		size_type con_st_nonzero_score = size_type();
		for (auto itr_v_c_analysis = par_v_c_analysis.begin(); itr_v_c_analysis != par_v_c_analysis.end(); ++itr_v_c_analysis) {
			if (itr_v_c_analysis->return_d_score() >= double{ 1.5 }) {
				con_d_score += itr_v_c_analysis->return_d_score();
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

	void create_v_c_analysis_d_score(std::vector<s_peptide_data>& par_v_s_peptide_data, std::vector<c_genefamily_data>& par_v_c_analysis, std::vector<c_genefamily_data>& par_v_c_analysis_distinct) {
		for (std::vector<c_genefamily_data>::iterator itr_v_c_analysis = par_v_c_analysis.begin(); itr_v_c_analysis != par_v_c_analysis.end(); ++itr_v_c_analysis) {
			string_type con_str_genefamily = itr_v_c_analysis->return_str_genefamily();
			size_type sw_peptideassociation_distinct = size_type();
			string_type con_str_genefamily_distinct = string_type();
			for (auto i = size_type(); i < con_str_genefamily.length(); ++i) {
				if (con_str_genefamily.at(i) == '*') {
					sw_peptideassociation_distinct = 1;
				}
				if (sw_peptideassociation_distinct == 0) {
					con_str_genefamily_distinct += con_str_genefamily.at(i);
				}
			}	
			std::vector<s_peptide_data*> con_v_s_peptide_data = itr_v_c_analysis->return_v_s_peptide_data();
			double con_d_score = double();
			double con_d_parameter_gene_family = double();
			double con_d_parameter_gene_family_train = double();
			double con_d_parameter_gene_family_weight = double();
			double con_d_parameter_IgP = double();
			for (auto itr_v_s_peptide_data = con_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
				con_d_parameter_gene_family_train = pow(((itr_v_c_analysis->return_d_score() + (double{ 100 } *IgFamily::SCORE_MEAN)) / (double{ 100 } *IgFamily::SCORE_MEAN)), double{ 1.6 });
				con_d_parameter_gene_family_train = 1;
				for (auto itr_v_p_peptideassociation_distinct = (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
					if (std::get<0>(*itr_v_p_peptideassociation_distinct)->return_str_genefamily() == con_str_genefamily_distinct) {
						con_d_parameter_gene_family_weight = std::get<1>(*itr_v_p_peptideassociation_distinct);
					}
				}
				con_d_parameter_gene_family = pow((con_d_parameter_gene_family_weight / (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.size()), double{ double{ 1 } / con_d_parameter_gene_family_train });
				con_d_parameter_IgP = (log_basechange(((double((*itr_v_s_peptide_data)->st_IgP) + double(IgFamily::PARSE_THRESHOLD_IgP)) / (double{ 2 } *double(IgFamily::PARSE_THRESHOLD_IgP))), double{ 9 }) + 1);
				con_d_score += (con_d_parameter_IgP * (*itr_v_s_peptide_data)->st_spectralcount * con_d_parameter_gene_family);
			}
			itr_v_c_analysis->set_d_score(con_d_score);
			con_str_genefamily_distinct.clear();
		}
		update_v_c_analysis_distinct(par_v_c_analysis_distinct);

		for (auto itr_v_s_peptide_data = par_v_s_peptide_data.begin(); itr_v_s_peptide_data != par_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {		
			double d_v_p_peptideassociation_d_score_sum = double();
			double d_v_p_peptideassociation_d_score_mean = double();
			for (auto itr_v_p_peptideassociation_distinct = itr_v_s_peptide_data->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != itr_v_s_peptide_data->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
				d_v_p_peptideassociation_d_score_sum += std::get<0>(*itr_v_p_peptideassociation_distinct)->return_d_score();
			}
			d_v_p_peptideassociation_d_score_mean = (d_v_p_peptideassociation_d_score_sum / itr_v_s_peptide_data->v_p_peptideassociation_distinct.size());
			for (auto itr_v_p_peptideassociation_distinct_2 = itr_v_s_peptide_data->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct_2 != itr_v_s_peptide_data->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct_2) {
				std::get<1>(*itr_v_p_peptideassociation_distinct_2) = (((std::get<0>(*itr_v_p_peptideassociation_distinct_2)->return_d_score() * double(1)) + (d_v_p_peptideassociation_d_score_mean)) / (double{ 2 } *d_v_p_peptideassociation_d_score_mean));
			}
		}
	}

	inline bool predicate_v_c_analysis_d_score(const c_genefamily_data& i, const c_genefamily_data& j) {
		return (i.return_d_score() > j.return_d_score());
	}

	inline void sort_v_c_analysis_d_score(std::vector<c_genefamily_data>& par_v_c_analysis) {
		std::sort(par_v_c_analysis.begin(), par_v_c_analysis.end(), predicate_v_c_analysis_d_score);
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

	std::vector<c_genefamily_data*> map_v_c_analysis_by_score(std::vector<c_genefamily_data>& par_v_c_analysis) {
		std::vector<c_genefamily_data*> con_v_c_analysis;
		c_genefamily_data* con_c_analysis;
		size_type st_count_map = size_type();
		while (st_count_map < par_v_c_analysis.size()) {
			double d_compare_map = double();
			for (std::vector<c_genefamily_data>::iterator itr_v_c_analysis = par_v_c_analysis.begin(); itr_v_c_analysis != par_v_c_analysis.end(); ++itr_v_c_analysis) {
				bool b_compare_map = ((std::find(con_v_c_analysis.begin(), con_v_c_analysis.end(), &(*itr_v_c_analysis))) == con_v_c_analysis.end());
				if ((d_compare_map <= itr_v_c_analysis->return_d_score()) && b_compare_map) {
					d_compare_map = itr_v_c_analysis->return_d_score();
					con_c_analysis = &(*itr_v_c_analysis);
				}
			}
			con_v_c_analysis.push_back(con_c_analysis);
			++st_count_map;
		}
		return con_v_c_analysis;
	}

	void train_v_c_analysis_d_score(std::vector<s_peptide_data>& par_v_s_peptide_data, std::vector<c_genefamily_data>& par_v_c_analysis, std::vector<c_genefamily_data>& par_v_c_analysis_distinct) {
		double delta_score_mean = double();
		for (IgFamily::GLOBAL_ITERATOR; IgFamily::GLOBAL_ITERATOR < IgFamily::ITERATE_TRAIN_SCORE; ++IgFamily::GLOBAL_ITERATOR) {
			delta_score_mean = IgFamily::SCORE_MEAN;
			create_global_score_mean(par_v_c_analysis_distinct);
			delta_score_mean -= IgFamily::SCORE_MEAN;
			std::cout << "\n\n iteration - " << IgFamily::GLOBAL_ITERATOR;
			std::cout << "   gene family score mean - " << std::fixed << std::setprecision(2) << IgFamily::SCORE_MEAN;
			if (IgFamily::GLOBAL_ITERATOR != 0) {
				std::cout << "   mean sample difference - " << std::fixed << std::setprecision(2) << abs(delta_score_mean);
			}
			fpf_data::create_v_c_analysis_d_score(par_v_s_peptide_data, par_v_c_analysis, par_v_c_analysis_distinct);
			fpf_data::create_v_c_analysis_v_s_peptide_data(par_v_c_analysis, par_v_s_peptide_data);
			//if (abs(delta_score_mean) < (IgFamily::SCORE_MEAN / double{ 10000 })) {
			//	break;
			//}
		}
		delta_score_mean = IgFamily::SCORE_MEAN;
		create_global_score_mean(par_v_c_analysis_distinct);
		delta_score_mean -= IgFamily::SCORE_MEAN;
		std::cout << "\n\n iteration - " << IgFamily::GLOBAL_ITERATOR;
		std::cout << "   gene family score mean - " << std::fixed << std::setprecision(2) << IgFamily::SCORE_MEAN;
		std::cout << "   mean sample difference - " << std::fixed << std::setprecision(2) << abs(delta_score_mean);
		fpf_data::create_v_c_analysis_v_s_peptide_data(par_v_c_analysis, par_v_s_peptide_data);
	}

	std::vector<c_genefamily_data*> map_v_c_analysis_by_genefamily(std::vector<c_genefamily_data>& par_v_c_analysis) {
		std::vector<c_genefamily_data*> con_v_c_analysis;
		c_genefamily_data* con_c_analysis;
		size_type st_count_map = size_type();
		while (st_count_map < par_v_c_analysis.size()) {
			string_type str_compare_map = string_type();
			for (auto itr_v_c_analysis = par_v_c_analysis.begin(); itr_v_c_analysis != par_v_c_analysis.end(); ++itr_v_c_analysis) {
				bool b_compare_map = ((std::find(con_v_c_analysis.begin(), con_v_c_analysis.end(), &(*itr_v_c_analysis))) == con_v_c_analysis.end());
				if ((str_compare_map < itr_v_c_analysis->return_str_genefamily()) && b_compare_map) {
					str_compare_map = itr_v_c_analysis->return_str_genefamily();
				}
			}
			con_v_c_analysis.push_back(con_c_analysis);
			str_compare_map.clear();
			++st_count_map;
		}
		return con_v_c_analysis;
	}

	void output(string_type par_str_fin, std::vector<c_genefamily_data*> par_v_c_analysis_distinct, std::vector<c_genefamily_data*> par_v_c_analysis) {
		
		std::cout << "\n\n ...streaming output";
		std::string output_v_c_analysis = par_str_fin + "_output.txt";
		std::ofstream fout_v_c_analysis;
		fout_v_c_analysis.open(output_v_c_analysis);
		fout_v_c_analysis << "-- IgFamily " << IgFamily::version << " --\n\n\n";
		fout_v_c_analysis << "Input file : " << IgFamily::INPUT_CSV;
		for (auto itr_v_c_analysis_distinct = par_v_c_analysis_distinct.begin(); itr_v_c_analysis_distinct != par_v_c_analysis_distinct.end(); ++itr_v_c_analysis_distinct) {
			if (itr_v_c_analysis_distinct == par_v_c_analysis_distinct.begin()) {
				fout_v_c_analysis << "\n\n\n\nTop gene families - ";
				fout_v_c_analysis << "     Score mean - " << std::fixed << std::setprecision(2) << IgFamily::SCORE_MEAN;
				fout_v_c_analysis << "\n\n";
			}
			if ((itr_v_c_analysis_distinct == par_v_c_analysis_distinct.end()) || ((*itr_v_c_analysis_distinct)->return_d_score() < double{ 0.01 })) {
				break;
			}
			fout_v_c_analysis << "\n  *  " << (*itr_v_c_analysis_distinct)->return_str_genefamily();
			fout_v_c_analysis << "   " << (*itr_v_c_analysis_distinct)->str_species;
			for (size_type itr_str_genefamily_output = (*itr_v_c_analysis_distinct)->return_str_genefamily().length() + (*itr_v_c_analysis_distinct)->str_species.length(); itr_str_genefamily_output < 35; ++itr_str_genefamily_output) {
				fout_v_c_analysis << " ";
			}
			fout_v_c_analysis << "   Score - " << std::fixed << std::setprecision(2) << (*itr_v_c_analysis_distinct)->return_d_score();
			std::vector<fpf_data::s_peptide_data*> con_itr_v_s_peptide_data = (*itr_v_c_analysis_distinct)->return_v_s_peptide_data();
			double d_peptide_1genefamily = double();
			double d_peptide_2genefamily = double();
			double d_peptide_3genefamily = double();
			for (auto itr_v_s_peptide_data = con_itr_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_itr_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
				string_type con_str_genefamily = (*itr_v_c_analysis_distinct)->return_str_genefamily();
				size_type sw_peptideassociation_distinct = size_type();
				string_type con_str_genefamily_distinct = string_type();
				for (auto i = size_type(); i < con_str_genefamily.length(); ++i) {
					if (con_str_genefamily.at(i) == '*') {
						sw_peptideassociation_distinct = 1;
					}
					if (sw_peptideassociation_distinct == 0) {
						con_str_genefamily_distinct += con_str_genefamily.at(i);
					}
				}
				double con_d_score = double();
				double con_d_parameter_gene_family = double();
				double con_d_parameter_gene_family_train = double();
				double con_d_parameter_gene_family_weight = double();
				double con_d_parameter_IgP = double();
				con_d_parameter_gene_family_train = pow((((*itr_v_c_analysis_distinct)->return_d_score() + (double{ 100 } *IgFamily::SCORE_MEAN)) / (double{ 100 } *IgFamily::SCORE_MEAN)), double{ 1.6 });
				con_d_parameter_gene_family_train = 1;
				for (auto itr_v_p_peptideassociation_distinct = (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
					if (std::get<0>(*itr_v_p_peptideassociation_distinct)->return_str_genefamily() == con_str_genefamily_distinct) {
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
				con_str_genefamily_distinct.clear();
			}
			d_peptide_1genefamily /= (*itr_v_c_analysis_distinct)->return_d_score();
			d_peptide_2genefamily /= (*itr_v_c_analysis_distinct)->return_d_score();
			d_peptide_3genefamily /= (*itr_v_c_analysis_distinct)->return_d_score();
			for (size_type itr_v_s_peptide_data_output = ((log10((*itr_v_c_analysis_distinct)->return_d_score()) >= double{ 1 }) ? (size_type(log10((*itr_v_c_analysis_distinct)->return_d_score())) + size_type(1)) : size_type(1)); itr_v_s_peptide_data_output < 6; ++itr_v_s_peptide_data_output) {
				fout_v_c_analysis << " ";
			}
			fout_v_c_analysis << "( 1GF - " << std::fixed << std::setprecision(2) << d_peptide_1genefamily;
			fout_v_c_analysis << ", 2GF - " << std::fixed << std::setprecision(2) << d_peptide_2genefamily;
			fout_v_c_analysis << ", 3GF - " << std::fixed << std::setprecision(2) << d_peptide_3genefamily;
			fout_v_c_analysis << " )";
			d_peptide_1genefamily = 0;
			d_peptide_2genefamily = 0;
			d_peptide_3genefamily = 0;
			//for (std::vector<fpf_data::c_genefamily_data*>::iterator itr_v_c_analysis_distinct_v_c_polyassociation = (itr_v_c_analysis_distinct + i)->ref_v_c_analysis_polyassociation().begin(); itr_v_c_analysis_distinct_v_c_polyassociation != (itr_v_c_analysis_distinct + i)->ref_v_c_analysis_polyassociation().end(); ++itr_v_c_analysis_distinct_v_c_polyassociation) {
			//	fout_v_c_analysis << "\n - - * " << (*itr_v_c_analysis_distinct_v_c_polyassociation)->return_str_genefamily();
			//	for (size_type itr_str_genefamily_output = (itr_v_c_analysis_distinct + i)->return_str_genefamily().length(); itr_str_genefamily_output < 20; ++itr_str_genefamily_output) {
			//		fout_v_c_analysis << " ";
			//	}
			//	fout_v_c_analysis << " - " << (*itr_v_c_analysis_distinct_v_c_polyassociation)->return_d_score();
			//}
		}

		//fout_v_c_analysis << "\n\n";
		//for (std::vector<fpf_data::c_genefamily_data*>::iterator itr_v_c_analysis_distinct = par_v_c_analysis_distinct.begin(); itr_v_c_analysis_distinct != par_v_c_analysis_distinct.end(); ++itr_v_c_analysis_distinct) {
		//	if ((*itr_v_c_analysis_distinct)->return_d_score() >= 5) {
		//		fout_v_c_analysis << "\n" << (*itr_v_c_analysis_distinct)->return_str_genefamily();
		//	}
		//}
		//fout_v_c_analysis << "\n\n";
		//for (std::vector<fpf_data::c_genefamily_data*>::iterator itr_v_c_analysis_distinct = par_v_c_analysis_distinct.begin(); itr_v_c_analysis_distinct != par_v_c_analysis_distinct.end(); ++itr_v_c_analysis_distinct) {
		//	if ((*itr_v_c_analysis_distinct)->return_d_score() >= 5) {
		//		fout_v_c_analysis << "\n" << (*itr_v_c_analysis_distinct)->return_d_score();
		//	}
		//}

		for (auto itr_v_c_analysis = par_v_c_analysis.begin(); itr_v_c_analysis != par_v_c_analysis.end(); ++itr_v_c_analysis) {
			sort_c_analysis_v_s_peptide_data_s_peptide((*itr_v_c_analysis)->return_v_s_peptide_data());
			IgFamily::SCORE_THRESHOLD = (IgFamily::SCORE_MEAN / double{ 2 });
			if (((*itr_v_c_analysis)->return_v_s_peptide_data().size() != 0) && ((*itr_v_c_analysis)->return_d_score() >= 3)) {
				fout_v_c_analysis << "\n\n\n\n\n" << (*itr_v_c_analysis)->return_str_genefamily();
				fout_v_c_analysis << "   " << (*itr_v_c_analysis)->str_species;
				fout_v_c_analysis << "   (top: ";
				double test2 = double();
				string_type test3 = string_type();
				for (auto test = (*itr_v_c_analysis)->ref_v_c_analysis_polyassociation().begin(); test != (*itr_v_c_analysis)->ref_v_c_analysis_polyassociation().end(); ++test) {
					if ((*test)->return_d_score() > test2) {
						test2 = (*test)->return_d_score();
						test3 = (*test)->return_str_genefamily();
					}
				}
				fout_v_c_analysis << test3 << ")";
				fout_v_c_analysis << "   #SC - " << (*itr_v_c_analysis)->return_st_totalspectralcount();
				fout_v_c_analysis << "   % - " << std::fixed << std::setprecision(2) << (*itr_v_c_analysis)->return_d_coverage();
				fout_v_c_analysis << "     Score - " << std::fixed << std::setprecision(2) << (*itr_v_c_analysis)->return_d_score();
				fout_v_c_analysis << "\n\n" << (*itr_v_c_analysis)->return_str_protein();
				fout_v_c_analysis << "\n" << (*itr_v_c_analysis)->return_str_alignment();
				fout_v_c_analysis << "\n";
				std::vector<fpf_data::s_peptide_data*> con_itr_v_s_peptide_data = (*itr_v_c_analysis)->return_v_s_peptide_data();
				for (auto itr_v_s_peptide_data = con_itr_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_itr_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
					//if ((*itr_v_s_peptide_data)->return_v_str_peptideassociation_distinct().size() <= 5) {
					fout_v_c_analysis << "\n * " << (*itr_v_s_peptide_data)->str_peptide;
					for (size_type itr_v_s_peptide_data_output = (*itr_v_s_peptide_data)->str_peptide.length(); itr_v_s_peptide_data_output < 50; ++itr_v_s_peptide_data_output) {
						fout_v_c_analysis << " ";
					}

					string_type con_s_genefamily = (*itr_v_c_analysis)->return_str_genefamily();
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
					//con_d_parameter_gene_family_train = pow((((*itr_v_c_analysis)->return_d_score() + (double{ 100 } *IgFamily::SCORE_MEAN)) / (double{ 100 } *IgFamily::SCORE_MEAN)), double{ 1.6 });
					con_d_parameter_gene_family_train = 1;
					for (auto itr_v_p_peptideassociation_distinct = (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
						if (std::get<0>(*itr_v_p_peptideassociation_distinct)->return_str_genefamily() == con_s_genefamily_distinct) {
							con_d_parameter_gene_family_weight = std::get<1>(*itr_v_p_peptideassociation_distinct);
						}
					}
					con_d_parameter_gene_family = pow((con_d_parameter_gene_family_weight / (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.size()), double{ double{ 1 } / con_d_parameter_gene_family_train });
					con_d_parameter_IgP = (log_basechange(((double((*itr_v_s_peptide_data)->st_IgP) + double(IgFamily::PARSE_THRESHOLD_IgP)) / (double{ 2 } *double(IgFamily::PARSE_THRESHOLD_IgP))), double{ 9 }) + 1);
					con_d_score += (con_d_parameter_IgP * (*itr_v_s_peptide_data)->st_spectralcount * con_d_parameter_gene_family);
					con_s_genefamily_distinct.clear();

					fout_v_c_analysis << "   Score - " << std::fixed << std::setprecision(2) << con_d_score;
					for (size_type itr_v_s_peptide_data_output = ((log10(con_d_score) >= double{ 1 }) ? (size_type(log10(con_d_score)) + size_type(1)) : size_type(1)); itr_v_s_peptide_data_output < 6; ++itr_v_s_peptide_data_output) {
						fout_v_c_analysis << " ";
					}
					//fout_v_c_analysis << "   n - " << std::fixed << std::setprecision(2) << con_d_score;
					//for (size_type itr_v_s_peptide_data_output = ((log10(con_d_score) >= double{ 1 }) ? (size_type(log10(con_d_score)) + size_type(1)) : size_type(1)); itr_v_s_peptide_data_output < 6; ++itr_v_s_peptide_data_output) {
					//	fout_v_c_analysis << " ";
					//}
					fout_v_c_analysis << "   SC - " << (*itr_v_s_peptide_data)->st_spectralcount;
					for (size_type itr_v_s_peptide_data_output = size_type(log10(double((*itr_v_s_peptide_data)->st_spectralcount))); itr_v_s_peptide_data_output < 2; ++itr_v_s_peptide_data_output) {
						fout_v_c_analysis << " ";
					}
					//fout_v_c_analysis << "   GF par - " << std::fixed << std::setprecision(2) << con_d_parameter_gene_family;
					//for (size_type itr_v_s_peptide_data_output = ((log10(con_d_score) >= double{ 1 }) ? size_type(log10(con_d_score)) : size_type(1)); itr_v_s_peptide_data_output < 4; ++itr_v_s_peptide_data_output) {
					//	fout_v_c_analysis << " ";
					//}
					//fout_v_c_analysis << "   Hom score - " << (*itr_v_s_peptide_data)->return_st_IgP() << "  ";
					//for (size_type itr_v_s_peptide_data_output = size_type(log10(double((*itr_v_s_peptide_data)->return_st_IgP()))); itr_v_s_peptide_data_output < 2; ++itr_v_s_peptide_data_output) {
					//	fout_v_c_analysis << " ";
					//}
					//fout_v_c_analysis << "   Train score - " << con_d_parameter_gene_family_train << "  ";
					//for (size_type itr_v_s_peptide_data_output = size_type(log10(con_d_parameter_gene_family_train)); itr_v_s_peptide_data_output < 3; ++itr_v_s_peptide_data_output) {
					//	fout_v_c_analysis << " ";
					//}
					//fout_v_c_analysis << "   Hom par - " << std::fixed << std::setprecision(2) << con_d_parameter_IgP;
					//for (size_type itr_v_s_peptide_data_output = size_type(log10(con_d_parameter_IgP)); itr_v_s_peptide_data_output < 3; ++itr_v_s_peptide_data_output) {
					//	fout_v_c_analysis << " ";
					//}
					fout_v_c_analysis << "   " << "Member of " << (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct.size() << " gene families - ";
					std::vector<std::pair<fpf_data::c_genefamily_data*, double>> con_itr_v_p_peptideassociation_distinct = (*itr_v_s_peptide_data)->v_p_peptideassociation_distinct;
					for (std::vector<std::pair<fpf_data::c_genefamily_data*, double>>::iterator itr_v_p_peptideassociation_distinct = con_itr_v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != con_itr_v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
						fout_v_c_analysis << std::get<0>(*itr_v_p_peptideassociation_distinct)->return_str_genefamily();
						fout_v_c_analysis << "(" << std::get<1>(*itr_v_p_peptideassociation_distinct) << ")";
						if ((itr_v_p_peptideassociation_distinct + 1) != con_itr_v_p_peptideassociation_distinct.end()) {
							fout_v_c_analysis << ", ";
						}
					}
				}
			}
		}
	}
}

#endif