// * * fpf_filesystem_analysis.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_FILESYSTEM_ANALYSIS
#define	FPF_FILESYSTEM_ANALYSIS
#include <cstdlib> // provides - size_t
#include <vector> // provides - vector
#include <iostream> // provides - std::ofstream
#include <algorithm> // provides - std::find, std::find_if, std::sort
#include "IgFamily.h"
#include "fpf_filesystem.h"



namespace fpf_filesystem_analysis {

	using std::string;
	using std::vector;
	using std::pair;

	struct s_filesystem_analysis;

	struct s_filesystem_analysis {
		fpf_filesystem::filesystem* filesystem;
		vector<fpf_data::peptide_data> v_s_peptide_data_analysis;
		vector<fpf_data::peptide_data> v_s_peptide_data_filtered_distinct_analysis;
	};

	vector<string> create_v_str_peptide_total_observed(vector<fpf_filesystem::filesystem> par_v_s_filesystem) {
		vector<string> con_v_str_peptide_total_observed;
		for (auto itr_v_s_filesystem = par_v_s_filesystem.begin(); itr_v_s_filesystem != par_v_s_filesystem.end(); ++itr_v_s_filesystem) {
			auto con_v_s_peptide_data = itr_v_s_filesystem->v_peptide_data;
			for (auto itr_v_s_peptide_data = con_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
				if (std::find(con_v_str_peptide_total_observed.begin(), con_v_str_peptide_total_observed.end(), itr_v_s_peptide_data->peptide_withmod)
					== con_v_str_peptide_total_observed.end()) {
					con_v_str_peptide_total_observed.push_back(itr_v_s_peptide_data->peptide_withmod);
				}
			}
		}
		std::sort(con_v_str_peptide_total_observed.begin(), con_v_str_peptide_total_observed.end());
		return con_v_str_peptide_total_observed;
	}

	vector<string> create_v_str_peptide_filtered_distinct_total_observed(vector<fpf_filesystem::filesystem> par_v_s_filesystem) {
		vector<string> con_v_str_peptide_filtered_distinct_total_observed;
		for (auto itr_v_s_filesystem = par_v_s_filesystem.begin(); itr_v_s_filesystem != par_v_s_filesystem.end(); ++itr_v_s_filesystem) {
			auto con_v_s_peptide_data = itr_v_s_filesystem->v_peptide_data;
			for (auto itr_v_s_peptide_data : con_v_s_peptide_data) {
				if (std::find(con_v_str_peptide_filtered_distinct_total_observed.begin(), con_v_str_peptide_filtered_distinct_total_observed.end(), itr_v_s_peptide_data.peptide_filtered)
					== con_v_str_peptide_filtered_distinct_total_observed.end()) {
					con_v_str_peptide_filtered_distinct_total_observed.push_back(itr_v_s_peptide_data.peptide_filtered);
				}
			}
		}
		std::sort(con_v_str_peptide_filtered_distinct_total_observed.begin(), con_v_str_peptide_filtered_distinct_total_observed.end());
		return con_v_str_peptide_filtered_distinct_total_observed;
	}

	s_filesystem_analysis create_s_filesystem_analysis(fpf_filesystem::filesystem& par_filesystem, vector<string>& par_v_str_peptide_total_observed, vector<string>& par_v_str_peptide_filtered_distinct_total_observed) {
		s_filesystem_analysis con_s_filesystem_analysis;
		static fpf_data::peptide_data static_default_s_peptide_data = fpf_data::peptide_data();
		for (auto itr_v_str_peptide_total_observed = par_v_str_peptide_total_observed.begin(); itr_v_str_peptide_total_observed != par_v_str_peptide_total_observed.end(); ++itr_v_str_peptide_total_observed) {
			auto find_v_s_peptide_data = std::find_if(par_filesystem.v_peptide_data.begin(), par_filesystem.v_peptide_data.end(), 
				[itr_v_str_peptide_total_observed](fpf_data::peptide_data& par_s_peptide_data) {
				return par_s_peptide_data.peptide_withmod == *itr_v_str_peptide_total_observed; });
			if (find_v_s_peptide_data != par_filesystem.v_peptide_data.end()) {
				con_s_filesystem_analysis.v_s_peptide_data_analysis.push_back(*find_v_s_peptide_data);
			}
			else {
				con_s_filesystem_analysis.v_s_peptide_data_analysis.push_back(static_default_s_peptide_data);
			}
		}
		static fpf_data::peptide_data static_default_s_peptide_data_filtered_distinct = fpf_data::peptide_data();
		for (auto itr_v_str_peptide_filtered_distinct_total_observed : par_v_str_peptide_filtered_distinct_total_observed) {
			auto find_v_s_peptide_data_filtered_distinct = std::find_if(par_filesystem.v_peptide_data.begin(), par_filesystem.v_peptide_data.end(),
				[itr_v_str_peptide_filtered_distinct_total_observed](fpf_data::peptide_data& par_s_peptide_data_filtered_distinct) {
				return par_s_peptide_data_filtered_distinct.peptide_filtered == itr_v_str_peptide_filtered_distinct_total_observed; });
			if (find_v_s_peptide_data_filtered_distinct != par_filesystem.v_peptide_data.end()) {
				con_s_filesystem_analysis.v_s_peptide_data_filtered_distinct_analysis.push_back(*find_v_s_peptide_data_filtered_distinct);
			}
			else {
				con_s_filesystem_analysis.v_s_peptide_data_filtered_distinct_analysis.push_back(static_default_s_peptide_data_filtered_distinct);
			}
		}
		con_s_filesystem_analysis.filesystem = &par_filesystem;
		return con_s_filesystem_analysis;
	}

	void fout_filesystem_peptide_data_summary(vector<s_filesystem_analysis>& par_v_s_filesystem_analysis, vector<string> par_v_str_peptide_total_observed, string par_str_output) {
		std::string output_v_s_filesystem_analysis_peptide_summary = "summary_data\\" + par_str_output;
		std::ofstream fout_v_s_filesystem_analysis_peptide_summary;
		fout_v_s_filesystem_analysis_peptide_summary.open(output_v_s_filesystem_analysis_peptide_summary);
		
		fout_v_s_filesystem_analysis_peptide_summary << ",";
		for (auto itr_v_s_filesystem_analysis = par_v_s_filesystem_analysis.begin(); itr_v_s_filesystem_analysis != par_v_s_filesystem_analysis.end(); ++itr_v_s_filesystem_analysis) {
			fout_v_s_filesystem_analysis_peptide_summary << itr_v_s_filesystem_analysis->filesystem->filename << ",";			
		}
		fout_v_s_filesystem_analysis_peptide_summary << std::endl;
		for (auto itr_v_s_filesystem_analysis = par_v_s_filesystem_analysis.begin(); itr_v_s_filesystem_analysis != par_v_s_filesystem_analysis.end(); ++itr_v_s_filesystem_analysis) {
			fout_v_s_filesystem_analysis_peptide_summary << itr_v_s_filesystem_analysis->filesystem->patientstatus << ",";
		}
		fout_v_s_filesystem_analysis_peptide_summary << std::endl;
		for (auto itr_v_s_filesystem_analysis = par_v_s_filesystem_analysis.begin(); itr_v_s_filesystem_analysis != par_v_s_filesystem_analysis.end(); ++itr_v_s_filesystem_analysis) {
			fout_v_s_filesystem_analysis_peptide_summary << itr_v_s_filesystem_analysis->filesystem->filesystem_replicate_count << ",";
		}
		fout_v_s_filesystem_analysis_peptide_summary << std::endl;
		for (auto itr_v_s_peptide_data_analysis = size_t(); itr_v_s_peptide_data_analysis < par_v_s_filesystem_analysis.begin()->v_s_peptide_data_analysis.size(); ++itr_v_s_peptide_data_analysis) {
			fout_v_s_filesystem_analysis_peptide_summary << par_v_str_peptide_total_observed[itr_v_s_peptide_data_analysis];
			fout_v_s_filesystem_analysis_peptide_summary << ",";
			for (auto itr_v_s_filesystem_analysis = par_v_s_filesystem_analysis.begin(); itr_v_s_filesystem_analysis != par_v_s_filesystem_analysis.end(); ++itr_v_s_filesystem_analysis) {				
				fout_v_s_filesystem_analysis_peptide_summary << itr_v_s_filesystem_analysis->v_s_peptide_data_analysis[itr_v_s_peptide_data_analysis].spectralcount;
				fout_v_s_filesystem_analysis_peptide_summary << ",";
			}
			fout_v_s_filesystem_analysis_peptide_summary << std::endl;
		}
	}

	void fout_filesystem_summary(vector<s_filesystem_analysis>& par_v_s_filesystem_analysis) {
		std::string output_v_s_filesystem_summary = "summary_data\\filesystem_summary.txt";
		std::ofstream fout_v_s_filesystem_summary;
		fout_v_s_filesystem_summary.open(output_v_s_filesystem_summary);

		for (auto itr_v_s_filesystem_analysis = par_v_s_filesystem_analysis.begin(); itr_v_s_filesystem_analysis != par_v_s_filesystem_analysis.end(); ++itr_v_s_filesystem_analysis) {
			fout_v_s_filesystem_summary << "ID: " << std::get<0>(itr_v_s_filesystem_analysis->filesystem->filesystem_id) << "," << std::get<1>(itr_v_s_filesystem_analysis->filesystem->filesystem_id) << ";\n";
			fout_v_s_filesystem_summary << "FILE: " << itr_v_s_filesystem_analysis->filesystem->filename << ";\n";
			fout_v_s_filesystem_summary << "VERSION: " << IgFamily::version << ";\n";
			fout_v_s_filesystem_summary << "REPLICATES: ";
			for (vector<pair<string, string>>::iterator itr_v_p_replicates = itr_v_s_filesystem_analysis->filesystem->v_filesystem_replicates.begin(); itr_v_p_replicates != itr_v_s_filesystem_analysis->filesystem->v_filesystem_replicates.end(); ++itr_v_p_replicates) {
				fout_v_s_filesystem_summary << std::get<0>(*itr_v_p_replicates) << "," << std::get<1>(*itr_v_p_replicates);
				if ((itr_v_p_replicates + 1) != itr_v_s_filesystem_analysis->filesystem->v_filesystem_replicates.end()) {
					fout_v_s_filesystem_summary << ",";
				}
			}
			fout_v_s_filesystem_summary << ";\n";
			fout_v_s_filesystem_summary << "STATUS: " << itr_v_s_filesystem_analysis->filesystem->patientstatus << ";";
			fout_v_s_filesystem_summary << std::endl << std::endl;
		}
	}

}

#endif