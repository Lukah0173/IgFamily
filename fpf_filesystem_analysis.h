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
	
	typedef fpf_filesystem::filesystem filesystem;
	typedef fpf_data::peptide_data peptide_data;

	struct filesystem_analysis;

	struct filesystem_analysis {
		fpf_filesystem::filesystem* filesystem;
		vector<fpf_data::peptide_data> v_peptide_data;
		vector<fpf_data::peptide_data> v_peptide_data_filtered;
	};

	vector<string> create_v_peptide_data_total_observed(vector<filesystem> par_v_filesystem) {
		vector<string> temp_v_peptide_data_total_observed = vector<string>();
		for (const auto& itr_v_filesystem : par_v_filesystem) { 
			const auto& temp_v_peptide_data = itr_v_filesystem.v_peptide_data;
			for (const auto& itr_v_peptide_data : temp_v_peptide_data) {
				if (std::find(temp_v_peptide_data_total_observed.begin(), temp_v_peptide_data_total_observed.end(), itr_v_peptide_data.peptide_withmod)
					== temp_v_peptide_data_total_observed.end()) {
					temp_v_peptide_data_total_observed.push_back(itr_v_peptide_data.peptide_withmod);
				}
			}
		}
		std::sort(temp_v_peptide_data_total_observed.begin(), temp_v_peptide_data_total_observed.end());
		return temp_v_peptide_data_total_observed;
	}

	vector<string> create_v_peptide_data_filtered_total_observed(vector<filesystem> par_v_filesystem) {
		vector<string> temp_v_peptide_data_total_observed = vector<string>();
		for (const auto& itr_v_filesystem : par_v_filesystem) {
			const auto& temp_v_peptide_data = itr_v_filesystem.v_peptide_data;
			for (const auto itr_v_peptide_data : temp_v_peptide_data) {
				if (std::find(temp_v_peptide_data_total_observed.begin(), temp_v_peptide_data_total_observed.end(), itr_v_peptide_data.peptide_withmod)
					== temp_v_peptide_data_total_observed.end()) {
					temp_v_peptide_data_total_observed.push_back(itr_v_peptide_data.peptide_filtered);
				}
			}
		}
		std::sort(temp_v_peptide_data_total_observed.begin(), temp_v_peptide_data_total_observed.end());
		return temp_v_peptide_data_total_observed;
	}

	filesystem_analysis create_filesystem_analysis(filesystem& par_filesystem, vector<string>& par_v_peptide_data_total_observed, vector<string>& par_v_peptide_data_filtered_total_observed) {
		filesystem_analysis temp_filesystem_analysis = filesystem_analysis();
		peptide_data temp_peptide_data = peptide_data();
		for (const auto itr_v_peptide_total_observed : par_v_peptide_data_total_observed) {
			const auto find_v_peptide_data = std::find_if(par_filesystem.v_peptide_data.begin(), par_filesystem.v_peptide_data.end(), 
				[itr_v_peptide_total_observed](fpf_data::peptide_data par_peptide_data) {
				return par_peptide_data.peptide_withmod == itr_v_peptide_total_observed; });
			if (find_v_peptide_data != par_filesystem.v_peptide_data.end()) {
				temp_filesystem_analysis.v_peptide_data.push_back(*find_v_peptide_data);
			}
			else {
				temp_filesystem_analysis.v_peptide_data.push_back(temp_peptide_data);
			}
		}
		peptide_data temp_peptide_data_filtered = peptide_data();
		for (const auto itr_v_peptide_total_observed : par_v_peptide_data_filtered_total_observed) {
			const auto find_v_peptide_data = std::find_if(par_filesystem.v_peptide_data.begin(), par_filesystem.v_peptide_data.end(),
				[itr_v_peptide_total_observed](fpf_data::peptide_data par_peptide_data) {
				return par_peptide_data.peptide_withmod == itr_v_peptide_total_observed; });
			if (find_v_peptide_data != par_filesystem.v_peptide_data.end()) {
				temp_filesystem_analysis.v_peptide_data.push_back(*find_v_peptide_data);
			}
			else {
				temp_filesystem_analysis.v_peptide_data.push_back(temp_peptide_data);
			}
		}
		temp_filesystem_analysis.filesystem = &par_filesystem;
		return temp_filesystem_analysis;
	}

	void fout_filesystem_peptide_data_summary(vector<filesystem_analysis>& par_v_filesystem_analysis, vector<string> par_v_peptide_data_total_observed, string par_output) {
		std::string output_v_filesystem_analysis_peptide_summary = "summary_data\\" + par_output;
		std::ofstream fout_v_filesystem_analysis_peptide_summary;
		fout_v_filesystem_analysis_peptide_summary.open(output_v_filesystem_analysis_peptide_summary);
		
		fout_v_filesystem_analysis_peptide_summary << ",";
		for (const auto& itr_v_filesystem_analysis : par_v_filesystem_analysis) {
			fout_v_filesystem_analysis_peptide_summary << itr_v_filesystem_analysis.filesystem->filename << ",";			
		}
		fout_v_filesystem_analysis_peptide_summary << std::endl;
		for (const auto& itr_v_filesystem_analysis : par_v_filesystem_analysis) {
			fout_v_filesystem_analysis_peptide_summary << itr_v_filesystem_analysis.filesystem->patientstatus << ",";
		}
		fout_v_filesystem_analysis_peptide_summary << std::endl;
		for (const auto& itr_v_filesystem_analysis : par_v_filesystem_analysis) {
			fout_v_filesystem_analysis_peptide_summary << itr_v_filesystem_analysis.filesystem->filesystem_replicate_count << ",";
		}
		fout_v_filesystem_analysis_peptide_summary << std::endl;
		for (auto itr_v_peptide_data_analysis = size_t(); itr_v_peptide_data_analysis < par_v_filesystem_analysis.begin()->v_peptide_data.size(); ++itr_v_peptide_data_analysis) {
			fout_v_filesystem_analysis_peptide_summary << par_v_peptide_data_total_observed[itr_v_peptide_data_analysis];
			fout_v_filesystem_analysis_peptide_summary << ",";
			fout_v_filesystem_analysis_peptide_summary << std::endl;
		}
	}

	void fout_filesystem_summary(vector<filesystem_analysis>& par_v_filesystem_analysis) {
		std::string output_v_filesystem_summary = "summary_data\\filesystem_summary.txt";
		std::ofstream fout_v_filesystem_summary;
		fout_v_filesystem_summary.open(output_v_filesystem_summary);

		for (auto itr_v_filesystem_analysis = par_v_filesystem_analysis.begin(); itr_v_filesystem_analysis != par_v_filesystem_analysis.end(); ++itr_v_filesystem_analysis) {
			fout_v_filesystem_summary << "ID: " << std::get<0>(itr_v_filesystem_analysis->filesystem->filesystem_id) << "," << std::get<1>(itr_v_filesystem_analysis->filesystem->filesystem_id) << ";\n";
			fout_v_filesystem_summary << "FILE: " << itr_v_filesystem_analysis->filesystem->filename << ";\n";
			fout_v_filesystem_summary << "VERSION: " << IgFamily::version << ";\n";
			fout_v_filesystem_summary << "REPLICATES: ";
			for (auto itr_v_p_replicates = itr_v_filesystem_analysis->filesystem->v_filesystem_replicates.begin(); itr_v_p_replicates != itr_v_filesystem_analysis->filesystem->v_filesystem_replicates.end(); ++itr_v_p_replicates) {
				fout_v_filesystem_summary << std::get<0>(*itr_v_p_replicates) << "," << std::get<1>(*itr_v_p_replicates);
				if ((itr_v_p_replicates + 1) != itr_v_filesystem_analysis->filesystem->v_filesystem_replicates.end()) {
					fout_v_filesystem_summary << ",";
				}
			}
			fout_v_filesystem_summary << ";\n";
			fout_v_filesystem_summary << "STATUS: " << itr_v_filesystem_analysis->filesystem->patientstatus << ";";
			fout_v_filesystem_summary << std::endl << std::endl;
		}
	}

}

#endif