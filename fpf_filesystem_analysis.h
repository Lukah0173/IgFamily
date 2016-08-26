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
	typedef fpf_filesystem::sample_analysis sample_analysis;
	typedef fpf_data::peptide_data peptide_data;

	struct filesystem_analysis;

	struct filesystem_analysis {
		fpf_filesystem::filesystem* filesystem;
		vector<fpf_data::peptide_data> v_replicate_peptide_data;
		vector<fpf_data::peptide_data> v_replicate_peptide_data_filtered;
	};

	void fout_filesystem_peptide_data_summary(vector<filesystem_analysis>& par_v_filesystem_analysis, vector<string> par_v_replicate_peptide_data_total_observed, string par_output) {
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
		for (auto itr_v_replicate_peptide_data_analysis = size_t(); itr_v_replicate_peptide_data_analysis < par_v_filesystem_analysis.begin()->v_replicate_peptide_data.size(); ++itr_v_replicate_peptide_data_analysis) {
			fout_v_filesystem_analysis_peptide_summary << par_v_replicate_peptide_data_total_observed[itr_v_replicate_peptide_data_analysis];
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