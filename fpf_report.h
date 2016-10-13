// * * fpf_report.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_REPORT
#define	FPF_REPORT

#include <cstdlib>
#include <iomanip>
#include <map>
#include <math.h>
#include <string>
#include <utility>
#include <vector>

#include "IgFamily.h"
#include "fpf_data.h"
#include "fpf_filesystem.h"
#include "fpf_multinomial.h"


namespace fpf_report {

	using std::map;
	using std::multimap;
	using std::pair;
	using std::string;
	using std::vector;

	using fpf_data::multinomial;
	using fpf_data::multinomial_frequency_type;
	using fpf_data::peptide_analysis;
	using fpf_filesystem::filesystem;
	using fpf_filesystem::sample_analysis;

	void fout_v_peptide_data(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::string output_v_peptide_data = par_filesystem.directory + par_filesystem.filename + "_peptide_data_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_v_peptide_data;
		fout_v_peptide_data.open(output_v_peptide_data);
		fout_v_peptide_data << "key,scan_ID,peptide_mz,peptide_z,peptide_rt,peptide_m,peptide_withmod,peptide_withoutmod,peptide_filtered,denovo_peptide,denovo_peptide_filtered\n";
		for (const auto& itr_v_peptide_data_map : par_sample_analysis.v_peptide_data_map) {
			fout_v_peptide_data << itr_v_peptide_data_map.second->key_peptide_data << ",";
			fout_v_peptide_data << itr_v_peptide_data_map.second->scan_ID << ",";
			fout_v_peptide_data << itr_v_peptide_data_map.second->peptide_mz << ",";
			fout_v_peptide_data << itr_v_peptide_data_map.second->peptide_z << ",";
			fout_v_peptide_data << itr_v_peptide_data_map.second->peptide_rt << ",";
			fout_v_peptide_data << itr_v_peptide_data_map.second->peptide_m << ",";
			fout_v_peptide_data << itr_v_peptide_data_map.second->peptide_withmod << ",";
			fout_v_peptide_data << itr_v_peptide_data_map.second->peptide_withoutmod << ",";
			fout_v_peptide_data << itr_v_peptide_data_map.second->peptide_filtered << ",";
			for (const auto& itr_v_denovo_aminoacid : itr_v_peptide_data_map.second->denovo_peptide_data.v_denovo_aminoacid) {
				fout_v_peptide_data << itr_v_denovo_aminoacid.aminoacid;
				fout_v_peptide_data << "[" << itr_v_denovo_aminoacid.aminoacid_localconfidence << "]";
			}
			fout_v_peptide_data << ",";
			for (const auto& itr_v_denovo_filtered_aminoacid : itr_v_peptide_data_map.second->denovo_peptide_data_filtered.v_denovo_aminoacid) {
				fout_v_peptide_data << itr_v_denovo_filtered_aminoacid.aminoacid;
				fout_v_peptide_data << "[" << itr_v_denovo_filtered_aminoacid.aminoacid_localconfidence << "]";
			}
			fout_v_peptide_data << ",\n";
		}
	}

	void fout_v_protein_data(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::string output_v_protein_data = par_filesystem.directory + par_filesystem.filename + "_protein_data_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_v_protein_data_map;
		fout_v_protein_data_map.open(output_v_protein_data);
		fout_v_protein_data_map << "key,protein_name,protein_type,protein_species,protein_protein,\n";
		for (const auto& itr_v_protein_data_map : par_sample_analysis.v_protein_data_map) {
			fout_v_protein_data_map << itr_v_protein_data_map.second->key_protein_data << ",";
			fout_v_protein_data_map << itr_v_protein_data_map.second->protein_name << ",";
			fout_v_protein_data_map << itr_v_protein_data_map.second->protein_type << ",";
			fout_v_protein_data_map << itr_v_protein_data_map.second->protein_species << ",";
			fout_v_protein_data_map << itr_v_protein_data_map.second->protein_protein << ",";
			fout_v_protein_data_map << "\n";
		}
	}

	void fout_v_peptide_analysis(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::string output_v_peptide_analysis = par_filesystem.directory + par_filesystem.filename + "_peptide_analysis_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_v_peptide_analysis;
		fout_v_peptide_analysis.open(output_v_peptide_analysis);
		fout_v_peptide_analysis << "key,peptide_filtered,replicate_count,v_peptide_withoutmod_mz,v_peptide_withoutmod_z,v_peptide_withoutmod_rt,v_peptide_withoutmod_m,\n";
		for (const auto& itr_v_peptide_anaylsis_map : par_sample_analysis.v_peptide_analysis_map) {
			fout_v_peptide_analysis << itr_v_peptide_anaylsis_map.second->key_peptide_analysis << ",";
			fout_v_peptide_analysis << itr_v_peptide_anaylsis_map.second->peptide_filtered << ",";
			fout_v_peptide_analysis << itr_v_peptide_anaylsis_map.second->replicate_count << ",";
			for (const auto& itr_v_peptide_withoutmod_mz : itr_v_peptide_anaylsis_map.second->v_peptide_withoutmod_mz) {
				fout_v_peptide_analysis << itr_v_peptide_withoutmod_mz << "|";
			}
			fout_v_peptide_analysis << ",";
			for (const auto& itr_v_peptide_withoutmod_z : itr_v_peptide_anaylsis_map.second->v_peptide_withoutmod_z) {
				fout_v_peptide_analysis << itr_v_peptide_withoutmod_z << "|";
			}
			fout_v_peptide_analysis << ",";
			for (const auto& itr_v_peptide_withoutmod_rt : itr_v_peptide_anaylsis_map.second->v_peptide_withoutmod_rt) {
				fout_v_peptide_analysis << itr_v_peptide_withoutmod_rt << "|";
			}
			fout_v_peptide_analysis << ",";
			for (const auto& itr_v_peptide_withoutmod_m : itr_v_peptide_anaylsis_map.second->v_peptide_withoutmod_m) {
				fout_v_peptide_analysis << itr_v_peptide_withoutmod_m << "|";
			}
			fout_v_peptide_analysis << ",\n";
		}
	}

	void fout_v_protein_analysis(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::string output_v_protein_analysis = par_filesystem.directory + par_filesystem.filename + "_protein_analysis_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_v_protein_analysis;
		fout_v_protein_analysis.open(output_v_protein_analysis);
		fout_v_protein_analysis << "key,protein_name,protein_score,protein_effective_spectral_count,proteinconstruct_sequencecoverage,protein_type,blastp_query,denovo_replicate_count,score_density,blastp_homology,blastp_homology_density,blastp_homology_density_conjugated,blastp_mismatch_count,alignment_coverage_delta,\n";
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			for (const auto& itr_v_homology_data_combined_by_protein : itr_v_protein_analysis.v_homology_data_combined_by_protein) {
				if (itr_v_homology_data_combined_by_protein.blastp_homology_density_conjugated >= IgFamily::REPORT_QUERY_PARAMETER_DENSITY_CONJUGATED_THRESHOLD) {
					fout_v_protein_analysis << itr_v_protein_analysis.key_protein_analysis << ",";
					fout_v_protein_analysis << itr_v_protein_analysis.p_protein_data->protein_name << ",";
					fout_v_protein_analysis << itr_v_protein_analysis.protein_score << ",";
					fout_v_protein_analysis << itr_v_protein_analysis.protein_effective_spectral_count << ",";
					fout_v_protein_analysis << itr_v_protein_analysis.proteinconstruct_sequencecoverage << ",";
					fout_v_protein_analysis << itr_v_protein_analysis.p_protein_data->protein_type << ",";
					fout_v_protein_analysis << itr_v_homology_data_combined_by_protein.blastp_query << ",";
					fout_v_protein_analysis << itr_v_homology_data_combined_by_protein.denovo_replicate_count << ",";
					fout_v_protein_analysis << itr_v_homology_data_combined_by_protein.score_density << ",";
					fout_v_protein_analysis << itr_v_homology_data_combined_by_protein.blastp_homology << ",";
					fout_v_protein_analysis << itr_v_homology_data_combined_by_protein.blastp_homology_density << ",";
					fout_v_protein_analysis << itr_v_homology_data_combined_by_protein.blastp_homology_density_conjugated << ",";
					fout_v_protein_analysis << itr_v_homology_data_combined_by_protein.blastp_mismatch_count << ",";
					fout_v_protein_analysis << itr_v_homology_data_combined_by_protein.alignment_coverage_delta << ",";
					fout_v_protein_analysis << ",\n";
				}
			}
		}
	}

	void fout_v_homology_data(const filesystem& par_filesystem, const sample_analysis& par_sample_analysis) {
		string output_blastp_summary = par_filesystem.directory + par_filesystem.filename;
		output_blastp_summary += "_blastp_summary.csv";
		std::ofstream fout_v_homology_data;
		fout_v_homology_data.open(output_blastp_summary);
		fout_v_homology_data << "key_query,";
		fout_v_homology_data << "query,";
		fout_v_homology_data << "key_subject_accession,";
		fout_v_homology_data << "subject_accession,";
		fout_v_homology_data << "mismatch_count,";
		fout_v_homology_data << "query_alignment_coverage,";
		fout_v_homology_data << "query_alignment_coverage_delta,";
		fout_v_homology_data << "score,";
		fout_v_homology_data << "score_transformed,";
		fout_v_homology_data << "par_dens,";
		fout_v_homology_data << "par_score,";
		fout_v_homology_data << "\n";
		for (auto itr_v_homology_data : par_sample_analysis.v_homology_data) {
			fout_v_homology_data << itr_v_homology_data.key_blastp_query << ",";
			fout_v_homology_data << itr_v_homology_data.blastp_query << ",";
			fout_v_homology_data << itr_v_homology_data.key_blastp_subject_accession << ",";
			fout_v_homology_data << itr_v_homology_data.blastp_subject_accession << ",";
			fout_v_homology_data << itr_v_homology_data.blastp_mismatch_count << ",";
			fout_v_homology_data << itr_v_homology_data.alignment_coverage << ",";
			fout_v_homology_data << itr_v_homology_data.alignment_coverage_delta << ",";
			fout_v_homology_data << itr_v_homology_data.blastp_homology << ",";
			fout_v_homology_data << itr_v_homology_data.blastp_homology_transformed_conjugated << ",";
			fout_v_homology_data << itr_v_homology_data.blastp_homology_density_conjugated << ",";
			fout_v_homology_data << itr_v_homology_data.score << ",";
			fout_v_homology_data << "\n";
		}
	}

	void fout_multinomial(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		string output_multinomial = par_filesystem.directory + par_filesystem.filename + "_multinomial_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_multinomial;
		fout_multinomial.open(output_multinomial);
		fout_multinomial << ",";
		fout_multinomial << "TOTAL,";
		for (const auto& itr_protein_data : par_sample_analysis.multinomial_data.v_p_protein_data) {
			fout_multinomial << itr_protein_data->protein_name << ",";
		}
		fout_multinomial << "\n";
		fout_multinomial << ",,";
		for (const auto& itr_protein_data : par_sample_analysis.multinomial_data.v_p_protein_data) {
			fout_multinomial << itr_protein_data->protein_type << ",";
		}
		fout_multinomial << "\n";
		for (auto i = 0; i < par_sample_analysis.multinomial_data.v_p_peptide_analysis.size(); ++i) {
			fout_multinomial << par_sample_analysis.multinomial_data.v_p_peptide_analysis[i] << ",";
			fout_multinomial << par_sample_analysis.multinomial_data.v_frequency_marginal_sum[i] << ",";
			for (auto j = 0; j < par_sample_analysis.multinomial_data.v_p_protein_data.size(); ++j) {
				fout_multinomial << par_sample_analysis.multinomial_data.v2_frequency[i][j] << ",";
			}
			fout_multinomial << "\n";
		}
	}

	void fout_multinomial_element(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		string output_multinomial_element = par_filesystem.directory + par_filesystem.filename + "_multinomial_peptide_" + par_sample_analysis.peptide_assignment_method + ".txt";
		std::ofstream fout_multinomial_element;
		fout_multinomial_element.open(output_multinomial_element);
		fout_multinomial_element << "\n";
		fout_multinomial_element << par_filesystem.filename;
		fout_multinomial_element << "\n\n\n";
		size_t format_ws_length{};
		for (auto i = 0; i < par_sample_analysis.multinomial_data.v_p_peptide_analysis.size(); ++i) {
			if (par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->peptide_filtered.length() > format_ws_length) {
				format_ws_length = par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->peptide_filtered.length();
			}
		}
		for (auto i = 0; i < par_sample_analysis.multinomial_data.v_p_peptide_analysis.size(); ++i) {
			fout_multinomial_element << par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->peptide_filtered;
			for (auto j = 0; j < (format_ws_length - par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->peptide_filtered.length() + 5); ++j) {
				fout_multinomial_element << " ";
			}
			vector<multinomial_frequency_type> temp_v_multinomial_frequency{};
			for (auto j = 0; j < par_sample_analysis.multinomial_data.v_p_protein_data.size(); ++j) {
				multinomial_frequency_type temp_multinomial_frequency{};
				temp_multinomial_frequency.p_protein_data = par_sample_analysis.multinomial_data.v_p_protein_data[j];
				temp_multinomial_frequency.multinomial_frequency = par_sample_analysis.multinomial_data.v2_frequency[i][j];
				temp_v_multinomial_frequency.push_back(temp_multinomial_frequency);
			}
			fpf_multinomial::sort_v_multinomial_frequency(temp_v_multinomial_frequency);
			for (auto j = 0; j < par_sample_analysis.multinomial_data.v_p_protein_data.size(); ++j) {
				if (temp_v_multinomial_frequency[j].multinomial_frequency > 0.2) {
					fout_multinomial_element << temp_v_multinomial_frequency[j].p_protein_data->protein_name << " (";
					fout_multinomial_element << temp_v_multinomial_frequency[j].multinomial_frequency << "), ";
				}
			}
			fout_multinomial_element << "\n";
		}
	}

	void fout_multinomial_element_nomatch(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		string output_multinomial_element_nomatch = par_filesystem.directory + par_filesystem.filename + "_multinomial_peptide_filtered_" + par_sample_analysis.peptide_assignment_method + ".txt";
		std::ofstream fout_multinomial_element_nomatch;
		fout_multinomial_element_nomatch.open(output_multinomial_element_nomatch);
		fout_multinomial_element_nomatch << "\n";
		fout_multinomial_element_nomatch << par_filesystem.filename;
		fout_multinomial_element_nomatch << "\n\n\n";
		size_t format_ws_length{};
		for (auto i = 0; i < par_sample_analysis.multinomial_data.v_p_peptide_analysis.size(); ++i) {
			if (par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->peptide_filtered.length() > format_ws_length) {
				format_ws_length = par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->peptide_filtered.length();
			}
		}
		for (auto i = 0; i < par_sample_analysis.multinomial_data.v_p_peptide_analysis.size(); ++i) {
			double format_frequency_threshold = double();
			for (auto j = 0; j < par_sample_analysis.multinomial_data.v_p_protein_data.size(); ++j) {
				if (par_sample_analysis.multinomial_data.v2_frequency[i][j] > format_frequency_threshold) {
					format_frequency_threshold = par_sample_analysis.multinomial_data.v2_frequency[i][j];
				}
			}
			if (format_frequency_threshold < 0.1) {
				fout_multinomial_element_nomatch << par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->peptide_filtered;
				for (auto j = 0; j < (format_ws_length - par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->peptide_filtered.length() + 5); ++j) {
					fout_multinomial_element_nomatch << " ";
				}
				vector<multinomial_frequency_type> temp_v_multinomial_frequency{};
				for (auto j = 0; j < par_sample_analysis.multinomial_data.v_p_protein_data.size(); ++j) {
					multinomial_frequency_type temp_multinomial_frequency{};
					temp_multinomial_frequency.p_protein_data = par_sample_analysis.multinomial_data.v_p_protein_data[j];
					temp_multinomial_frequency.multinomial_frequency = par_sample_analysis.multinomial_data.v2_frequency[i][j];
					temp_v_multinomial_frequency.push_back(temp_multinomial_frequency);
				}
				fpf_multinomial::sort_v_multinomial_frequency(temp_v_multinomial_frequency);
				for (auto j = 0; j < par_sample_analysis.multinomial_data.v_p_protein_data.size(); ++j) {
					if (temp_v_multinomial_frequency[j].multinomial_frequency > 0.1) {
						fout_multinomial_element_nomatch << temp_v_multinomial_frequency[j].p_protein_data->protein_name << " (";
						fout_multinomial_element_nomatch << temp_v_multinomial_frequency[j].multinomial_frequency << "), ";
					}
				}
				fout_multinomial_element_nomatch << "\n";
			}
		}
	}

	void fout_multinomial_contaminants_report(filesystem& par_filesystem, sample_analysis& par_sample_analysis, const string& par_filename_tag) {
		struct output_report {
		public:
			double peptide_theoretical_mz;
			double peptide_observed_mz;
			double peptide_z;
			double peptide_rt;
			string peptide_withmod;
			string peptide_filtered;
			string peptide_associated_protein_type;
			vector<multinomial_frequency_type> v_peptide_associated_proteins;
		};
		output_report temp_contaminants_report{};
		multimap<double, output_report> v_contaminants_report_map{};
		string output_multinomial_element = par_filesystem.directory + par_filesystem.filename + "_" + par_filename_tag + "_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_multinomial_element;
		fout_multinomial_element.open(output_multinomial_element);
		for (auto i = 0; i < par_sample_analysis.multinomial_data.v_p_peptide_analysis.size(); ++i) {
			vector<multinomial_frequency_type> temp_v_multinomial_frequency{};
			for (auto j = 0; j < par_sample_analysis.multinomial_data.v_p_protein_data.size(); ++j) {
				multinomial_frequency_type temp_multinomial_frequency{};
				temp_multinomial_frequency.p_protein_data = par_sample_analysis.multinomial_data.v_p_protein_data[j];
				temp_multinomial_frequency.multinomial_frequency = par_sample_analysis.multinomial_data.v2_frequency[i][j];
				temp_v_multinomial_frequency.push_back(temp_multinomial_frequency);
			}
			fpf_multinomial::sort_v_multinomial_frequency(temp_v_multinomial_frequency);
			if ((!temp_v_multinomial_frequency.empty()) && (temp_v_multinomial_frequency.begin()->multinomial_frequency > 0.2)) {
				//if ((temp_v_multinomial_frequency.begin()->p_protein_data->protein_type == "UNIPROT")
					//|| (temp_v_multinomial_frequency.begin()->p_protein_data->protein_type == "CONT")) {
				for (const auto& itr_peptide_data : par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->v_peptide_data) {
					temp_contaminants_report.peptide_theoretical_mz = ((itr_peptide_data->peptide_m + itr_peptide_data->peptide_z) / itr_peptide_data->peptide_z);
					temp_contaminants_report.peptide_observed_mz = itr_peptide_data->peptide_mz;
					temp_contaminants_report.peptide_z = itr_peptide_data->peptide_z;
					temp_contaminants_report.peptide_rt = itr_peptide_data->peptide_rt;
					temp_contaminants_report.peptide_withmod = itr_peptide_data->peptide_withmod;
					if (temp_v_multinomial_frequency.begin()->p_protein_data->protein_type == "UNIPROT") {
						temp_contaminants_report.peptide_associated_protein_type = "CONT";
					}
					else {
						temp_contaminants_report.peptide_associated_protein_type = temp_v_multinomial_frequency.begin()->p_protein_data->protein_type;
					}
					temp_contaminants_report.peptide_filtered = par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->peptide_filtered;
					temp_contaminants_report.v_peptide_associated_proteins = temp_v_multinomial_frequency;
					v_contaminants_report_map.insert(std::make_pair(temp_contaminants_report.peptide_theoretical_mz, temp_contaminants_report));
				}
				//}
			}
			else {
				for (const auto& itr_peptide_data : par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->v_peptide_data) {
					temp_contaminants_report.peptide_theoretical_mz = ((itr_peptide_data->peptide_m + itr_peptide_data->peptide_z) / itr_peptide_data->peptide_z);
					temp_contaminants_report.peptide_observed_mz = itr_peptide_data->peptide_mz;
					temp_contaminants_report.peptide_z = itr_peptide_data->peptide_z;
					temp_contaminants_report.peptide_rt = itr_peptide_data->peptide_rt;
					temp_contaminants_report.peptide_withmod = itr_peptide_data->peptide_withmod;
					temp_contaminants_report.peptide_associated_protein_type = "NA";
					temp_contaminants_report.peptide_filtered = par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->peptide_filtered;
					temp_contaminants_report.v_peptide_associated_proteins = temp_v_multinomial_frequency;
					v_contaminants_report_map.insert(std::make_pair(temp_contaminants_report.peptide_theoretical_mz, temp_contaminants_report));
				}
			}
		}
		fout_multinomial_element << "peptide_theoretical_mz,";
		fout_multinomial_element << "peptide_observed_mz,";
		fout_multinomial_element << "peptide_z,";
		fout_multinomial_element << "peptide_rt,";
		fout_multinomial_element << "peptide_withmod,";
		fout_multinomial_element << "peptide_filtered,";
		fout_multinomial_element << "associated_protein_type,";
		fout_multinomial_element << "associated_protein_name,";
		fout_multinomial_element << "\n";
		for (const auto& itr_contaminants_report_map : v_contaminants_report_map) {
			fout_multinomial_element << std::fixed << std::setprecision(3) << itr_contaminants_report_map.second.peptide_theoretical_mz << ","; 
			fout_multinomial_element << std::fixed << std::setprecision(3) << itr_contaminants_report_map.second.peptide_observed_mz << ",";
			fout_multinomial_element << std::fixed << std::setprecision(0) << itr_contaminants_report_map.second.peptide_z << ",";
			fout_multinomial_element << std::fixed << std::setprecision(2) << itr_contaminants_report_map.second.peptide_rt << ",";
			fout_multinomial_element << itr_contaminants_report_map.second.peptide_withmod << ",";
			fout_multinomial_element << itr_contaminants_report_map.second.peptide_filtered << ",";
			fout_multinomial_element << itr_contaminants_report_map.second.peptide_associated_protein_type << ",";
			for (const auto& itr_peptide_associated_proteins : itr_contaminants_report_map.second.v_peptide_associated_proteins) {
				if (itr_peptide_associated_proteins.multinomial_frequency > 0.2) {
					fout_multinomial_element << itr_peptide_associated_proteins.p_protein_data->protein_name << " (";
					fout_multinomial_element << itr_peptide_associated_proteins.multinomial_frequency << ")|";
				}
			}
			fout_multinomial_element << ",\n";
		}
	}

	multimap<double, double> fout_multinomial_contaminants_list(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		string output_multinomial_element = par_filesystem.directory + par_filesystem.filename + "_contaminants_list_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_multinomial_element;
		fout_multinomial_element.open(output_multinomial_element);
		multimap<double, double> v_contaminants_list_map{};
		for (auto i = 0; i < par_sample_analysis.multinomial_data.v_p_peptide_analysis.size(); ++i) {
			vector<multinomial_frequency_type> temp_v_multinomial_frequency{};
			for (auto j = 0; j < par_sample_analysis.multinomial_data.v_p_protein_data.size(); ++j) {
				multinomial_frequency_type temp_multinomial_frequency{};
				temp_multinomial_frequency.p_protein_data = par_sample_analysis.multinomial_data.v_p_protein_data[j];
				temp_multinomial_frequency.multinomial_frequency = par_sample_analysis.multinomial_data.v2_frequency[i][j];
				temp_v_multinomial_frequency.push_back(temp_multinomial_frequency);
			}
			fpf_multinomial::sort_v_multinomial_frequency(temp_v_multinomial_frequency);
			if (!temp_v_multinomial_frequency.empty()) {
				if ((temp_v_multinomial_frequency.begin()->p_protein_data->protein_type == "UNIPROT")
					|| (temp_v_multinomial_frequency.begin()->p_protein_data->protein_type == "CONT")) {
					if (temp_v_multinomial_frequency.begin()->multinomial_frequency > 0.2) {
						for (auto k = 0; k < par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->v_peptide_withoutmod_mz.size(); ++k) {
							double contaminants_list_map_first = par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->v_peptide_withoutmod_m[k] / par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->v_peptide_withoutmod_z[k];
							double contaminants_list_map_second = par_sample_analysis.multinomial_data.v_p_peptide_analysis[i]->v_peptide_withoutmod_rt[k];
							v_contaminants_list_map.insert(std::make_pair(contaminants_list_map_first, contaminants_list_map_second));
						}
					}
				}
			}
		}
		for (const auto& itr_v_contaminants_list_map : v_contaminants_list_map) {
			fout_multinomial_element << itr_v_contaminants_list_map.first << ",";
			fout_multinomial_element << itr_v_contaminants_list_map.second << ",";
			fout_multinomial_element << "\n";
		}
		return v_contaminants_list_map;
	}

	void temp_fout_multinomial_contaminants_list(filesystem& par_filesystem, sample_analysis& par_sample_analysis, map<double, double>& par_v_selected_combined_contaminants_list) {
		string output_multinomial_element = par_filesystem.directory + par_filesystem.filename + "_contaminants_list_combined_" + par_sample_analysis.peptide_assignment_method + ".txt";
		std::ofstream fout_multinomial_element;
		fout_multinomial_element.open(output_multinomial_element);
		for (const auto& itr_selected_combined_contaminants_list : par_v_selected_combined_contaminants_list) {
			fout_multinomial_element << itr_selected_combined_contaminants_list.first << ",";
			fout_multinomial_element << itr_selected_combined_contaminants_list.second << ",";
			fout_multinomial_element << "\n";
		}
	}

	void fout_multinomial_protein_score(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::string output_multinomial_comparison = par_filesystem.directory + par_filesystem.filename + "_protein_score_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_multinomial_protein_score;
		fout_multinomial_protein_score.open(output_multinomial_comparison);
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if (itr_v_protein_analysis.p_protein_data->protein_type == "IGV") {
				fout_multinomial_protein_score << itr_v_protein_analysis.p_protein_data->protein_name << ",";
				fout_multinomial_protein_score << (itr_v_protein_analysis.protein_score / IgFamily::REPORT_PROTEIN_SCORE_OUTPUT_SCALE) << "\n";
			}
		}
	}

	void fout_multinomial_protein_density(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::string output_multinomial_comparison = par_filesystem.directory + par_filesystem.filename + "_protein_density_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_multinomial_protein_density;
		fout_multinomial_protein_density.open(output_multinomial_comparison);
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if (itr_v_protein_analysis.p_protein_data->protein_type == "IGV") {
				fout_multinomial_protein_density << itr_v_protein_analysis.p_protein_data->protein_name << ",";
				fout_multinomial_protein_density << itr_v_protein_analysis.protein_density << "\n";
			}
		}
	}

	void fout_protein_pseudoabundance_score(const filesystem& par_filesystem, const sample_analysis& par_sample_analysis, const double& par_parameter_density_threshold, const string& par_filename_tag) {
		std::string output_protein_pseudoabundance_score_uniquepeptides = par_filesystem.directory + par_filesystem.filename + "_protein_pseudoabundance_score_" + par_filename_tag + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_protein_pseudoabundance_score_uniquepeptides;
		fout_protein_pseudoabundance_score_uniquepeptides.open(output_protein_pseudoabundance_score_uniquepeptides);
		std::string output_protein_pseudoabundance_peptides_uniquepeptides = par_filesystem.directory + par_filesystem.filename + "_protein_pseudoabundance_peptides_" + par_filename_tag + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_protein_pseudoabundance_peptides_uniquepeptides;
		fout_protein_pseudoabundance_peptides_uniquepeptides.open(output_protein_pseudoabundance_peptides_uniquepeptides);
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			double score_model0{};
			double score_model1{};
			double score_model2{};
			for (const auto& itr_v_homology_data_combined_by_protein : itr_v_protein_analysis.v_homology_data_combined_by_protein) {
				if (itr_v_homology_data_combined_by_protein.blastp_homology_density > par_parameter_density_threshold) {
					for (auto i = 0; i < itr_v_homology_data_combined_by_protein.p_peptide_analysis->v_peptide_data.size(); ++i) {
						score_model0 += double{ 1 };
						score_model1 += std::pow(itr_v_homology_data_combined_by_protein.blastp_homology_density_conjugated, IgFamily::PARAMETER_SCORE_CONJUGATION_WEIGHT);
						score_model2 += std::pow(itr_v_homology_data_combined_by_protein.blastp_homology_density_conjugated, double(2.5));
					}
				}
			}
			if ((itr_v_protein_analysis.p_protein_data->protein_type == "IGV")
				&& (score_model1 != double(0))) {
				fout_protein_pseudoabundance_score_uniquepeptides << itr_v_protein_analysis.p_protein_data->protein_name << ",";
				fout_protein_pseudoabundance_score_uniquepeptides << score_model0 << ",";
				fout_protein_pseudoabundance_score_uniquepeptides << score_model1 << ",";
				fout_protein_pseudoabundance_score_uniquepeptides << score_model2 << ",";
				if (&itr_v_protein_analysis != &par_sample_analysis.v_protein_analysis.back()) {
					fout_protein_pseudoabundance_score_uniquepeptides << "\n";
				}
			}
			for (const auto& itr_v_homology_data_combined_by_protein : itr_v_protein_analysis.v_homology_data_combined_by_protein) {
				if (itr_v_homology_data_combined_by_protein.blastp_homology_density > par_parameter_density_threshold) {
					for (auto i = 0; i < itr_v_homology_data_combined_by_protein.p_peptide_analysis->v_peptide_data.size(); ++i) {
						if (itr_v_protein_analysis.p_protein_data->protein_type == "IGV") {
							fout_protein_pseudoabundance_peptides_uniquepeptides << itr_v_homology_data_combined_by_protein.blastp_query << ",";
							fout_protein_pseudoabundance_peptides_uniquepeptides << itr_v_homology_data_combined_by_protein.blastp_subject_accession << ",";;
							if (&itr_v_protein_analysis != &par_sample_analysis.v_protein_analysis.back()) {
								fout_protein_pseudoabundance_peptides_uniquepeptides << "\n";
							}
						}
					}
				}
			}
		}
	}

	void fout_html_report(filesystem& par_filesystem, sample_analysis& par_sample_analysis, bool par_alloutput, bool par_summary, const double& par_protein_score_threshold, const string& par_filename_tag) {
		std::string output_html_report{};
		if (par_alloutput && !par_summary) {
			output_html_report = par_filesystem.directory + par_filesystem.filename + "_" + par_filename_tag + "_" + par_sample_analysis.peptide_assignment_method + ".html";
		}
		if (!par_alloutput && !par_summary) {
			output_html_report = par_filesystem.directory + par_filesystem.filename + "_" + par_filename_tag + "_" + par_sample_analysis.peptide_assignment_method + ".html";
		}
		if (par_alloutput && par_summary) {
			output_html_report = par_filesystem.directory + par_filesystem.filename + "_" + par_filename_tag + "_" + par_sample_analysis.peptide_assignment_method + ".html";
		}
		if (!par_alloutput && par_summary) {
			output_html_report = par_filesystem.directory + par_filesystem.filename + "_" + par_filename_tag + "_" + par_sample_analysis.peptide_assignment_method + ".html";
		}
		std::ofstream fout_html_report;
		fout_html_report.open(output_html_report);
		vector<string> dummy;
		fout_html_report << "\
							<!DOCTYPE html>\n\
							<head>\n\
							<meta charset = \"UTF-8\">\n\
							<title>HTML Template</title>\n\
							</head>\n\
							<body>\n\
							<p><font face=\"Lucida Console\" size=\"3\" color=\"black\">";
		fout_html_report << "<style> \
						.menu{ \
						white-space: nowrap; \
						} \
						.mismatch { \
						color: black; \
						border-bottom: 2px solid black; \
						} \
						</style>";
		fout_html_report << "\n<br>" << par_filesystem.filename;
		fout_html_report << "\n\n<br><br>" << par_filesystem.fileversion;
		if (IgFamily::FILESYSTEM_MODE) {
			fout_html_report << "&nbsp&nbsp&nbsp" << par_filesystem.enzyme;
			fout_html_report << "&nbsp&nbsp&nbsp" << par_filesystem.denono_deltamass << "&nbspDa";
		}
		double protein_analysis_score_max{};
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if ((itr_v_protein_analysis.p_protein_data->protein_type == "IGV") && (itr_v_protein_analysis.protein_score > protein_analysis_score_max)) {
				protein_analysis_score_max = itr_v_protein_analysis.protein_score;
			}
		}
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			if ((par_alloutput) || (itr_v_protein_analysis.p_protein_data->protein_type == "IGV")) {
				if (((itr_v_protein_analysis.protein_score / IgFamily::REPORT_PROTEIN_SCORE_OUTPUT_SCALE) >= par_protein_score_threshold)
				&& (itr_v_protein_analysis.protein_density >= IgFamily::REPORT_PROTEIN_DENSITY_THRESHOLD)) {
					fout_html_report << "\n\n\n<br><br><br> " << itr_v_protein_analysis.p_protein_data->protein_name;
					fout_html_report << "&nbsp&nbsp&nbspScore: " << std::fixed << std::setprecision(2) << (itr_v_protein_analysis.protein_score / IgFamily::REPORT_PROTEIN_SCORE_OUTPUT_SCALE);
					fout_html_report << "&nbsp&nbsp&nbspDensity: " << std::fixed << std::setprecision(3) << itr_v_protein_analysis.protein_density;
					fout_html_report << "&nbsp&nbsp&nbspEffective SC: " << std::fixed << std::setprecision(2) << itr_v_protein_analysis.protein_effective_spectral_count;
					fout_html_report << "&nbsp&nbsp&nbspCoverage: " << std::fixed << std::setprecision(0) << itr_v_protein_analysis.proteinconstruct_sequencecoverage << "%";
					fout_html_report << "\n\n<br><br> " << itr_v_protein_analysis.p_protein_data->protein_protein;
					fout_html_report << "\n\n<br> ";
					size_t i{};
					for (const auto& itr_proteinconstruct : itr_v_protein_analysis.proteinconstruct) {
						if (itr_proteinconstruct.aminoacid != '.') {
							if ((itr_proteinconstruct.aminoacid != itr_v_protein_analysis.p_protein_data->protein_protein.at(i)) && !((itr_proteinconstruct.aminoacid == 'L') && (itr_v_protein_analysis.p_protein_data->protein_protein.at(i) == 'I'))) {
								fout_html_report << "<span class=\"mismatch\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density_conjugated) >= 0.80) {
								fout_html_report << "<font color=\"#4c62d6\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density_conjugated < 0.80) && (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density_conjugated >= 0.50)) {
								fout_html_report << "<font color=\"#239B56\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density_conjugated < 0.50) && (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density_conjugated >= 0.20)) {
								fout_html_report << "<font color=\"#E67E22\">";
							}
							if (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density_conjugated < 0.20) {
								fout_html_report << "<font color=\"red\">";
							}
						}
						fout_html_report << itr_proteinconstruct.aminoacid;
						if (itr_proteinconstruct.aminoacid != itr_v_protein_analysis.p_protein_data->protein_protein.at(i)) {
							fout_html_report << "</span>";
						}
						if (itr_proteinconstruct.aminoacid != '.') {
							fout_html_report << "</font>";
						}
						++i;
					}
					fout_html_report << "&nbsp&nbsp&nbsp" << "Conjugated density";
					fout_html_report << "\n\n<br> ";
					i = size_t();
					for (const auto& itr_proteinconstruct : itr_v_protein_analysis.proteinconstruct) {
						if (itr_proteinconstruct.aminoacid != '.') {
							if ((itr_proteinconstruct.aminoacid != itr_v_protein_analysis.p_protein_data->protein_protein.at(i)) && !((itr_proteinconstruct.aminoacid == 'L') && (itr_v_protein_analysis.p_protein_data->protein_protein.at(i) == 'I'))) {
								fout_html_report << "<span class=\"mismatch\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density) >= 0.80) {
								fout_html_report << "<font color=\"#4c62d6\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density < 0.80) && (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density >= 0.50)) {
								fout_html_report << "<font color=\"#239B56\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density < 0.50) && (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density >= 0.20)) {
								fout_html_report << "<font color=\"#E67E22\">";
							}
							if (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_density < 0.20) {
								fout_html_report << "<font color=\"red\">";
							}
						}
						fout_html_report << itr_proteinconstruct.aminoacid;
						if (itr_proteinconstruct.aminoacid != itr_v_protein_analysis.p_protein_data->protein_protein.at(i)) {
							fout_html_report << "</span>";
						}
						if (itr_proteinconstruct.aminoacid != '.') {
							fout_html_report << "</font>";
						}
						++i;
					}
					fout_html_report << "&nbsp&nbsp&nbsp" << "Density";
					fout_html_report << "\n\n<br> ";
					i = size_t();
					for (const auto& itr_proteinconstruct : itr_v_protein_analysis.proteinconstruct) {
						if (itr_proteinconstruct.aminoacid != '.') {
							if ((itr_proteinconstruct.aminoacid != itr_v_protein_analysis.p_protein_data->protein_protein.at(i)) && !((itr_proteinconstruct.aminoacid == 'L') && (itr_v_protein_analysis.p_protein_data->protein_protein.at(i) == 'I'))) {
								fout_html_report << "<span class=\"mismatch\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed_conjugated) >= std::pow(45, IgFamily::PARAMETER_HOMOLOGY_WEIGHT)) {
								fout_html_report << "<font color=\"#4c62d6\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed_conjugated < std::pow(45, IgFamily::PARAMETER_HOMOLOGY_WEIGHT)) && (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed_conjugated >= std::pow(35, IgFamily::PARAMETER_HOMOLOGY_WEIGHT))) {
								fout_html_report << "<font color=\"#239B56\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed_conjugated < std::pow(35, IgFamily::PARAMETER_HOMOLOGY_WEIGHT)) && (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed_conjugated >= std::pow(25, IgFamily::PARAMETER_HOMOLOGY_WEIGHT))) {
								fout_html_report << "<font color=\"#E67E22\">";
							}
							if (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed_conjugated < std::pow(25, IgFamily::PARAMETER_HOMOLOGY_WEIGHT)) {
								fout_html_report << "<font color=\"red\">";
							}
						}
						fout_html_report << itr_proteinconstruct.aminoacid;
						if (itr_proteinconstruct.aminoacid != itr_v_protein_analysis.p_protein_data->protein_protein.at(i)) {
							fout_html_report << "</span>";
						}
						if (itr_proteinconstruct.aminoacid != '.') {
							fout_html_report << "</font>";
						}
						++i;
					}
					fout_html_report << "&nbsp&nbsp&nbsp" << "Conjugated homology";
					fout_html_report << "\n\n<br> ";
					i = size_t();
					for (const auto& itr_proteinconstruct : itr_v_protein_analysis.proteinconstruct) {
						if (itr_proteinconstruct.aminoacid != '.') {
							if ((itr_proteinconstruct.aminoacid != itr_v_protein_analysis.p_protein_data->protein_protein.at(i)) && !((itr_proteinconstruct.aminoacid == 'L') && (itr_v_protein_analysis.p_protein_data->protein_protein.at(i) == 'I'))) {
								fout_html_report << "<span class=\"mismatch\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed) >= std::pow(45, IgFamily::PARAMETER_HOMOLOGY_WEIGHT)) {
								fout_html_report << "<font color=\"#4c62d6\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed < std::pow(45, IgFamily::PARAMETER_HOMOLOGY_WEIGHT)) && (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed >= std::pow(35, IgFamily::PARAMETER_HOMOLOGY_WEIGHT))) {
								fout_html_report << "<font color=\"#239B56\">";
							}
							if ((itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed < std::pow(35, IgFamily::PARAMETER_HOMOLOGY_WEIGHT)) && (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed >= std::pow(25, IgFamily::PARAMETER_HOMOLOGY_WEIGHT))) {
								fout_html_report << "<font color=\"#E67E22\">";
							}
							if (itr_proteinconstruct.proteinconstruct_homology_data.blastp_homology_transformed < std::pow(25, IgFamily::PARAMETER_HOMOLOGY_WEIGHT)) {
								fout_html_report << "<font color=\"red\">";
							}
						}
						fout_html_report << itr_proteinconstruct.aminoacid;
						if (itr_proteinconstruct.aminoacid != itr_v_protein_analysis.p_protein_data->protein_protein.at(i)) {
							fout_html_report << "</span>";
						}
						if (itr_proteinconstruct.aminoacid != '.') {
							fout_html_report << "</font>";
						}
						++i;
					}
					fout_html_report << "&nbsp&nbsp&nbsp" << "Homology";
					fout_html_report << "\n\n<br> ";
					i = size_t();
					for (const auto& itr_proteinconstruct : itr_v_protein_analysis.proteinconstruct) {
						if (itr_proteinconstruct.aminoacid_localconfidence > 0) {
							if (itr_proteinconstruct.aminoacid != '.') {
								if ((itr_proteinconstruct.aminoacid != itr_v_protein_analysis.p_protein_data->protein_protein.at(i)) && !((itr_proteinconstruct.aminoacid == 'L') && (itr_v_protein_analysis.p_protein_data->protein_protein.at(i) == 'I'))) {
									fout_html_report << "<span class=\"mismatch\">";
								}
								if (itr_proteinconstruct.aminoacid_localconfidence >= 90) {
									fout_html_report << "<font color=\"#4c62d6\">";
								}
								if ((itr_proteinconstruct.aminoacid_localconfidence < 90) && (itr_proteinconstruct.aminoacid_localconfidence >= 80)) {
									fout_html_report << "<font color=\"#239B56\">";
								}
								if ((itr_proteinconstruct.aminoacid_localconfidence < 80) && (itr_proteinconstruct.aminoacid_localconfidence >= 60)) {
									fout_html_report << "<font color=\"#E67E22\">";
								}
								if (itr_proteinconstruct.aminoacid_localconfidence < 60) {
									fout_html_report << "<font color=\"red\">";
								}
							}
							fout_html_report << itr_proteinconstruct.aminoacid;
							if (itr_proteinconstruct.aminoacid != itr_v_protein_analysis.p_protein_data->protein_protein.at(i)) {
								fout_html_report << "</span>";
							}
							if (itr_proteinconstruct.aminoacid != '.') {
								fout_html_report << "</font>";
							}
						}
						else {
							fout_html_report << '.';
						}
						++i;
					}
					fout_html_report << "&nbsp&nbsp&nbsp" << "Local&nbspconfidence";
					fout_html_report << "\n\n<br><br><br>";
					if (!par_summary) {
						for (const auto& itr_proteinconstruct : itr_v_protein_analysis.proteinconstruct) {
							fout_html_report << "&nbsp";
						}						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
						//fout_html_report << "Score ";
						//fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
						fout_html_report << "Sc dens";
						fout_html_report << "&nbsp&nbsp&nbsp";
						fout_html_report << "Spec";
						fout_html_report << "&nbsp&nbsp";
						fout_html_report << "Homo";
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
						fout_html_report << "Conj";
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp";
						fout_html_report << "Dens";
						fout_html_report << "<br>";
						for (const auto& itr_homology_data : itr_v_protein_analysis.v_homology_data_combined_by_protein) {
							if ((itr_homology_data.blastp_homology_density_conjugated >= IgFamily::REPORT_QUERY_PARAMETER_DENSITY_CONJUGATED_THRESHOLD)
								&& (itr_homology_data.score_density >= IgFamily::REPORT_QUERY_PARAMETER_SCORE_DENSITY_THRESHOLD)) {
								fout_html_report << "\n<br> ";
								size_t st_mismatch{};
								for (auto i = 0; i < itr_homology_data.alignment.length(); i) {
									if (itr_homology_data.alignment.at(i) == '.') {
										if (i < itr_homology_data.p_protein_data->protein_protein.length()) {
											fout_html_report << ".";
										}
										++i;
									}
									else {
										for (const auto& itr_denovo_aminoacid : itr_homology_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->v_denovo_aminoacid) {
											if (((i >= itr_homology_data.alignment.length()) || (itr_denovo_aminoacid.aminoacid != itr_homology_data.p_protein_data->protein_protein.at(i)) && !((itr_denovo_aminoacid.aminoacid == 'L') && (itr_v_protein_analysis.p_protein_data->protein_protein.at(i) == 'I')))) {
												fout_html_report << "<span class=\"mismatch\">";
											}
											if (itr_denovo_aminoacid.aminoacid_localconfidence >= 90) {
												fout_html_report << "<font color=\"#4c62d6\">" << itr_denovo_aminoacid.aminoacid << "</font>";
											}
											if ((itr_denovo_aminoacid.aminoacid_localconfidence < 90) && (itr_denovo_aminoacid.aminoacid_localconfidence >= 80)) {
												fout_html_report << "<font color=\"#239B56\">" << itr_denovo_aminoacid.aminoacid << "</font>";
											}
											if ((itr_denovo_aminoacid.aminoacid_localconfidence < 80) && (itr_denovo_aminoacid.aminoacid_localconfidence >= 60)) {
												fout_html_report << "<font color=\"#E67E22\">" << itr_denovo_aminoacid.aminoacid << "</font>";
											}
											if (itr_denovo_aminoacid.aminoacid_localconfidence < 60) {
												fout_html_report << "<font color=\"red\">" << itr_denovo_aminoacid.aminoacid << "</font>";
											}
											if ((i >= itr_homology_data.alignment.length()) || (itr_denovo_aminoacid.aminoacid != itr_homology_data.p_protein_data->protein_protein.at(i))) {
												fout_html_report << "</span>";
											}
											++i;
										}
									}
									st_mismatch = i;
								}
								st_mismatch = (st_mismatch - itr_homology_data.alignment.length());
								for (auto j = 0; j < (5 - st_mismatch); ++j) {
									fout_html_report << "&nbsp";
								}
								fout_html_report << std::fixed << std::setprecision(3) << itr_homology_data.score_density;
								fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
								//fout_html_report << std::fixed << std::setprecision(0) << (itr_homology_data.score / std::pow(itr_homology_data.score, (IgFamily::PARAMETER_HOMOLOGY_WEIGHT - 1)));
								//fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
								fout_html_report << itr_homology_data.denovo_replicate_count;
								if (itr_homology_data.denovo_replicate_count < 1) {
									fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
								}
								else {
									for (auto j = 0; j < (5 - std::floor(std::log10(itr_homology_data.denovo_replicate_count))); ++j) {
										fout_html_report << "&nbsp";
									}
								}
								fout_html_report << std::fixed << std::setprecision(2) << itr_homology_data.blastp_homology;
								if (itr_homology_data.blastp_homology_transformed < 1) {
									fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
								}
								else {
									for (auto j = 0; j < (5 - std::floor(std::log10(itr_homology_data.blastp_homology))); ++j) {
										fout_html_report << "&nbsp";
									}
								}
								if ((itr_homology_data.blastp_homology_density_conjugated) >= 0.80) {
									fout_html_report << "<font color=\"#4c62d6\">";
								}
								if ((itr_homology_data.blastp_homology_density_conjugated < 0.80) && (itr_homology_data.blastp_homology_density_conjugated >= 0.50)) {
									fout_html_report << "<font color=\"#239B56\">";
								}
								if ((itr_homology_data.blastp_homology_density_conjugated < 0.50) && (itr_homology_data.blastp_homology_density_conjugated >= 0.20)) {
									fout_html_report << "<font color=\"#E67E22\">";
								}
								if (itr_homology_data.blastp_homology_density_conjugated < 0.20) {
									fout_html_report << "<font color=\"red\">";
								}
								fout_html_report << std::fixed << std::setprecision(3) << itr_homology_data.blastp_homology_density_conjugated;
								fout_html_report << "</font>";
								fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
								if ((itr_homology_data.blastp_homology_density) >= 0.80) {
									fout_html_report << "<font color=\"#4c62d6\">";
								}
								if ((itr_homology_data.blastp_homology_density < 0.80) && (itr_homology_data.blastp_homology_density >= 0.50)) {
									fout_html_report << "<font color=\"#239B56\">";
								}
								if ((itr_homology_data.blastp_homology_density < 0.50) && (itr_homology_data.blastp_homology_density >= 0.20)) {
									fout_html_report << "<font color=\"#E67E22\">";
								}
								if (itr_homology_data.blastp_homology_density < 0.20) {
									fout_html_report << "<font color=\"red\">";
								}
								fout_html_report << std::fixed << std::setprecision(3) << itr_homology_data.blastp_homology_density;
								fout_html_report << "</font>";
								fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
								for (const auto& itr_homology_data_aggregated_by_homology_distribution : itr_homology_data.v_homology_data_aggregated_by_homology_distribution) {
									if (itr_homology_data_aggregated_by_homology_distribution->blastp_homology_density_conjugated >= IgFamily::REPORT_V_HOMOLOGY_DATA_AGGREGATED_BY_PROTEIN_COJUGATED_DENSITY_THRESHOLD) {
										fout_html_report << "&nbsp&nbsp";
										fout_html_report << itr_homology_data_aggregated_by_homology_distribution->blastp_subject_accession;
										fout_html_report << "(";
										fout_html_report << std::setprecision(2) << itr_homology_data_aggregated_by_homology_distribution->blastp_homology_density_conjugated;
										fout_html_report << "),";
									}
								}
							}
						}
					}
				}
			}
		}
		fout_html_report << "\
		</p></font>\n \
	</body>\n \
</html>\n ";
	}

	void fout_v_sample_analysis_comparison(vector<filesystem>& par_v_filesystem, filesystem& par_filesystem) {
		string sample_analysis_comparison{};
		for (const auto& itr_v_filesystem : par_v_filesystem) {

		}
		string output_v_sample_analysis_comparison = par_filesystem.directory + par_filesystem.filename + "_sample_comparison_" + par_filesystem.v_sample_analysis.begin()->peptide_assignment_method + ".txt";
		std::ofstream fout_v_sample_analysis_comparison;
		fout_v_sample_analysis_comparison.open(output_v_sample_analysis_comparison);
		std::cout << "\n ...outputting sample comparison for " << par_filesystem.filename;
		fout_v_sample_analysis_comparison << sample_analysis_comparison;
	}

}
#endif

