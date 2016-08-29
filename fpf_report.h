// * * fpf_report.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_REPORT
#define	FPF_REPORT

#include <cstdlib>
#include <iomanip>
#include <math.h>
#include <string>
#include <vector>

#include "IgFamily.h"
#include "fpf_filesystem.h"


namespace fpf_report {

	using std::string;
	using std::vector;

	typedef fpf_filesystem::filesystem filesystem;
	typedef fpf_filesystem::sample_analysis sample_analysis;

	void fout_v_peptide_data(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::string output_v_peptide_data = par_filesystem.directory + par_filesystem.filename + "_peptide_summary_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_v_peptide_data;
		fout_v_peptide_data.open(output_v_peptide_data);
		fout_v_peptide_data << "key,scan_ID,peptide_filtered,peptide_withmod,denovo_peptide,\n";
		for (const auto& itr_v_peptide_data : par_sample_analysis.v_peptide_data) {
			fout_v_peptide_data << itr_v_peptide_data.key_peptide_data << ",";
			fout_v_peptide_data << itr_v_peptide_data.scan_ID << ",";
			fout_v_peptide_data << itr_v_peptide_data.peptide_filtered << ",";
			fout_v_peptide_data << itr_v_peptide_data.peptide_withmod << ",";
			for (const auto& itr_v_denovo_aminoacid : itr_v_peptide_data.denovo_peptide_data.v_denovo_aminoacid) {
				fout_v_peptide_data << itr_v_denovo_aminoacid.aminoacid;
				fout_v_peptide_data << "[" << itr_v_denovo_aminoacid.aminoacid_localconfidence << "]";
			}
			fout_v_peptide_data << ",\n";
		}
	}

	void fout_v_peptide_analysis(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::string output_v_peptide_analysis = par_filesystem.directory + par_filesystem.filename + "_peptide_analysis_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_v_peptide_analysis;
		fout_v_peptide_analysis.open(output_v_peptide_analysis);
		fout_v_peptide_analysis << "key,peptide_filtered,replicate_count,\n";
		for (const auto& itr_v_peptide_anaylsis : par_sample_analysis.v_peptide_analysis) {
			fout_v_peptide_analysis << itr_v_peptide_anaylsis.key_peptide_analysis << ",";
			fout_v_peptide_analysis << itr_v_peptide_anaylsis.peptide_filtered << ",";
			fout_v_peptide_analysis << itr_v_peptide_anaylsis.replicate_count << ",";
			fout_v_peptide_analysis << "\n";
		}
	}

	void fout_v_protein_analysis(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::string output_v_protein_analysis = par_filesystem.directory + par_filesystem.filename + "_protein_analysis_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_v_protein_analysis;
		fout_v_protein_analysis.open(output_v_protein_analysis);
		fout_v_protein_analysis << "key,protein_name,protein_protein,protein_score,proteinconstruct_sequencecoverage,\n";
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			fout_v_protein_analysis << itr_v_protein_analysis.key_protein_analysis << ",";
			fout_v_protein_analysis << itr_v_protein_analysis.p_protein_data->protein_name << ",";
			fout_v_protein_analysis << itr_v_protein_analysis.p_protein_data->protein_protein << ",";
			fout_v_protein_analysis << itr_v_protein_analysis.protein_score << ",";
			fout_v_protein_analysis << itr_v_protein_analysis.proteinconstruct_sequencecoverage << ",";
			for (const auto& itr_v_proteinconstruct_from_denovos : itr_v_protein_analysis.proteinconstruct_from_denovo) {
				fout_v_protein_analysis << itr_v_proteinconstruct_from_denovos.aminoacid;
				fout_v_protein_analysis << "[" << itr_v_proteinconstruct_from_denovos.aminoacid_localconfidence << "]";
			}
			fout_v_protein_analysis << ",\n";
		}
	}

	void fout_multinomial_comparison(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::string output_multinomial_comparison = par_filesystem.directory + par_filesystem.filename + "_protein_score_" + par_sample_analysis.peptide_assignment_method + ".csv";
		std::ofstream fout_multinomial_comparison;
		fout_multinomial_comparison.open(output_multinomial_comparison);
		for (const auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis_selected_by_polymorphism) {
			if (itr_v_protein_analysis.p_protein_data->protein_type == "IG") {
				fout_multinomial_comparison << itr_v_protein_analysis.p_protein_data->protein_name << ",";
				fout_multinomial_comparison << itr_v_protein_analysis.protein_score << "\n";
			}
		}
	}

	void fout_html_report(filesystem& par_filesystem, sample_analysis& par_sample_analysis, bool par_alloutput, bool par_summary) {
		std::string output_html_report{};
		if (par_alloutput && !par_summary) {
			output_html_report = par_filesystem.directory + par_filesystem.filename + "_report_" + par_sample_analysis.peptide_assignment_method + ".html";
		}
		if (!par_alloutput && !par_summary) {
			output_html_report = par_filesystem.directory + par_filesystem.filename + "_report_IG_" + par_sample_analysis.peptide_assignment_method + ".html";
		}
		if (par_alloutput && par_summary) {
			output_html_report = par_filesystem.directory + par_filesystem.filename + "_report_summary_" + par_sample_analysis.peptide_assignment_method + ".html";
		}
		if (!par_alloutput && par_summary) {
			output_html_report = par_filesystem.directory + par_filesystem.filename + "_report_IG_summary_" + par_sample_analysis.peptide_assignment_method + ".html";
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
		for (const auto& itr_protein_analysis : par_sample_analysis.v_protein_analysis_selected_by_polymorphism) {
			if ((par_alloutput) || (itr_protein_analysis.p_protein_data->protein_type == "IG")) {
				if (itr_protein_analysis.protein_score > REPORT_SCORE_THRESHOLD) {
					fout_html_report << "\n\n\n<br><br><br> " << itr_protein_analysis.p_protein_data->protein_name;
					fout_html_report << "&nbsp&nbsp&nbspScore: " << std::fixed << std::setprecision(2) << itr_protein_analysis.protein_score;
					fout_html_report << "&nbsp&nbsp&nbspCoverage: " << std::fixed << std::setprecision(0) << itr_protein_analysis.proteinconstruct_sequencecoverage << "%";
					fout_html_report << "\n\n<br><br> " << itr_protein_analysis.p_protein_data->protein_protein;
					fout_html_report << "\n\n<br> ";
					size_t i{};
					for (const auto& itr_proteinconstruct_from_denovo : itr_protein_analysis.proteinconstruct_from_denovo) {
						if (itr_proteinconstruct_from_denovo.aminoacid_localconfidence > 0) {
							if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
								if ((itr_proteinconstruct_from_denovo.aminoacid != itr_protein_analysis.p_protein_data->protein_protein.at(i)) && !((itr_proteinconstruct_from_denovo.aminoacid == 'L') && (itr_protein_analysis.p_protein_data->protein_protein.at(i) == 'I'))) {
									fout_html_report << "<span class=\"mismatch\">";
								}
								if (itr_proteinconstruct_from_denovo.aminoacid_localconfidence > 90) {
									fout_html_report << "<font color=\"#4c62d6\">";
								}
								if ((itr_proteinconstruct_from_denovo.aminoacid_localconfidence <= 90) && (itr_proteinconstruct_from_denovo.aminoacid_localconfidence > 80)) {
									fout_html_report << "<font color=\"#239B56\">";
								}
								if ((itr_proteinconstruct_from_denovo.aminoacid_localconfidence <= 80) && (itr_proteinconstruct_from_denovo.aminoacid_localconfidence > 60)) {
									fout_html_report << "<font color=\"#E67E22\">";
								}
								if (itr_proteinconstruct_from_denovo.aminoacid_localconfidence <= 60) {
									fout_html_report << "<font color=\"red\">";
								}
							}
							fout_html_report << itr_proteinconstruct_from_denovo.aminoacid;
							if (itr_proteinconstruct_from_denovo.aminoacid != itr_protein_analysis.p_protein_data->protein_protein.at(i)) {
								fout_html_report << "</span>";
							}
							if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
								fout_html_report << "</font>";
							}
						}
						else {
							fout_html_report << '.';
						}
						++i;
					}
					fout_html_report << "&nbsp&nbsp&nbsp" << "Local&nbspconfidence";
					fout_html_report << "\n\n<br> ";
					i = size_t();
					for (const auto& itr_proteinconstruct_from_denovo : itr_protein_analysis.proteinconstruct_from_denovo) {
						if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
							if ((itr_proteinconstruct_from_denovo.aminoacid != itr_protein_analysis.p_protein_data->protein_protein.at(i)) && !((itr_proteinconstruct_from_denovo.aminoacid == 'L') && (itr_protein_analysis.p_protein_data->protein_protein.at(i) == 'I'))) {
								fout_html_report << "<span class=\"mismatch\">";
							}
							if (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed > 50) {
								fout_html_report << "<font color=\"#4c62d6\">";
							}
							if ((itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed <= 50) && (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed > 35)) {
								fout_html_report << "<font color=\"#239B56\">";
							}
							if ((itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed <= 35) && (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed > 20)) {
								fout_html_report << "<font color=\"#E67E22\">";
							}
							if (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed <= 20) {
								fout_html_report << "<font color=\"red\">";
							}
						}
						fout_html_report << itr_proteinconstruct_from_denovo.aminoacid;
						if (itr_proteinconstruct_from_denovo.aminoacid != itr_protein_analysis.p_protein_data->protein_protein.at(i)) {
							fout_html_report << "</span>";
						}
						if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
							fout_html_report << "</font>";
						}
						++i;
					}
					fout_html_report << "&nbsp&nbsp&nbsp" << "Homology";
					fout_html_report << "\n\n<br><br><br>";
					if (!par_summary) {
						for (const auto& itr_proteinconstruct_from_denovo : itr_protein_analysis.proteinconstruct_from_denovo) {
							fout_html_report << "&nbsp";
						}
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
						fout_html_report << "Total";
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp";
						fout_html_report << "Spec";
						fout_html_report << "&nbsp&nbsp";
						fout_html_report << "Score";
						fout_html_report << "&nbsp&nbsp&nbsp";
						fout_html_report << "Homo";
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
						fout_html_report << "Dist";
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp";
						fout_html_report << "DN&nbspconf";
						fout_html_report << "&nbsp&nbsp";
						fout_html_report << "DN&nbspconf&nbspavr";
						fout_html_report << "<br>";
						for (const auto& itr_homology_data : itr_protein_analysis.v_homology_data_combined_by_protein) {
							if (itr_homology_data.blastp_evalue_transformed > BLASTP_EVALUETRANSFORMED_THRESHOLD) {
								if ((itr_homology_data.blastp_parameter_score > REPORT_QUERY_ALIGNMENT_PARSCORE_OUTPUT_THRESHOLD)
									&& ((itr_homology_data.blastp_parameter_score * itr_homology_data.denovo_replicate_count) > REPORT_QUERY_ALIGNMENT_TOTALSCORE_OUTPUT_THRESHOLD)) {
									fout_html_report << "\n<br> ";
									int st_mismatch = int();
									for (auto i = 0; i < itr_homology_data.query_alignment.length(); i) {
										if (itr_homology_data.query_alignment.at(i) == '.') {
											if (i < itr_homology_data.p_protein_data->protein_protein.length()) {
												fout_html_report << ".";
											}
											++i;
										}
										else {
											for (const auto& itr_denovo_aminoacid : itr_homology_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->v_denovo_aminoacid) {
												if (((i >= itr_homology_data.query_alignment.length()) || (itr_denovo_aminoacid.aminoacid != itr_homology_data.p_protein_data->protein_protein.at(i)) && !((itr_denovo_aminoacid.aminoacid == 'L') && (itr_protein_analysis.p_protein_data->protein_protein.at(i) == 'I')))) {
													fout_html_report << "<span class=\"mismatch\">";
												}
												if (itr_denovo_aminoacid.aminoacid_localconfidence > 90) {
													fout_html_report << "<font color=\"#4c62d6\">" << itr_denovo_aminoacid.aminoacid << "</font>";
												}
												if ((itr_denovo_aminoacid.aminoacid_localconfidence <= 90) && (itr_denovo_aminoacid.aminoacid_localconfidence > 80)) {
													fout_html_report << "<font color=\"#239B56\">" << itr_denovo_aminoacid.aminoacid << "</font>";
												}
												if ((itr_denovo_aminoacid.aminoacid_localconfidence <= 80) && (itr_denovo_aminoacid.aminoacid_localconfidence > 60)) {
													fout_html_report << "<font color=\"#E67E22\">" << itr_denovo_aminoacid.aminoacid << "</font>";
												}
												if (itr_denovo_aminoacid.aminoacid_localconfidence <= 60) {
													fout_html_report << "<font color=\"red\">" << itr_denovo_aminoacid.aminoacid << "</font>";
												}
												if ((i >= itr_homology_data.query_alignment.length()) || (itr_denovo_aminoacid.aminoacid != itr_homology_data.p_protein_data->protein_protein.at(i))) {
													fout_html_report << "</span>";
												}
												++i;
											}
										}
										st_mismatch = i;
									}
									st_mismatch = (st_mismatch - itr_homology_data.query_alignment.length());
									for (auto j = 0; j < (5 - st_mismatch); ++j) {
										fout_html_report << "&nbsp";
									}
									fout_html_report << std::fixed << std::setprecision(2) << (itr_homology_data.blastp_parameter_score * itr_homology_data.denovo_replicate_count);
									if ((itr_homology_data.blastp_parameter_score * itr_homology_data.denovo_replicate_count) < 1) {
										fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
									}
									else {
										for (auto j = 0; j < (5 - std::floor(std::log10(itr_homology_data.blastp_parameter_score * itr_homology_data.denovo_replicate_count))); ++j) {
											fout_html_report << "&nbsp";
										}
									}
									fout_html_report << itr_homology_data.denovo_replicate_count;
									if (itr_homology_data.denovo_replicate_count < 1) {
										fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
									}
									else {
										for (auto j = 0; j < (5 - std::floor(std::log10(itr_homology_data.denovo_replicate_count))); ++j) {
											fout_html_report << "&nbsp";
										}
									}
									fout_html_report << std::fixed << std::setprecision(2) << itr_homology_data.blastp_parameter_score;
									if (itr_homology_data.blastp_parameter_score < 1) {
										fout_html_report << "&nbsp&nbsp&nbsp&nbsp";
									}
									else {
										for (auto j = 0; j < (4 - std::floor(std::log10(itr_homology_data.blastp_parameter_score))); ++j) {
											fout_html_report << "&nbsp";
										}
									}
									fout_html_report << std::fixed << std::setprecision(2) << itr_homology_data.blastp_evalue_transformed;
									if (itr_homology_data.blastp_evalue_transformed < 1) {
										fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
									}
									else {
										for (auto j = 0; j < (5 - std::floor(std::log10(itr_homology_data.blastp_evalue_transformed))); ++j) {
											fout_html_report << "&nbsp";
										}
									}
									fout_html_report << std::fixed << std::setprecision(3) << itr_homology_data.blastp_parameter_density;
									fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
									fout_html_report << std::fixed << std::setprecision(2) << itr_homology_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average;
									if (itr_homology_data.blastp_evalue_transformed < 1) {
										fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
									}
									else {
										for (auto j = 0; j < (5 - std::floor(std::log10(itr_homology_data.p_peptide_analysis->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average))); ++j) {
											fout_html_report << "&nbsp";
										}
									}
									fout_html_report << std::fixed << std::setprecision(2) << itr_homology_data.p_peptide_analysis->v_denovo_peptide_averagescore;
									if (itr_homology_data.blastp_evalue_transformed < 1) {
										fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
									}
									else {
										for (auto j = 0; j < (5 - std::floor(std::log10(itr_homology_data.p_peptide_analysis->v_denovo_peptide_averagescore))); ++j) {
											fout_html_report << "&nbsp";
										}
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
}
#endif

