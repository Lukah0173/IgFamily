// * * fpf_report.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_REPORT
#define	FPF_REPORT

#include <cstdlib>						// provides - size_t
#include <string>						// provides - std::string
#include <vector>						// provides - std::vector
#include <iomanip>						// provides - std::setprecision
#include <math.h>						// provides - std::log10, std::floor

#include "IgFamily.h"
#include "fpf_filesystem.h"


namespace fpf_report {

	using std::string;
	using std::vector;

	typedef fpf_filesystem::filesystem filesystem;

	void fout_multinomial_comparison(filesystem& par_filesystem) {
		std::string output_multinomial_comparison = par_filesystem.directory + "multinomial_comparison.csv";
		std::ofstream fout_multinomial_comparison;
		fout_multinomial_comparison.open(output_multinomial_comparison);
		for (const auto& itr_v_category_analysis : par_filesystem.v_category_analysis_selected_by_polymorphism) {
			if (itr_v_category_analysis.p_FASTA_category->category_type == "IG") {
				fout_multinomial_comparison << itr_v_category_analysis.p_FASTA_category->category_name << ",";
				fout_multinomial_comparison << itr_v_category_analysis.category_score << "\n";
			}
		}
	}

	void fout_html_report(filesystem& par_filesystem) {
		std::string output_html_report{};
		if (IgFamily::NOVOR_DENOVO) {
			output_html_report = par_filesystem.directory + "NOVOR_report.html";
		}
		else {
			output_html_report = par_filesystem.directory + "report.html";
		}
		std::ofstream fout_html_report;
		fout_html_report.open(output_html_report);
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
		fout_html_report << "&nbsp&nbsp&nbsp" << par_filesystem.enzyme;
		fout_html_report << "&nbsp&nbsp&nbsp" << par_filesystem.denono_deltamass << "&nbspDa";
		for (const auto& itr_category_analysis : par_filesystem.v_category_analysis_selected_by_polymorphism) {
			fout_html_report << "\n\n\n<br><br><br> " << itr_category_analysis.p_FASTA_category->category_name;
			fout_html_report << "&nbsp&nbsp&nbspScore: " << std::fixed << std::setprecision(2) << itr_category_analysis.category_score;
			fout_html_report << "&nbsp&nbsp&nbspCoverage: " << std::fixed << std::setprecision(0) << itr_category_analysis.proteinconstruct_sequencecoverage << "%";
			fout_html_report << "\n\n<br><br> " << itr_category_analysis.p_FASTA_category->category_protein;
			fout_html_report << "\n\n<br> ";
			size_t i{};
			for (const auto& itr_proteinconstruct_from_denovo : itr_category_analysis.v_proteinconstruct_from_denovo) {
				if (itr_proteinconstruct_from_denovo.aminoacid_localconfidence > 0) {
					if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
						if ((itr_proteinconstruct_from_denovo.aminoacid != itr_category_analysis.p_FASTA_category->category_protein.at(i)) && !((itr_proteinconstruct_from_denovo.aminoacid == 'L') && (itr_category_analysis.p_FASTA_category->category_protein.at(i) == 'I'))) {
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
					if (itr_proteinconstruct_from_denovo.aminoacid != itr_category_analysis.p_FASTA_category->category_protein.at(i)) {
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
			for (const auto& itr_proteinconstruct_from_denovo : itr_category_analysis.v_proteinconstruct_from_denovo) {
				if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
					if ((itr_proteinconstruct_from_denovo.aminoacid != itr_category_analysis.p_FASTA_category->category_protein.at(i)) && !((itr_proteinconstruct_from_denovo.aminoacid == 'L') && (itr_category_analysis.p_FASTA_category->category_protein.at(i) == 'I'))) {
						fout_html_report << "<span class=\"mismatch\">";
					}
					if (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed > 35) {
						fout_html_report << "<font color=\"#4c62d6\">";
					}
					if ((itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed <= 35) && (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed > 25)) {
						fout_html_report << "<font color=\"#239B56\">";
					}
					if ((itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed <= 25) && (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed > 15)) {
						fout_html_report << "<font color=\"#E67E22\">";
					}
					if (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed <= 15) {
						fout_html_report << "<font color=\"red\">";
					}
				}
				fout_html_report << itr_proteinconstruct_from_denovo.aminoacid;
				if (itr_proteinconstruct_from_denovo.aminoacid != itr_category_analysis.p_FASTA_category->category_protein.at(i)) {
					fout_html_report << "</span>";
				}
				if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
					fout_html_report << "</font>";
				}
				++i;
			}
			fout_html_report << "&nbsp&nbsp&nbsp" << "Homology";
			fout_html_report << "\n\n<br><br><br>";
			for (const auto& itr_proteinconstruct_from_denovo : itr_category_analysis.v_proteinconstruct_from_denovo) {
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
			fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
			fout_html_report << "DN&nbspconf";
			fout_html_report << "&nbsp&nbsp";
			fout_html_report << "DN&nbspconf&nbspavr";
			fout_html_report << "<br>";
			for (const auto& itr_blastp_data : itr_category_analysis.v_blastp_data_combined_by_category) {
				if (itr_blastp_data.blastp_evalue_transformed > BLASTP_EVALUETRANSFORMED_THRESHOLD) {
					fout_html_report << "\n<br> ";
					int st_mismatch = int();
					for (auto i = 0; i < itr_blastp_data.query_alignment.length(); i) {
						if (itr_blastp_data.query_alignment.at(i) == '.') {
							if (i < itr_blastp_data.p_FASTA_category->category_protein.length()) {
								fout_html_report << ".";
							}
							++i;
						}
						else {
							for (const auto& itr_denovo_aminoacid : itr_blastp_data.p_peptide_data->p_denovo_peptide_best_by_averagelocalconfidence->v_denovo_aminoacid) {
								if (((i >= itr_blastp_data.query_alignment.length()) || (itr_denovo_aminoacid.aminoacid != itr_blastp_data.p_FASTA_category->category_protein.at(i)) && !((itr_denovo_aminoacid.aminoacid == 'L') && (itr_category_analysis.p_FASTA_category->category_protein.at(i) == 'I')))) {
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
								if ((i >= itr_blastp_data.query_alignment.length()) || (itr_denovo_aminoacid.aminoacid != itr_blastp_data.p_FASTA_category->category_protein.at(i))) {
									fout_html_report << "</span>";
								}
								++i;
							}
						}
						st_mismatch = i;
					}
					st_mismatch = (st_mismatch - itr_blastp_data.query_alignment.length());
					for (auto j = 0; j < (5 - st_mismatch); ++j) {
						fout_html_report << "&nbsp";
					}
					fout_html_report << std::fixed << std::setprecision(2) << (itr_blastp_data.blastp_parameter_score * itr_blastp_data.denovo_replicate_count);
					if ((itr_blastp_data.blastp_parameter_score * itr_blastp_data.denovo_replicate_count) < 1) {
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
					}
					else {
						for (auto j = 0; j < (5 - std::floor(std::log10(itr_blastp_data.blastp_parameter_score * itr_blastp_data.denovo_replicate_count))); ++j) {
							fout_html_report << "&nbsp";
						}
					}
					fout_html_report << itr_blastp_data.denovo_replicate_count;
					if (itr_blastp_data.denovo_replicate_count < 1) {
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
					}
					else {
						for (auto j = 0; j < (5 - std::floor(std::log10(itr_blastp_data.denovo_replicate_count))); ++j) {
							fout_html_report << "&nbsp";
						}
					}
					fout_html_report << std::fixed << std::setprecision(2) << itr_blastp_data.blastp_parameter_score;
					if (itr_blastp_data.blastp_parameter_score < 1) {
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp";
					}
					else {
						for (auto j = 0; j < (4 - std::floor(std::log10(itr_blastp_data.blastp_parameter_score))); ++j) {
							fout_html_report << "&nbsp";
						}
					}
					fout_html_report << std::fixed << std::setprecision(2) << itr_blastp_data.blastp_evalue_transformed;
					if (itr_blastp_data.blastp_evalue_transformed < 1) {
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
					}
					else {
						for (auto j = 0; j < (5 - std::floor(std::log10(itr_blastp_data.blastp_evalue_transformed))); ++j) {
							fout_html_report << "&nbsp";
						}
					}
					fout_html_report << std::fixed << std::setprecision(3) << itr_blastp_data.blastp_parameter_density;
					fout_html_report << "&nbsp&nbsp&nbsp&nbsp";
					fout_html_report << std::fixed << std::setprecision(2) << itr_blastp_data.p_peptide_data->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average;
					if (itr_blastp_data.blastp_evalue_transformed < 1) {
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
					}
					else {
						for (auto j = 0; j < (5 - std::floor(std::log10(itr_blastp_data.p_peptide_data->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average))); ++j) {
							fout_html_report << "&nbsp";
						}
					}
					fout_html_report << std::fixed << std::setprecision(2) << itr_blastp_data.p_peptide_data->v_denovo_peptide_averagescore;
					if (itr_blastp_data.blastp_evalue_transformed < 1) {
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
					}
					else {
						for (auto j = 0; j < (5 - std::floor(std::log10(itr_blastp_data.p_peptide_data->v_denovo_peptide_averagescore))); ++j) {
							fout_html_report << "&nbsp";
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

	void fout_html_report_filtered(filesystem& par_filesystem) {
		std::string output_html_report{};
		if (IgFamily::NOVOR_DENOVO) {
			output_html_report = par_filesystem.directory + "NOVOR_report_IG.html";
		}
		else {
			output_html_report = par_filesystem.directory + "report_IG.html";
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
		fout_html_report << "&nbsp&nbsp&nbsp" << par_filesystem.enzyme;
		fout_html_report << "&nbsp&nbsp&nbsp" << par_filesystem.denono_deltamass << "&nbspDa";
		for (const auto& itr_category_analysis : par_filesystem.v_category_analysis_selected_by_polymorphism) {
			if (itr_category_analysis.p_FASTA_category->category_type == "IG") {
				fout_html_report << "\n\n\n<br><br><br> " << itr_category_analysis.p_FASTA_category->category_name;
				fout_html_report << "&nbsp&nbsp&nbspScore: " << std::fixed << std::setprecision(2) << itr_category_analysis.category_score;
				fout_html_report << "&nbsp&nbsp&nbspCoverage: " << std::fixed << std::setprecision(0) << itr_category_analysis.proteinconstruct_sequencecoverage << "%";
				fout_html_report << "\n\n<br><br> " << itr_category_analysis.p_FASTA_category->category_protein;
				fout_html_report << "\n\n<br> ";
				size_t i{};
				for (const auto& itr_proteinconstruct_from_denovo : itr_category_analysis.v_proteinconstruct_from_denovo) {
					if (itr_proteinconstruct_from_denovo.aminoacid_localconfidence > 0) {
						if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
							if ((itr_proteinconstruct_from_denovo.aminoacid != itr_category_analysis.p_FASTA_category->category_protein.at(i)) && !((itr_proteinconstruct_from_denovo.aminoacid == 'L') && (itr_category_analysis.p_FASTA_category->category_protein.at(i) == 'I'))) {
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
						if (itr_proteinconstruct_from_denovo.aminoacid != itr_category_analysis.p_FASTA_category->category_protein.at(i)) {
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
				for (const auto& itr_proteinconstruct_from_denovo : itr_category_analysis.v_proteinconstruct_from_denovo) {
					if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
						if ((itr_proteinconstruct_from_denovo.aminoacid != itr_category_analysis.p_FASTA_category->category_protein.at(i)) && !((itr_proteinconstruct_from_denovo.aminoacid == 'L') && (itr_category_analysis.p_FASTA_category->category_protein.at(i) == 'I'))) {
							fout_html_report << "<span class=\"mismatch\">";
						}
						if (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed > 35) {
							fout_html_report << "<font color=\"#4c62d6\">";
						}
						if ((itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed <= 35) && (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed > 25)) {
							fout_html_report << "<font color=\"#239B56\">";
						}
						if ((itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed <= 25) && (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed > 15)) {
							fout_html_report << "<font color=\"#E67E22\">";
						}
						if (itr_proteinconstruct_from_denovo.aminoacid_evalue_transformed <= 15) {
							fout_html_report << "<font color=\"red\">";
						}
					}
					fout_html_report << itr_proteinconstruct_from_denovo.aminoacid;
					if (itr_proteinconstruct_from_denovo.aminoacid != itr_category_analysis.p_FASTA_category->category_protein.at(i)) {
						fout_html_report << "</span>";
					}
					if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
						fout_html_report << "</font>";
					}
					++i;
				}
				fout_html_report << "&nbsp&nbsp&nbsp" << "Homology";
				fout_html_report << "\n\n<br><br><br>";
				for (const auto& itr_proteinconstruct_from_denovo : itr_category_analysis.v_proteinconstruct_from_denovo) {
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
				for (const auto& itr_blastp_data : itr_category_analysis.v_blastp_data_combined_by_category) {
					if (itr_blastp_data.blastp_evalue_transformed > BLASTP_EVALUETRANSFORMED_THRESHOLD) {
						fout_html_report << "\n<br> ";
						int st_mismatch = int();
						for (auto i = 0; i < itr_blastp_data.query_alignment.length(); i) {
							if (itr_blastp_data.query_alignment.at(i) == '.') {
								if (i < itr_blastp_data.p_FASTA_category->category_protein.length()) {
									fout_html_report << ".";
								}
								++i;
							}
							else {
								for (const auto& itr_denovo_aminoacid : itr_blastp_data.p_peptide_data->p_denovo_peptide_best_by_averagelocalconfidence->v_denovo_aminoacid) {
									if (((i >= itr_blastp_data.query_alignment.length()) || (itr_denovo_aminoacid.aminoacid != itr_blastp_data.p_FASTA_category->category_protein.at(i)) && !((itr_denovo_aminoacid.aminoacid == 'L') && (itr_category_analysis.p_FASTA_category->category_protein.at(i) == 'I')))) {
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
									if ((i >= itr_blastp_data.query_alignment.length()) || (itr_denovo_aminoacid.aminoacid != itr_blastp_data.p_FASTA_category->category_protein.at(i))) {
										fout_html_report << "</span>";
									}
									++i;
								}
							}
							st_mismatch = i;
						}
						st_mismatch = (st_mismatch - itr_blastp_data.query_alignment.length());
						for (auto j = 0; j < (5 - st_mismatch); ++j) {
							fout_html_report << "&nbsp";
						}
						fout_html_report << std::fixed << std::setprecision(2) << (itr_blastp_data.blastp_parameter_score * itr_blastp_data.denovo_replicate_count);
						if ((itr_blastp_data.blastp_parameter_score * itr_blastp_data.denovo_replicate_count) < 1) {
							fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
						}
						else {
							for (auto j = 0; j < (5 - std::floor(std::log10(itr_blastp_data.blastp_parameter_score * itr_blastp_data.denovo_replicate_count))); ++j) {
								fout_html_report << "&nbsp";
							}
						}
						fout_html_report << itr_blastp_data.denovo_replicate_count;
						if (itr_blastp_data.denovo_replicate_count < 1) {
							fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
						}
						else {
							for (auto j = 0; j < (5 - std::floor(std::log10(itr_blastp_data.denovo_replicate_count))); ++j) {
								fout_html_report << "&nbsp";
							}
						}
						fout_html_report << std::fixed << std::setprecision(2) << itr_blastp_data.blastp_parameter_score;
						if (itr_blastp_data.blastp_parameter_score < 1) {
							fout_html_report << "&nbsp&nbsp&nbsp&nbsp";
						}
						else {
							for (auto j = 0; j < (4 - std::floor(std::log10(itr_blastp_data.blastp_parameter_score))); ++j) {
								fout_html_report << "&nbsp";
							}
						}
						fout_html_report << std::fixed << std::setprecision(2) << itr_blastp_data.blastp_evalue_transformed;
						if (itr_blastp_data.blastp_evalue_transformed < 1) {
							fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
						}
						else {
							for (auto j = 0; j < (5 - std::floor(std::log10(itr_blastp_data.blastp_evalue_transformed))); ++j) {
								fout_html_report << "&nbsp";
							}
						}
						fout_html_report << std::fixed << std::setprecision(3) << itr_blastp_data.blastp_parameter_density;
						fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
						fout_html_report << std::fixed << std::setprecision(2) << itr_blastp_data.p_peptide_data->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average;
						if (itr_blastp_data.blastp_evalue_transformed < 1) {
							fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
						}
						else {
							for (auto j = 0; j < (5 - std::floor(std::log10(itr_blastp_data.p_peptide_data->p_denovo_peptide_best_by_averagelocalconfidence->localconfidence_average))); ++j) {
								fout_html_report << "&nbsp";
							}
						}
						fout_html_report << std::fixed << std::setprecision(2) << itr_blastp_data.p_peptide_data->v_denovo_peptide_averagescore;
						if (itr_blastp_data.blastp_evalue_transformed < 1) {
							fout_html_report << "&nbsp&nbsp&nbsp&nbsp&nbsp";
						}
						else {
							for (auto j = 0; j < (5 - std::floor(std::log10(itr_blastp_data.p_peptide_data->v_denovo_peptide_averagescore))); ++j) {
								fout_html_report << "&nbsp";
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

