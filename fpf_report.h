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
#include "fpf_data.h"
#include "fpf_filesystem.h"



namespace fpf_report {

	typedef std::string string_type;
	typedef size_t size_type;
	typedef fpf_data::multinomial_category_data_type multinomial_category_data_type;
	typedef fpf_data::peptide_data_type peptide_data_type;
	typedef fpf_data::denovo_peptide_type denovo_peptide_type;
	typedef fpf_data::blastp_type blastp_type;
	typedef fpf_data::proteinconstruct_from_denovo_type proteinconstruct_from_denovo_type;
	typedef fpf_data::report_type report_type;
	typedef fpf_filesystem::filesystem_type filesystem_type;

	void create_s_report(filesystem_type& par_s_filesystem) {
		report_type con_s_report;
		std::vector<report_type> con_v_s_report;
		for (auto& itr_v_s_blastp : par_s_filesystem.v_s_blastp) {
			auto find_s_report = std::find_if(con_v_s_report.begin(), con_v_s_report.end(),
				[itr_v_s_blastp](report_type par_s_report) {
				return par_s_report.str_protein_accession == itr_v_s_blastp.str_blastp_subject_accession;
			});
			if (find_s_report != con_v_s_report.end()) {
				find_s_report->v_s_blastp_genefamily_combined.push_back(itr_v_s_blastp);
			}
			else {
				double con_d_score = double();
				con_s_report.v_s_blastp_genefamily_combined.push_back(itr_v_s_blastp);
				con_s_report.str_report_multinomial_category_class = itr_v_s_blastp.str_blastp_subject_accession_class;
				con_s_report.str_protein_accession = itr_v_s_blastp.str_blastp_subject_accession;
				con_s_report.str_protein = itr_v_s_blastp.str_protein;
				con_s_report.d_score = con_d_score;
				con_v_s_report.push_back(con_s_report);
				con_s_report.v_s_blastp_genefamily_combined.clear();
			}
		}
		for (auto& itr_v_s_report : con_v_s_report) {
			for (auto& itr_v_s_blastp : itr_v_s_report.v_s_blastp_genefamily_combined) {
				auto& find_s_peptide_data = std::find_if(par_s_filesystem.v_s_peptide_data.begin(), par_s_filesystem.v_s_peptide_data.end(),
					[itr_v_s_blastp](peptide_data_type par_s_peptide_data) {
					return par_s_peptide_data.str_peptide_filtered == itr_v_s_blastp.str_blastp_query;
				});
				if (find_s_peptide_data == par_s_filesystem.v_s_peptide_data.end()) {
					std::cout << "\n\n ERROR: ";
					std::cout << "\n\n find_s_peptide_data == par_s_filesystem.v_s_peptide_data.end()";
					string_type str_catch_error;
					std::cin >> str_catch_error;
				}
				else {
					itr_v_s_blastp.s_denovo_peptide_best = denovo_peptide_type();
					for (auto itr_v_s_denovo_peptide : find_s_peptide_data->v_s_denovo_peptide) {
						if (find_s_peptide_data->v_s_denovo_peptide.size() == 0) {
							std::cout << "\n\n ERROR: ";
							std::cout << "\n\n find_s_peptide_data->v_s_denovo_peptide.size() == 0";
							string_type str_catch_error;
							std::cin >> str_catch_error;
						}
						else {
							++itr_v_s_blastp.st_count_denovo_replicates;
							if (itr_v_s_blastp.s_denovo_peptide_best.d_denovo_peptide_localconfidence_average < itr_v_s_denovo_peptide.d_denovo_peptide_localconfidence_average) {
								itr_v_s_blastp.s_denovo_peptide_best = itr_v_s_denovo_peptide;
							}
						}
					}
				}
				if (itr_v_s_report.str_report_multinomial_category_class != "UNIPROT") {
					itr_v_s_report.d_score += (itr_v_s_blastp.d_blastp_par_prop * itr_v_s_blastp.st_count_denovo_replicates * 5);
				}
				else {
					itr_v_s_report.d_score += (itr_v_s_blastp.d_blastp_par_prop * itr_v_s_blastp.st_count_denovo_replicates);
				}
			}
		}

		par_s_filesystem.v_s_report = con_v_s_report;
	}

	void create_str_proteinconstruct_from_denovo(filesystem_type& par_s_filesystem) {
		for (auto& itr_v_s_report : par_s_filesystem.v_s_report) {
			for (size_type i = 0; i < itr_v_s_report.str_protein.length(); ++i) {
				proteinconstruct_from_denovo_type con_s_proteinconstruct_from_denovo;
				con_s_proteinconstruct_from_denovo.ch_aminoacid = '.';
				con_s_proteinconstruct_from_denovo.d_score = 0;
				itr_v_s_report.proteinconstruct_from_denovo_type.push_back(con_s_proteinconstruct_from_denovo);
			}
			for (size_type i = 0; i < itr_v_s_report.str_protein.length(); ++i) {
				for (auto itr_v_s_blastp : itr_v_s_report.v_s_blastp_genefamily_combined) {
					if ((itr_v_s_blastp.str_blastp_query_alignment.at(i) != '.') && ((itr_v_s_blastp.d_blastp_par_prop * itr_v_s_blastp.st_count_denovo_replicates) > itr_v_s_report.proteinconstruct_from_denovo_type[i].d_score)) {
						itr_v_s_report.proteinconstruct_from_denovo_type[i].ch_aminoacid = itr_v_s_blastp.str_blastp_query_alignment.at(i);
						itr_v_s_report.proteinconstruct_from_denovo_type[i].d_score = (itr_v_s_blastp.d_blastp_par_prop * itr_v_s_blastp.st_count_denovo_replicates);
					}
				}
			}
		}
	}

	inline bool predicate_v_s_blastp(const blastp_type& i, const blastp_type& j) {
		return (i.d_blastp_par_prop > j.d_blastp_par_prop);
	}

	inline void sort_v_s_blastp(std::vector<blastp_type>& par_v_s_blastp) {
		std::sort(par_v_s_blastp.begin(), par_v_s_blastp.end(), predicate_v_s_blastp);
	}

	inline bool predicate_v_s_blastp_withspectralcount(const blastp_type& i, const blastp_type& j) {
		return ((i.d_blastp_par_prop * i.st_count_denovo_replicates) > (j.d_blastp_par_prop * j.st_count_denovo_replicates));
	}

	inline void sort_v_s_blastp_withspectralcount(std::vector<blastp_type>& par_v_s_blastp) {
		std::sort(par_v_s_blastp.begin(), par_v_s_blastp.end(), predicate_v_s_blastp_withspectralcount);
	}

	inline bool predicate_v_s_report(const report_type& i, const report_type& j) {
		return (i.d_score > j.d_score);
	}

	inline void sort_v_s_report(std::vector<report_type>& par_v_s_report) {
		std::sort(par_v_s_report.begin(), par_v_s_report.end(), predicate_v_s_report);
	}

	void fout_s_report(filesystem_type& par_s_filesystem) {
		std::string output_s_report = par_s_filesystem.str_directory + "report.txt";
		std::ofstream fout_s_report;
		fout_s_report.open(output_s_report);
		std::vector<string_type> dummy;
		std::cout << "\n\n ...generating report for " << par_s_filesystem.str_filename;
		fout_s_report << "\n " << par_s_filesystem.str_filename;
		for (auto itr_v_s_report : par_s_filesystem.v_s_report) {
			fout_s_report << "\n\n\n " << itr_v_s_report.str_protein_accession;
			fout_s_report << "     Score: " << std::fixed << std::setprecision(2) << itr_v_s_report.d_score;
			fout_s_report << "\n\n " << itr_v_s_report.str_protein;
			fout_s_report << "\n";
			for (auto itr_v_s_blastp : itr_v_s_report.v_s_blastp_genefamily_combined) {
				if (std::find(dummy.begin(), dummy.end(), itr_v_s_blastp.str_blastp_query_alignment) == dummy.end()) {
					fout_s_report << "\n " << itr_v_s_blastp.str_blastp_query_alignment;
					dummy.push_back(itr_v_s_blastp.str_blastp_query_alignment);
					fout_s_report << "     " << std::fixed << std::setprecision(2) << itr_v_s_blastp.d_blastp_par_prop;
				}
			}
			dummy.clear();
		}
	}

	void fout_html_report(filesystem_type& par_s_filesystem) {
		std::string output_html_report = par_s_filesystem.str_directory + "report.html";
		std::ofstream fout_html_report;
		fout_html_report.open(output_html_report);
		std::vector<string_type> dummy;
		std::cout << "\n\n ...generating html report for " << par_s_filesystem.str_filename;
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
		fout_html_report << "\n<br>" << par_s_filesystem.str_filename;
		fout_html_report << "\n\n<br><br>" << par_s_filesystem.str_fileversion;
		fout_html_report << "&nbsp&nbsp&nbsp" << par_s_filesystem.str_enzyme;
		fout_html_report << "&nbsp&nbsp&nbsp" << par_s_filesystem.str_denono_deltamass << "&nbspDa";
		for (auto itr_v_s_report : par_s_filesystem.v_s_report) {
			fout_html_report << "\n\n\n<br><br><br> " << itr_v_s_report.str_protein_accession;
			fout_html_report << "     Score: " << std::fixed << std::setprecision(2) << itr_v_s_report.d_score;
			fout_html_report << "\n\n<br><br> " << itr_v_s_report.str_protein;
			fout_html_report << "\n\n<br> ";
			size_type foo = size_type();
			for (auto itr_s_proteinconstruct_from_denovo : itr_v_s_report.proteinconstruct_from_denovo_type) {
				if (itr_s_proteinconstruct_from_denovo.ch_aminoacid != '.') {
					if (itr_s_proteinconstruct_from_denovo.ch_aminoacid != itr_v_s_report.str_protein.at(foo)) {
						fout_html_report << "<span class=\"mismatch\">";
					}
					if (itr_s_proteinconstruct_from_denovo.d_score > 5) {
						fout_html_report << "<font color=\"#239B56\">";
					}
					if ((itr_s_proteinconstruct_from_denovo.d_score <= 5) && (itr_s_proteinconstruct_from_denovo.d_score > 2)) {
						fout_html_report << "<font color=\"#E67E22\">";
					}
					if (itr_s_proteinconstruct_from_denovo.d_score < 2) {
						fout_html_report << "<font color=\"red\">";
					}					
				}
				fout_html_report << itr_s_proteinconstruct_from_denovo.ch_aminoacid;
				if (itr_s_proteinconstruct_from_denovo.ch_aminoacid != itr_v_s_report.str_protein.at(foo)) {
					fout_html_report << "</span>";
				}
				if (itr_s_proteinconstruct_from_denovo.ch_aminoacid != '.') {
					fout_html_report << "</font>";
				}
				++foo;
			}
			fout_html_report << "\n\n<br><br>";
			for (auto itr_v_s_blastp : itr_v_s_report.v_s_blastp_genefamily_combined) {
				fout_html_report << "\n<br> ";
				int st_mismatch = int();
				for (int i = 0; i < itr_v_s_blastp.str_blastp_query_alignment.length(); i) {
					if (itr_v_s_blastp.str_blastp_query_alignment.at(i) == '.') {
						if (i < itr_v_s_blastp.str_protein.length()) {
							fout_html_report << ".";
						}
						++i;
					}
					else {
						for (auto itr_v_s_denovo_aminoacid : itr_v_s_blastp.s_denovo_peptide_best.v_s_denovo_aminoacid) {
							if ((i >= itr_v_s_blastp.str_blastp_query_alignment.length()) || (itr_v_s_denovo_aminoacid.ch_aminoacid != itr_v_s_blastp.str_protein.at(i))) {
								fout_html_report << "<span class=\"mismatch\">";
							}
							if (itr_v_s_denovo_aminoacid.d_denovo_localconfidence > 80) {
								fout_html_report << "<font color=\"#239B56\">" << itr_v_s_denovo_aminoacid.ch_aminoacid << "</font>";
							}
							if ((itr_v_s_denovo_aminoacid.d_denovo_localconfidence <= 80) && (itr_v_s_denovo_aminoacid.d_denovo_localconfidence > 60)) {
								fout_html_report << "<font color=\"#E67E22\">" << itr_v_s_denovo_aminoacid.ch_aminoacid << "</font>";
							}
							if (itr_v_s_denovo_aminoacid.d_denovo_localconfidence <= 60) {
								fout_html_report << "<font color=\"red\">" << itr_v_s_denovo_aminoacid.ch_aminoacid << "</font>";
							}
							if ((i >= itr_v_s_blastp.str_blastp_query_alignment.length()) || (itr_v_s_denovo_aminoacid.ch_aminoacid != itr_v_s_blastp.str_protein.at(i))) {
								fout_html_report << "</span>";
							}
							++i;
						}
					}
					st_mismatch = i;
				}
				st_mismatch = (st_mismatch - itr_v_s_blastp.str_blastp_query_alignment.length());
				for (auto j = 0; j < (5 - st_mismatch); ++j) {
					fout_html_report << "&nbsp";
				}
				fout_html_report << std::fixed << std::setprecision(2) << itr_v_s_blastp.d_blastp_par_prop;
				fout_html_report << "&nbsp&nbsp&nbsp" << itr_v_s_blastp.st_count_denovo_replicates;
			}
		}

		fout_html_report << "\
		</p></font>\n \
	</body>\n \
</html>\n ";
	}
}
#endif

