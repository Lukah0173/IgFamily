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

int test;

namespace fpf_report {

	typedef std::string string_type;
	typedef size_t size_type;
	typedef fpf_data::multinomial_element_data_type multinomial_element_data_type;
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
				find_s_report->v_s_blastp.push_back(itr_v_s_blastp);
			}
			else {
				double con_d_score = double();
				con_s_report.v_s_blastp.push_back(itr_v_s_blastp);
				con_s_report.str_protein_accession = itr_v_s_blastp.str_blastp_subject_accession;
				con_s_report.str_protein = itr_v_s_blastp.str_protein;
				con_s_report.d_score = con_d_score;
				con_v_s_report.push_back(con_s_report);
				con_s_report.v_s_blastp.clear();
			}
		}
		for (auto& itr_v_s_report : con_v_s_report) {
			for (auto& itr_v_s_blastp : itr_v_s_report.v_s_blastp) {
				itr_v_s_report.d_score += itr_v_s_blastp.d_blastp_par_prop;
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
							if (itr_v_s_blastp.s_denovo_peptide_best.d_denovo_peptide_localconfidence_average < itr_v_s_denovo_peptide.d_denovo_peptide_localconfidence_average) {
								itr_v_s_blastp.s_denovo_peptide_best = itr_v_s_denovo_peptide;
							}
						}
					}
				}
			}
		}

		par_s_filesystem.v_s_report = con_v_s_report;
	}

	void create_str_protein_from_denovo(filesystem_type& par_s_filesystem) {
		for (auto& itr_v_s_report : par_s_filesystem.v_s_report) {
			for (size_type i = 0; i < itr_v_s_report.str_protein.length(); ++i) {
				proteinconstruct_from_denovo_type con_s_proteinconstruct_from_denovo;
				con_s_proteinconstruct_from_denovo.ch_aminoacid = '.';
				con_s_proteinconstruct_from_denovo.d_score = 0;
				itr_v_s_report.proteinconstruct_from_denovo_type.push_back(con_s_proteinconstruct_from_denovo);
			}
			for (size_type i = 0; i < itr_v_s_report.str_protein.length(); ++i) {
				for (auto itr_v_s_blastp : itr_v_s_report.v_s_blastp) {
					if ((itr_v_s_blastp.str_blastp_query_alignment.at(i) != '.') && (itr_v_s_blastp.d_blastp_par_prop > itr_v_s_report.proteinconstruct_from_denovo_type[i].d_score)) {
						itr_v_s_report.proteinconstruct_from_denovo_type[i].ch_aminoacid = itr_v_s_blastp.str_blastp_query_alignment.at(i);
						itr_v_s_report.proteinconstruct_from_denovo_type[i].d_score = itr_v_s_blastp.d_blastp_par_prop;
					}
				}
			}
		}
	}

	inline bool predicate_v_s_blastp(const blastp_type& i, const blastp_type& j) {
		return (i.d_blastp_evalue < j.d_blastp_evalue);
	}

	inline void sort_v_s_blastp(std::vector<blastp_type>& par_v_s_blastp) {
		std::sort(par_v_s_blastp.begin(), par_v_s_blastp.end(), predicate_v_s_blastp);
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
			for (auto itr_v_s_blastp : itr_v_s_report.v_s_blastp) {
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
		for (auto itr_v_s_report : par_s_filesystem.v_s_report) {
			fout_html_report << "\n\n\n<br><br><br> " << itr_v_s_report.str_protein_accession;
			fout_html_report << "     Score: " << std::fixed << std::setprecision(2) << itr_v_s_report.d_score;
			fout_html_report << "\n\n<br><br> " << itr_v_s_report.str_protein;
			fout_html_report << "\n\n<br> ";
			for (auto itr_s_proteinconstruct_from_denovo : itr_v_s_report.proteinconstruct_from_denovo_type) {
				if (itr_s_proteinconstruct_from_denovo.ch_aminoacid != '.') {
					if (itr_s_proteinconstruct_from_denovo.d_score > 3) {
						fout_html_report << "<font color=\"#3498DB\">";
					}
					if ((itr_s_proteinconstruct_from_denovo.d_score <= 3) && (itr_s_proteinconstruct_from_denovo.d_score > 1)) {
						fout_html_report << "<font color=\"#239B56\">";
					}
					if (itr_s_proteinconstruct_from_denovo.d_score < 1) {
						fout_html_report << "<font color=\"#E67E22\">";
					}
				}
				fout_html_report << itr_s_proteinconstruct_from_denovo.ch_aminoacid;
				if (itr_s_proteinconstruct_from_denovo.ch_aminoacid != '.') {
					fout_html_report << "</font>";
				}
			}
			fout_html_report << "\n\n<br><br>";
			for (auto itr_v_s_blastp : itr_v_s_report.v_s_blastp) {
				if (std::find(dummy.begin(), dummy.end(), itr_v_s_blastp.str_blastp_query_alignment) == dummy.end()) {
					fout_html_report << "\n<br> ";
					for (size_type i = 0; i < itr_v_s_blastp.str_blastp_query_alignment.length(); ++i) {
						if (itr_v_s_blastp.str_blastp_query_alignment.at(i) == '.') {
							fout_html_report << itr_v_s_blastp.str_blastp_query_alignment.at(i);
						}
						else {
							size_type st_count_query = size_type();
							//if ((itr_v_s_blastp.str_blastp_query_alignment.at(i) != itr_v_s_blastp.str_protein.at(i)) && !((itr_v_s_blastp.str_blastp_query_alignment.at(i) != 'L') && (itr_v_s_blastp.str_protein.at(i) == 'I')) && ((itr_v_s_blastp.str_protein.at(i) != 'L') && (itr_v_s_blastp.str_blastp_query_alignment.at(i) == 'I'))) {
							//	fout_html_report << "<font color=\"red\">" << itr_v_s_blastp.str_blastp_query_alignment.at(i) << "</font>";
							//}
							//++st_count_query;
							//if (itr_v_s_blastp.)
							//else {
							//	fout_html_report << itr_v_s_blastp.str_blastp_query_alignment.at(i);
							//}
						}
					}
					dummy.push_back(itr_v_s_blastp.str_blastp_query_alignment);
					fout_html_report << "     " << std::fixed << std::setprecision(2) << itr_v_s_blastp.d_blastp_par_prop;
				}
				dummy.clear();
			}
		}

		fout_html_report << "\
		</p></font>\n \
	</body>\n \
</html>\n ";
	}
}
#endif






//		for (auto itr_v_s_filesystem_mnom : itr_v_blastp_filesystem_data.v_s_filesystem_mnom) {
//			auto find_v_str_mnom_xlabel = std::find_if(con_s_model_data.v_str_mnom_xlabel.begin(),
//													   con_s_model_data.v_str_mnom_xlabel.end(),
//													   [itr_v_s_filesystem_mnom](string_type par_str_mnom_xlabel) {
//				return par_str_mnom_xlabel == itr_v_s_filesystem_mnom.str_mnom_comp;
//			});