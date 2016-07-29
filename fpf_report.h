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
#include "fpf_data.h"
#include "fpf_filesystem.h"


namespace fpf_report {
	typedef std::string string_type;
	typedef size_t size_type;
	typedef fpf_data::c_genefamily_data c_genefamily_data_type;
	typedef fpf_filesystem::s_filesystem s_filesystem_type;
	typedef fpf_filesystem::s_report s_report_type;
	typedef fpf_filesystem::s_blastp s_blastp_type;
	typedef fpf_filesystem::s_proteinconstruct_from_denovo s_proteinconstruct_from_denovo_type;

	void create_s_report(s_filesystem_type& par_s_filesystem) {
		s_report_type con_s_report;
		std::vector<s_report_type> con_v_s_report;
		for (auto& itr_v_s_blastp : par_s_filesystem.v_s_blastp) {
			auto find_s_report = std::find_if(con_v_s_report.begin(), con_v_s_report.end(),
				[itr_v_s_blastp](s_report_type par_s_report) {
				return par_s_report.str_protein_accession == itr_v_s_blastp.str_blastp_subject_accession;
			});
			if (find_s_report != con_v_s_report.end()) {
				find_s_report->v_s_blastp.push_back(itr_v_s_blastp);
			} else {
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
			for (auto itr_v_s_blastp : itr_v_s_report.v_s_blastp) {
				itr_v_s_report.d_score += itr_v_s_blastp.d_blastp_par_prop;
			}
		}
		par_s_filesystem.v_s_report = con_v_s_report;
	}

	void create_str_protein_from_denovo(s_filesystem_type& par_s_filesystem) {
		for (auto& itr_v_s_report : par_s_filesystem.v_s_report) {
			for (auto i = 0; i < itr_v_s_report.str_protein.length(); ++i) {
				fpf_filesystem::s_proteinconstruct_from_denovo con_s_proteinconstruct_from_denovo;
				con_s_proteinconstruct_from_denovo.ch_aminoacid = '.';
				con_s_proteinconstruct_from_denovo.d_score = 0;
				itr_v_s_report.s_proteinconstruct_from_denovo.push_back(con_s_proteinconstruct_from_denovo);
			}
			for (auto i = 0; i < itr_v_s_report.str_protein.length(); ++i) {
				for (auto itr_v_s_blastp : itr_v_s_report.v_s_blastp) {
					if ((itr_v_s_blastp.str_blastp_query_alignment.at(i) != '.') && (itr_v_s_blastp.d_blastp_par_prop > itr_v_s_report.s_proteinconstruct_from_denovo[i].d_score)) {
						itr_v_s_report.s_proteinconstruct_from_denovo[i].ch_aminoacid = itr_v_s_blastp.str_blastp_query_alignment.at(i);
						itr_v_s_report.s_proteinconstruct_from_denovo[i].d_score = itr_v_s_blastp.d_blastp_par_prop;
					}
				}
			}
		}
	}

	inline bool predicate_v_s_blastp(const s_blastp_type& i, const s_blastp_type& j) {
		return (i.d_blastp_evalue < j.d_blastp_evalue);
	}

	inline void sort_v_s_blastp(std::vector<s_blastp_type>& par_v_s_blastp) {
		std::sort(par_v_s_blastp.begin(), par_v_s_blastp.end(), predicate_v_s_blastp);
	}

	inline bool predicate_v_s_report(const s_report_type& i, const s_report_type& j) {
		return (i.d_score > j.d_score);
	}

	inline void sort_v_s_report(std::vector<s_report_type>& par_v_s_report) {
		std::sort(par_v_s_report.begin(), par_v_s_report.end(), predicate_v_s_report);
	}

	void fout_s_report(s_filesystem_type& par_s_filesystem) {
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

	void fout_html_report(s_filesystem_type& par_s_filesystem) {
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
			for (auto itr_s_proteinconstruct_from_denovo : itr_v_s_report.s_proteinconstruct_from_denovo) {
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
			//fout_html_report << "\n<br>";
			//for (auto itr_s_proteinconstruct_from_denovo_ident : itr_v_s_report.s_proteinconstruct_from_denovo_ident) {
			//	fout_html_report << itr_s_proteinconstruct_from_denovo_ident.ch_aminoacid;
			//}
			fout_html_report << "\n\n<br><br>";
			for (auto itr_v_s_blastp : itr_v_s_report.v_s_blastp) {
				if (std::find(dummy.begin(), dummy.end(), itr_v_s_blastp.str_blastp_query_alignment) == dummy.end()) {
					fout_html_report << "\n<br> ";
					for (auto i = 0; i < itr_v_s_blastp.str_blastp_query_alignment.length(); ++i) {
						if (itr_v_s_blastp.str_blastp_query_alignment.at(i) == '.') {
							fout_html_report << itr_v_s_blastp.str_blastp_query_alignment.at(i);
						}
						else {
							if ((itr_v_s_blastp.str_blastp_query_alignment.at(i) != itr_v_s_blastp.str_protein.at(i)) && !((itr_v_s_blastp.str_blastp_query_alignment.at(i) != 'L') && (itr_v_s_blastp.str_protein.at(i) == 'I')) && ((itr_v_s_blastp.str_protein.at(i) != 'L') && (itr_v_s_blastp.str_blastp_query_alignment.at(i) == 'I'))) {
								fout_html_report << "<font color=\"red\">" << itr_v_s_blastp.str_blastp_query_alignment.at(i) << "</font>";
							}
							else {
								fout_html_report << itr_v_s_blastp.str_blastp_query_alignment.at(i);
							}
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