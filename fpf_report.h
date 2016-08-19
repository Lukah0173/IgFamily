// * * fpf_report.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_REPORT
#define	FPF_REPORT

#include <cstdlib>						// provides - size_t
#include <string>						// provides - std::string
#include <vector>						// provides - vector
#include <iomanip>						// provides - std::setprecision
#include "fpf_data.h"
#include "fpf_filesystem.h"


namespace fpf_report {

	using std::string;
	using std::vector;

	typedef fpf_filesystem::filesystem filesystem;
	typedef fpf_data::multinomial_category multinomial_category;
	typedef fpf_data::peptide_data peptide_data;
	typedef fpf_data::denovo_peptide denovo_peptide;
	typedef fpf_data::blastp_data blastp_data;
	typedef fpf_data::proteinconstruct_from_denovo proteinconstruct_from_denovo;
	typedef fpf_data::category_report category_report;

	void create_category_report(filesystem& par_filesystem) {
		category_report temp_category_report;
		vector<category_report> temp_v_category_report;
		for (const auto itr_blastp_data : par_filesystem.v_blastp_data) {
			auto find_category_report = std::find_if(temp_v_category_report.begin(), temp_v_category_report.end(),
				[itr_blastp_data](const category_report par_category_report) {
				return par_category_report.category_name == itr_blastp_data.blastp_subject_accession;
			});
			if (find_category_report != temp_v_category_report.end()) {
				find_category_report->v_blastp_combined_by_category.push_back(itr_blastp_data);
			}
			else {
				double temp_category_score = double();
				temp_category_report.v_blastp_combined_by_category.push_back(itr_blastp_data);
				temp_category_report.category_class = itr_blastp_data.blastp_subject_accession_class;
				temp_category_report.category_name = itr_blastp_data.blastp_subject_accession;
				temp_category_report.category_protein = itr_blastp_data.category_protein;
				temp_category_report.category_score = temp_category_score;
				temp_v_category_report.push_back(temp_category_report);
				temp_category_report.v_blastp_combined_by_category.clear();
			}
		}
		for (auto& itr_category_report : temp_v_category_report) {
			for (auto& itr_blastp_data : itr_category_report.v_blastp_combined_by_category) {
				const auto find_peptide_data = std::find_if(par_filesystem.v_peptide_data.begin(), par_filesystem.v_peptide_data.end(),
					[itr_blastp_data](const peptide_data par_peptide_data) {
					return par_peptide_data.peptide_filtered == itr_blastp_data.blastp_query;
				});
				if (find_peptide_data == par_filesystem.v_peptide_data.end()) {
					std::cout << "\n\n ERROR: ";
					std::cout << "\n\n find_peptide_data == par_filesystem.v_peptide_data.end()";
					string str_catch_error;
					std::cin >> str_catch_error;
				}
				else {
					itr_blastp_data.denovo_peptide_best_averagelocalconfidence = denovo_peptide();
					for (const auto itr_denovo_peptide : find_peptide_data->v_denovo_peptide_data) {
						if (find_peptide_data->v_denovo_peptide_data.size() == 0) {
							std::cout << "\n\n ERROR: ";
							std::cout << "\n\n find_peptide_data->v_s_denovo_peptide.size() == 0";
							string str_catch_error;
							std::cin >> str_catch_error;
						}
						else {
							++itr_blastp_data.denovo_replicate_count;
							if (itr_blastp_data.denovo_peptide_best_averagelocalconfidence.localconfidence_average < itr_denovo_peptide.localconfidence_average) {
								itr_blastp_data.denovo_peptide_best_averagelocalconfidence = itr_denovo_peptide;
							}
						}
					}
				}
				if (itr_category_report.category_class != "UNIPROT") {
					itr_category_report.category_score += (itr_blastp_data.blastp_evalue_transformed * itr_blastp_data.denovo_replicate_count * 5);
				}
				else {
					itr_category_report.category_score += (itr_blastp_data.blastp_evalue_transformed * itr_blastp_data.denovo_replicate_count);
				}
			}
		}

		par_filesystem.v_category_report = temp_v_category_report;
	}

	void create_proteinconstruct_from_denovo(filesystem& par_filesystem) {
		for (auto& itr_category_report : par_filesystem.v_category_report) {
			for (size_t i = 0; i < itr_category_report.category_protein.length(); ++i) {
				proteinconstruct_from_denovo temp_proteinconstruct_from_denovo;
				temp_proteinconstruct_from_denovo.aminoacid = '.';
				temp_proteinconstruct_from_denovo.aminoacid_score = 0;
				itr_category_report.proteinconstruct_from_denovo.push_back(temp_proteinconstruct_from_denovo);
			}
			vector<blastp_data> v_blastp_query_alignment_selected = vector<blastp_data>();
			vector<blastp_data> v_blastp_query_alignment_rejected = vector<blastp_data>();
			for (const auto itr_v_blastp : itr_category_report.v_blastp_combined_by_category) {		
				const auto find_str_blastp_query_alignment_rejected = std::find_if(v_blastp_query_alignment_rejected.begin(), v_blastp_query_alignment_rejected.end(),
					[itr_v_blastp](const blastp_data par_s_blastp) {
					return par_s_blastp.query_alignment == itr_v_blastp.query_alignment;
				});
				if (find_str_blastp_query_alignment_rejected == v_blastp_query_alignment_rejected.end()) {
					blastp_data temp_blastp_query_alignment = blastp_data();
					for (const auto itr_v_s_blastp_2 : itr_category_report.v_blastp_combined_by_category) {
						for (auto i = 0; i < itr_v_blastp.query_alignment.length(); ++i) {
							if ((itr_v_blastp.query_alignment.at(i) != '.') && (itr_v_s_blastp_2.query_alignment.at(i) != '.')) {
								if (temp_blastp_query_alignment.query_alignment == "") {
									temp_blastp_query_alignment.query_alignment = itr_v_blastp.query_alignment;
									temp_blastp_query_alignment.blastp_evalue_transformed = itr_v_blastp.blastp_evalue_transformed;
									temp_blastp_query_alignment.denovo_replicate_count = itr_v_blastp.denovo_replicate_count;
								}
								if ((itr_v_s_blastp_2.blastp_evalue_transformed * itr_v_s_blastp_2.denovo_replicate_count) >= (temp_blastp_query_alignment.blastp_evalue_transformed * temp_blastp_query_alignment.denovo_replicate_count)) {
									if (temp_blastp_query_alignment.query_alignment != "") {
										v_blastp_query_alignment_rejected.push_back(temp_blastp_query_alignment);
									}
									temp_blastp_query_alignment.query_alignment = itr_v_s_blastp_2.query_alignment;
									temp_blastp_query_alignment.blastp_evalue_transformed = itr_v_s_blastp_2.blastp_evalue_transformed;
									temp_blastp_query_alignment.denovo_replicate_count = itr_v_s_blastp_2.denovo_replicate_count;
								}
								else {
									v_blastp_query_alignment_rejected.push_back(itr_v_s_blastp_2);
								}
								break;
							}
						}
					}
					if (temp_blastp_query_alignment.query_alignment != "") {
						v_blastp_query_alignment_selected.push_back(temp_blastp_query_alignment);
					}
				}
			}
			for (auto i = 0; i < itr_category_report.category_protein.length(); ++i) {
				for (const auto& itr_blastp_query_alignment_selected : v_blastp_query_alignment_selected) {
					if (itr_blastp_query_alignment_selected.query_alignment.at(i) != '.') {
						itr_category_report.proteinconstruct_from_denovo[i].aminoacid = itr_blastp_query_alignment_selected.query_alignment.at(i);
						itr_category_report.proteinconstruct_from_denovo[i].aminoacid_score = itr_blastp_query_alignment_selected.blastp_evalue_transformed;
					}
				}
			}
		}
	}

	inline bool predicate_blastp_data(const blastp_data& i, const blastp_data& j) {
		return (i.blastp_evalue_transformed > j.blastp_evalue_transformed);
	}

	inline void sort_v_blastp_data(vector<blastp_data>& par_v_blastp_data) {
		std::sort(par_v_blastp_data.begin(), par_v_blastp_data.end(), predicate_blastp_data);
	}

	inline bool predicate_blastp_data_with_spectralcount(const blastp_data& i, const blastp_data& j) {
		return ((i.blastp_evalue_transformed * i.denovo_replicate_count) > (j.blastp_evalue_transformed * j.denovo_replicate_count));
	}

	inline void sort_v_blastp_data_with_spectralcount(vector<blastp_data>& par_v_blastp_data) {
		std::sort(par_v_blastp_data.begin(), par_v_blastp_data.end(), predicate_blastp_data_with_spectralcount);
	}

	inline bool predicate_category_report(const category_report& i, const category_report& j) {
		return (i.category_score > j.category_score);
	}

	inline void sort_v_category_report(vector<category_report>& par_v_blastp_data) {
		std::sort(par_v_blastp_data.begin(), par_v_blastp_data.end(), predicate_category_report);
	}

	void fout_html_report(filesystem& par_filesystem) {
		std::string output_html_report = par_filesystem.directory + "report.html";
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
		for (const auto& itr_category_report : par_filesystem.v_category_report) {
			fout_html_report << "\n\n\n<br><br><br> " << itr_category_report.category_name;
			fout_html_report << "     Score: " << std::fixed << std::setprecision(2) << itr_category_report.category_score;
			fout_html_report << "\n\n<br><br> " << itr_category_report.category_protein;
			fout_html_report << "\n\n<br> ";
			size_t foo = size_t();
			for (const auto& itr_proteinconstruct_from_denovo : itr_category_report.proteinconstruct_from_denovo) {
				if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
					if (itr_proteinconstruct_from_denovo.aminoacid != itr_category_report.category_protein.at(foo)) {
						fout_html_report << "<span class=\"mismatch\">";
					}
					if (itr_proteinconstruct_from_denovo.aminoacid_score > 5) {
						fout_html_report << "<font color=\"#239B56\">";
					}
					if ((itr_proteinconstruct_from_denovo.aminoacid_score <= 5) && (itr_proteinconstruct_from_denovo.aminoacid_score > 2)) {
						fout_html_report << "<font color=\"#E67E22\">";
					}
					if (itr_proteinconstruct_from_denovo.aminoacid_score < 2) {
						fout_html_report << "<font color=\"red\">";
					}					
				}
				fout_html_report << itr_proteinconstruct_from_denovo.aminoacid;
				if (itr_proteinconstruct_from_denovo.aminoacid != itr_category_report.category_protein.at(foo)) {
					fout_html_report << "</span>";
				}
				if (itr_proteinconstruct_from_denovo.aminoacid != '.') {
					fout_html_report << "</font>";
				}
				++foo;
			}
			fout_html_report << "\n\n<br><br>";
			for (const auto& itr_blastp_data : itr_category_report.v_blastp_combined_by_category) {
				fout_html_report << "\n<br> ";
				int st_mismatch = int();
				for (auto i = 0; i < itr_blastp_data.query_alignment.length(); i) {
					if (itr_blastp_data.query_alignment.at(i) == '.') {
						if (i < itr_blastp_data.category_protein.length()) {
							fout_html_report << ".";
						}
						++i;
					}
					else {
						for (const auto& itr_denovo_aminoacid : itr_blastp_data.denovo_peptide_best_averagelocalconfidence.v_denovo_aminoacid) {
							if ((i >= itr_blastp_data.query_alignment.length()) || (itr_denovo_aminoacid.aminoacid != itr_blastp_data.category_protein.at(i))) {
								fout_html_report << "<span class=\"mismatch\">";
							}
							if (itr_denovo_aminoacid.aminoacid_score > 80) {
								fout_html_report << "<font color=\"#239B56\">" << itr_denovo_aminoacid.aminoacid << "</font>";
							}
							if ((itr_denovo_aminoacid.aminoacid_score <= 80) && (itr_denovo_aminoacid.aminoacid_score > 60)) {
								fout_html_report << "<font color=\"#E67E22\">" << itr_denovo_aminoacid.aminoacid << "</font>";
							}
							if (itr_denovo_aminoacid.aminoacid_score <= 60) {
								fout_html_report << "<font color=\"red\">" << itr_denovo_aminoacid.aminoacid << "</font>";
							}
							if ((i >= itr_blastp_data.query_alignment.length()) || (itr_denovo_aminoacid.aminoacid != itr_blastp_data.category_protein.at(i))) {
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
				fout_html_report << std::fixed << std::setprecision(2) << itr_blastp_data.blastp_evalue_transformed;
				fout_html_report << "&nbsp&nbsp&nbsp" << itr_blastp_data.denovo_replicate_count;
			}
		}

		fout_html_report << "\
		</p></font>\n \
	</body>\n \
</html>\n ";
	}
}
#endif

