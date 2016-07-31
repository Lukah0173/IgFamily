// * * fpf_filesystem.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_FILESYSTEM_ANALYSIS
#define	FPF_FILESYSTEM_ANALYSIS
#include <cstdlib> // provides - size_t
#include <vector> // provides - std::vector
#include <iostream> // provides - std::ofstream
#include <algorithm> // provides - std::find, std::find_if, std::sort
#include "IgFamily.h"
#include "fpf_filesystem.h"



namespace fpf_filesystem_analysis {

	struct s_filesystem_analysis;

	typedef std::string string_type;
	typedef size_t size_type;

	struct s_filesystem_analysis {
		fpf_filesystem::s_filesystem* s_filesystem;
		std::vector<fpf_data::peptide_data_type> v_s_peptide_data_analysis;
		std::vector<fpf_data::peptide_data_type> v_s_peptide_data_filtered_distinct_analysis;
	};

	void fout_file_genefamily_summary(fpf_filesystem::s_filesystem par_s_filesystem) {
		std::string output_v_c_analysis = par_s_filesystem.str_directory + "\\" + par_s_filesystem.str_filename + "_genefamily_summary.txt";
		std::ofstream fout_v_c_analysis;
		fout_v_c_analysis.open(output_v_c_analysis);
		fout_v_c_analysis << "-- IgFamily " << IgFamily::version << " --\n\n\n";
		fout_v_c_analysis << "Input file : " << IgFamily::INPUT_CSV << "\n\n\n";

		std::vector<fpf_data::s_multinomial_element_data*> con_map_main_v_c_analysis;
		std::vector<fpf_data::s_multinomial_element_data*> map_main_v_c_analysis = map_v_s_multinomial_element_data_by_score(par_s_filesystem.v_c_analysis_data);
		std::vector<fpf_data::s_multinomial_element_data*> map_main_v_c_analysis_distict = map_v_s_multinomial_element_data_by_score(par_s_filesystem.v_c_analysis_distinct_data);
		if (IgFamily::MAP_FOUT_BY_SCORE == 0) {
			if (IgFamily::MAP_FOUT_BY_DISTINCT == 1) {
				//con_map_main_v_c_analysis = par_s_filesystem.v_c_analysis_data;
			}
			else {
				//con_map_main_v_c_analysis = par_s_filesystem.v_c_analysis_distinct_data;
			}
		}
		if (IgFamily::MAP_FOUT_BY_SCORE == 1) {
			if (IgFamily::MAP_FOUT_BY_DISTINCT == 1) {
				con_map_main_v_c_analysis = map_main_v_c_analysis_distict;
			}
			else {
				con_map_main_v_c_analysis = map_main_v_c_analysis;
			}
		}

		for (std::vector<fpf_data::s_multinomial_element_data*>::iterator itr_v_c_analysis_distinct = con_map_main_v_c_analysis.begin(); itr_v_c_analysis_distinct != con_map_main_v_c_analysis.end(); ++itr_v_c_analysis_distinct) {
			if (itr_v_c_analysis_distinct == con_map_main_v_c_analysis.begin()) {
				fout_v_c_analysis << "\nTop gene families - " << "     Score mean - " << std::fixed << std::setprecision(2) << IgFamily::SCORE_MEAN;
				fout_v_c_analysis << "\n\n";
			}
			if ((itr_v_c_analysis_distinct == con_map_main_v_c_analysis.end()) || ((*itr_v_c_analysis_distinct)->d_score < double{ 0.01 })) {
				break;
			}
			fout_v_c_analysis << "\n  *  " << (*itr_v_c_analysis_distinct)->str_multinomial_element_name;
			for (size_type itr_s_genefamily_output = (*itr_v_c_analysis_distinct)->str_multinomial_element_name.length(); itr_s_genefamily_output < 20; ++itr_s_genefamily_output) {
				fout_v_c_analysis << " ";
			}
			fout_v_c_analysis << "   Score - " << std::fixed << std::setprecision(2) << (*itr_v_c_analysis_distinct)->d_score;
			std::vector<fpf_data::peptide_data_type> con_itr_v_s_peptide_data = (*itr_v_c_analysis_distinct)->v_s_peptide_data;
			double d_peptide_1genefamily = double();
			double d_peptide_2genefamily = double();
			double d_peptide_3genefamily = double();
			for (std::vector<fpf_data::peptide_data_type>::const_iterator itr_v_s_peptide_data = con_itr_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_itr_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
				string_type con_s_genefamily = (*itr_v_c_analysis_distinct)->str_multinomial_element_name;
				size_type sw_s_peptideassociation_distinct = size_type();
				string_type con_s_genefamily_distinct = string_type();
				for (unsigned itr_s_genefamily = 0; itr_s_genefamily < con_s_genefamily.length(); ++itr_s_genefamily) {
					if (con_s_genefamily.at(itr_s_genefamily) == '*') {
						sw_s_peptideassociation_distinct = 1;
					}
					if (sw_s_peptideassociation_distinct == 0) {
						con_s_genefamily_distinct += con_s_genefamily.at(itr_s_genefamily);
					}
				}
				double con_d_score = double();
				double con_d_parameter_gene_family = double();
				double con_d_parameter_gene_family_train = double();
				double con_d_parameter_gene_family_weight = double();
				double con_d_parameter_IgP = double();
				con_d_parameter_gene_family_train = pow((((*itr_v_c_analysis_distinct)->d_score + (double{ 100 } *IgFamily::SCORE_MEAN)) / (double{ 100 } *IgFamily::SCORE_MEAN)), double{ 1.6 });
				con_d_parameter_gene_family_train = 1;
				for (auto itr_v_p_peptideassociation_distinct = itr_v_s_peptide_data->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != itr_v_s_peptide_data->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
					if (std::get<0>(*itr_v_p_peptideassociation_distinct)->str_multinomial_element_name == con_s_genefamily_distinct) {
						con_d_parameter_gene_family_weight = std::get<1>(*itr_v_p_peptideassociation_distinct);
					}
				}
				con_d_parameter_gene_family = pow((con_d_parameter_gene_family_weight / itr_v_s_peptide_data->v_p_peptideassociation_distinct.size()), double{ double{ 1 } / con_d_parameter_gene_family_train });
				con_d_parameter_IgP = (fpf_data::log_basechange(((double(itr_v_s_peptide_data->st_IgP) + double(IgFamily::PARSE_THRESHOLD_IgP)) / (double{ 2 } *double(IgFamily::PARSE_THRESHOLD_IgP))), double{ 9 }) + 1);
				con_d_score += (con_d_parameter_IgP * itr_v_s_peptide_data->st_spectralcount * con_d_parameter_gene_family);
				if (itr_v_s_peptide_data->v_p_peptideassociation_distinct.size() == 1) {
					d_peptide_1genefamily += con_d_score;
				}
				if (itr_v_s_peptide_data->v_p_peptideassociation_distinct.size() == 2) {
					d_peptide_2genefamily += con_d_score;
				}
				if (itr_v_s_peptide_data->v_p_peptideassociation_distinct.size() == 3) {
					d_peptide_3genefamily += con_d_score;
				}
				con_s_genefamily_distinct.clear();
			}
			d_peptide_1genefamily /= (*itr_v_c_analysis_distinct)->d_score;
			d_peptide_2genefamily /= (*itr_v_c_analysis_distinct)->d_score;
			d_peptide_3genefamily /= (*itr_v_c_analysis_distinct)->d_score;
			for (size_type itr_v_s_peptide_data_output = ((log10((*itr_v_c_analysis_distinct)->d_score) >= double{ 1 }) ? (size_type(log10((*itr_v_c_analysis_distinct)->d_score)) + size_type(1)) : size_type(1)); itr_v_s_peptide_data_output < 6; ++itr_v_s_peptide_data_output) {
				fout_v_c_analysis << " ";
			}
			fout_v_c_analysis << "( 1GF - " << std::fixed << std::setprecision(2) << d_peptide_1genefamily;
			fout_v_c_analysis << ", 2GF - " << std::fixed << std::setprecision(2) << d_peptide_2genefamily;
			fout_v_c_analysis << ", 3GF - " << std::fixed << std::setprecision(2) << d_peptide_3genefamily;
			fout_v_c_analysis << " )";
			d_peptide_1genefamily = 0;
			d_peptide_2genefamily = 0;
			d_peptide_3genefamily = 0;
			//for (std::vector<fpf_data::s_multinomial_element_data*>::iterator itr_v_c_analysis_distinct_v_c_polyassociation = (itr_v_c_analysis_distinct + i)->ref_v_c_analysis_polyassociation().begin(); itr_v_c_analysis_distinct_v_c_polyassociation != (itr_v_c_analysis_distinct + i)->ref_v_c_analysis_polyassociation().end(); ++itr_v_c_analysis_distinct_v_c_polyassociation) {
			//	fout_v_c_analysis << "\n - - * " << (*itr_v_c_analysis_distinct_v_c_polyassociation)->str_multinomial_element_name;
			//	for (size_type itr_s_genefamily_output = (itr_v_c_analysis_distinct + i)->str_multinomial_element_name.length(); itr_s_genefamily_output < 20; ++itr_s_genefamily_output) {
			//		fout_v_c_analysis << " ";
			//	}
			//	fout_v_c_analysis << " - " << (*itr_v_c_analysis_distinct_v_c_polyassociation)->d_score;
			//}
		}

		for (std::vector<fpf_data::s_multinomial_element_data*>::iterator itr_v_c_analysis = map_main_v_c_analysis_distict.begin(); itr_v_c_analysis != map_main_v_c_analysis_distict.end(); ++itr_v_c_analysis) {
			fpf_data::sort_c_analysis_v_s_peptide_data_s_peptide((*itr_v_c_analysis)->v_s_peptide_data);
			IgFamily::SCORE_THRESHOLD = (IgFamily::SCORE_MEAN / double{ 2 });
			if (((*itr_v_c_analysis)->v_s_peptide_data.size() != 0) && ((*itr_v_c_analysis)->d_score >= 0.01)) {
				fout_v_c_analysis << "\n\n\n\n\n" << (*itr_v_c_analysis)->str_multinomial_element_name;
				fout_v_c_analysis << "   (top: ";
				double test2 = double();
				string_type test3 = string_type();
				for (std::vector<fpf_data::s_multinomial_element_data*>::iterator itr_v_c_polyassociation = (*itr_v_c_analysis)->v_s_multinomial_element_polyassociation.begin(); itr_v_c_polyassociation != (*itr_v_c_analysis)->v_s_multinomial_element_polyassociation.end(); ++itr_v_c_polyassociation) {
					if ((*itr_v_c_polyassociation)->d_score > test2) {
						test2 = (*itr_v_c_polyassociation)->d_score;
						test3 = (*itr_v_c_polyassociation)->str_multinomial_element_name;
					}
				}
				fout_v_c_analysis << test3 << ")";
				fout_v_c_analysis << "   #SC - " << (*itr_v_c_analysis)->st_totalspectralcount;
				fout_v_c_analysis << "   % - " << std::fixed << std::setprecision(2) << (*itr_v_c_analysis)->d_coverage;
				fout_v_c_analysis << "     Score - " << std::fixed << std::setprecision(2) << (*itr_v_c_analysis)->d_score;
				fout_v_c_analysis << "\n\n" << (*itr_v_c_analysis)->str_protein;
				fout_v_c_analysis << "\n" << (*itr_v_c_analysis)->str_alignment;
				fout_v_c_analysis << "\n";
				std::vector<fpf_data::peptide_data_type> con_itr_v_s_peptide_data = (*itr_v_c_analysis)->v_s_peptide_data;
				for (std::vector<fpf_data::peptide_data_type>::const_iterator itr_v_s_peptide_data = con_itr_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_itr_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
					fout_v_c_analysis << "\n * " << itr_v_s_peptide_data->str_peptide;
					for (size_type itr_v_s_peptide_data_output = itr_v_s_peptide_data->str_peptide.length(); itr_v_s_peptide_data_output < 50; ++itr_v_s_peptide_data_output) {
						fout_v_c_analysis << " ";
					}

					string_type con_s_genefamily = (*itr_v_c_analysis)->str_multinomial_element_name;
					size_type sw_s_peptideassociation_distinct = size_type();
					string_type con_s_genefamily_distinct = string_type();
					for (unsigned itr_s_genefamily = 0; itr_s_genefamily < con_s_genefamily.length(); ++itr_s_genefamily) {
						if (con_s_genefamily.at(itr_s_genefamily) == '*') {
							sw_s_peptideassociation_distinct = 1;
						}
						if (sw_s_peptideassociation_distinct == 0) {
							con_s_genefamily_distinct += con_s_genefamily.at(itr_s_genefamily);
						}
					}
					double con_d_score = double();
					double con_d_parameter_gene_family = double();
					double con_d_parameter_gene_family_weight = double();
					double con_d_parameter_IgP = double();
					for (auto itr_v_p_peptideassociation_distinct = itr_v_s_peptide_data->v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != itr_v_s_peptide_data->v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
						if (std::get<0>(*itr_v_p_peptideassociation_distinct)->str_multinomial_element_name == con_s_genefamily_distinct) {
							con_d_parameter_gene_family_weight = std::get<1>(*itr_v_p_peptideassociation_distinct);
						}
					}
					con_d_parameter_gene_family = (con_d_parameter_gene_family_weight / itr_v_s_peptide_data->v_p_peptideassociation_distinct.size());
					con_d_parameter_IgP = (fpf_data::log_basechange(((double(itr_v_s_peptide_data->st_IgP) + double(IgFamily::PARSE_THRESHOLD_IgP)) / (double{ 2 } *double(IgFamily::PARSE_THRESHOLD_IgP))), double{ 9 }) + 1);
					con_d_score += (con_d_parameter_IgP * itr_v_s_peptide_data->st_spectralcount * con_d_parameter_gene_family);
					con_s_genefamily_distinct.clear();

					fout_v_c_analysis << "   Score - " << std::fixed << std::setprecision(2) << con_d_score;
					for (size_type itr_v_s_peptide_data_output = ((log10(con_d_score) >= double{ 1 }) ? (size_type(log10(con_d_score)) + size_type(1)) : size_type(1)); itr_v_s_peptide_data_output < 6; ++itr_v_s_peptide_data_output) {
						fout_v_c_analysis << " ";
					}
					string_type str_fout_ws = "              ";
					for (auto itr_v_b_replicate_found = itr_v_s_peptide_data->v_p_replicate_data.begin(); itr_v_b_replicate_found != itr_v_s_peptide_data->v_p_replicate_data.end(); ++itr_v_b_replicate_found) {
						if (itr_v_b_replicate_found == itr_v_s_peptide_data->v_p_replicate_data.begin()) {
							fout_v_c_analysis << " (";
						}
						fout_v_c_analysis << "[" << std::get<1>(*itr_v_b_replicate_found) << "]";
						if ((itr_v_b_replicate_found + 1) != itr_v_s_peptide_data->v_p_replicate_data.end()) {
							fout_v_c_analysis << ",";
						}
						if ((itr_v_b_replicate_found + 1) == itr_v_s_peptide_data->v_p_replicate_data.end()) {
							fout_v_c_analysis << ")";
						}
						if (std::get<1>(*itr_v_b_replicate_found) >= 10) {
							for (auto i = double(1); i < log10(double(std::get<1>(*itr_v_b_replicate_found) + size_type(1))); ++i) {
								if (str_fout_ws != "") {
									str_fout_ws.pop_back();
								}
							}
						}
						if ((itr_v_b_replicate_found + 1) == itr_v_s_peptide_data->v_p_replicate_data.end()) {
							fout_v_c_analysis << str_fout_ws;
						}
					}
					for (size_type itr_v_s_peptide_data_output = ((log10(itr_v_s_peptide_data->st_filesystem_replicate) >= double{ 1 }) ? (size_type(log10(itr_v_s_peptide_data->st_filesystem_replicate)) + size_type(1)) : size_type(1)); itr_v_s_peptide_data_output < 3; ++itr_v_s_peptide_data_output) {
						fout_v_c_analysis << " ";
					}
					fout_v_c_analysis << "   SC - " << itr_v_s_peptide_data->st_spectralcount;
					for (size_type itr_v_s_peptide_data_output = ((itr_v_s_peptide_data->st_spectralcount >= size_type{ 1 }) ? size_type(log10(double(itr_v_s_peptide_data->st_spectralcount))) : size_type(0)); itr_v_s_peptide_data_output < 5; ++itr_v_s_peptide_data_output) {
						fout_v_c_analysis << " ";
					}
					fout_v_c_analysis << "   GF par - " << std::fixed << std::setprecision(2) << con_d_parameter_gene_family;
					for (size_type itr_v_s_peptide_data_output = ((log10(con_d_score) >= double{ 1 }) ? size_type(log10(con_d_score)) : size_type(1)); itr_v_s_peptide_data_output < 4; ++itr_v_s_peptide_data_output) {
						fout_v_c_analysis << " ";
					}
					str_fout_ws = "              ";
					for (auto itr_v_b_replicate_found = itr_v_s_peptide_data->v_p_replicate_data.begin(); itr_v_b_replicate_found != itr_v_s_peptide_data->v_p_replicate_data.end(); ++itr_v_b_replicate_found) {
						if (itr_v_b_replicate_found == itr_v_s_peptide_data->v_p_replicate_data.begin()) {
							fout_v_c_analysis << " (";
						}
						fout_v_c_analysis << "[" << std::get<2>(*itr_v_b_replicate_found) << "]";
						if ((itr_v_b_replicate_found + 1) != itr_v_s_peptide_data->v_p_replicate_data.end()) {
							fout_v_c_analysis << ",";
						}
						if ((itr_v_b_replicate_found + 1) == itr_v_s_peptide_data->v_p_replicate_data.end()) {
							fout_v_c_analysis << ")";
						}
						if (std::get<2>(*itr_v_b_replicate_found) >= 10) {
							for (auto i = double(1); i < log10(double(std::get<2>(*itr_v_b_replicate_found) + size_type(1))); ++i) {
								if (str_fout_ws != "") {
									str_fout_ws.pop_back();
								}
							}
						}
						if ((itr_v_b_replicate_found + 1) == itr_v_s_peptide_data->v_p_replicate_data.end()) {
							fout_v_c_analysis << str_fout_ws;
						}
					}
					fout_v_c_analysis << "   Hom score - " << itr_v_s_peptide_data->st_IgP << "  ";
					for (size_type itr_v_s_peptide_data_output = size_type(log10(double(itr_v_s_peptide_data->st_IgP))); itr_v_s_peptide_data_output < 2; ++itr_v_s_peptide_data_output) {
						fout_v_c_analysis << " ";
					}
					fout_v_c_analysis << "   Hom par - " << std::fixed << std::setprecision(2) << con_d_parameter_IgP;
					for (size_type itr_v_s_peptide_data_output = size_type(log10(con_d_parameter_IgP)); itr_v_s_peptide_data_output < 3; ++itr_v_s_peptide_data_output) {
						fout_v_c_analysis << " ";
					}
					fout_v_c_analysis << "   " << "Member of " << itr_v_s_peptide_data->v_p_peptideassociation_distinct.size() << " gene families - ";
					std::vector<std::pair<fpf_data::s_multinomial_element_data*, double>> con_itr_v_p_peptideassociation_distinct = itr_v_s_peptide_data->v_p_peptideassociation_distinct;
					for (std::vector<std::pair<fpf_data::s_multinomial_element_data*, double>>::iterator itr_v_p_peptideassociation_distinct = con_itr_v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != con_itr_v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
						fout_v_c_analysis << std::get<0>(*itr_v_p_peptideassociation_distinct)->str_multinomial_element_name;
						fout_v_c_analysis << "(" << std::get<1>(*itr_v_p_peptideassociation_distinct) << ")";
						if ((itr_v_p_peptideassociation_distinct + 1) != con_itr_v_p_peptideassociation_distinct.end()) {
							fout_v_c_analysis << ", ";
						}
					}
				}
			}
		}
	}

	void fout_file_peptide_summary(fpf_filesystem::s_filesystem par_s_filesystem) {
		string_type str_peptide_summary = par_s_filesystem.str_directory + "\\" + par_s_filesystem.str_filename + "_peptide_summary.txt";
		std::ofstream fout_file_peptide_summary;
		fout_file_peptide_summary.open(str_peptide_summary);
		fout_file_peptide_summary << "-- IgFamily " << IgFamily::version << " --\n\n\n";
		fout_file_peptide_summary << "Input file : " << IgFamily::INPUT_CSV << "\n\n\n";

		if (IgFamily::MAP_FOUT_PEPTIDE_SUMMARY_BY_SPECTRALCOUNT == 1) {
			fpf_data::sort_v_s_peptide_data_st_spectralcount(par_s_filesystem.v_s_peptide_data);
		}

		for (auto itr_v_s_peptide_data = par_s_filesystem.v_s_peptide_data.begin(); itr_v_s_peptide_data != par_s_filesystem.v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
			fout_file_peptide_summary << "\n * " << itr_v_s_peptide_data->str_peptide;
			for (size_type itr_v_s_peptide_data_output = itr_v_s_peptide_data->str_peptide.length(); itr_v_s_peptide_data_output < 50; ++itr_v_s_peptide_data_output) {
				fout_file_peptide_summary << " ";
			}

			string_type str_fout_ws = "                   ";
			for (auto itr_v_b_replicate_found = itr_v_s_peptide_data->v_p_replicate_data.begin(); itr_v_b_replicate_found != itr_v_s_peptide_data->v_p_replicate_data.end(); ++itr_v_b_replicate_found) {
				if (itr_v_b_replicate_found == itr_v_s_peptide_data->v_p_replicate_data.begin()) {
					fout_file_peptide_summary << " (";
				}
				fout_file_peptide_summary << "[" << std::get<1>(*itr_v_b_replicate_found) << "]";
				if ((itr_v_b_replicate_found + 1) != itr_v_s_peptide_data->v_p_replicate_data.end()) {
					fout_file_peptide_summary << ",";
				}
				if ((itr_v_b_replicate_found + 1) == itr_v_s_peptide_data->v_p_replicate_data.end()) {
					fout_file_peptide_summary << ")";
				}
				if (std::get<1>(*itr_v_b_replicate_found) >= 10) {
					for (auto i = double(1); i < log10(double(std::get<1>(*itr_v_b_replicate_found) + size_type(1))); ++i) {
						if (str_fout_ws != "") {
							str_fout_ws.pop_back();
						}
					}
				}
				if ((itr_v_b_replicate_found + 1) == itr_v_s_peptide_data->v_p_replicate_data.end()) {
					fout_file_peptide_summary << str_fout_ws;
				}
			}
			for (size_type itr_v_s_peptide_data_output = ((log10(itr_v_s_peptide_data->st_filesystem_replicate) >= double{ 1 }) ? (size_type(log10(itr_v_s_peptide_data->st_filesystem_replicate)) + size_type(1)) : size_type(1)); itr_v_s_peptide_data_output < 3; ++itr_v_s_peptide_data_output) {
				fout_file_peptide_summary << " ";
			}
			fout_file_peptide_summary << "   SC - " << itr_v_s_peptide_data->st_spectralcount;
			for (size_type itr_v_s_peptide_data_output = ((itr_v_s_peptide_data->st_spectralcount >= size_type{ 1 }) ? size_type(log10(double(itr_v_s_peptide_data->st_spectralcount))) : size_type(0)); itr_v_s_peptide_data_output < 5; ++itr_v_s_peptide_data_output) {
				fout_file_peptide_summary << " ";
			}
			str_fout_ws = "              ";
			for (auto itr_v_b_replicate_found = itr_v_s_peptide_data->v_p_replicate_data.begin(); itr_v_b_replicate_found != itr_v_s_peptide_data->v_p_replicate_data.end(); ++itr_v_b_replicate_found) {
				if (itr_v_b_replicate_found == itr_v_s_peptide_data->v_p_replicate_data.begin()) {
					fout_file_peptide_summary << " (";
				}
				fout_file_peptide_summary << "[" << std::get<2>(*itr_v_b_replicate_found) << "]";
				if ((itr_v_b_replicate_found + 1) != itr_v_s_peptide_data->v_p_replicate_data.end()) {
					fout_file_peptide_summary << ",";
				}
				if ((itr_v_b_replicate_found + 1) == itr_v_s_peptide_data->v_p_replicate_data.end()) {
					fout_file_peptide_summary << ")";
				}
				if (std::get<2>(*itr_v_b_replicate_found) >= 10) {
					for (auto i = double(1); i < log10(double(std::get<2>(*itr_v_b_replicate_found) + size_type(1))); ++i) {
						if (str_fout_ws != "") {
							str_fout_ws.pop_back();
						}
					}
				}
				if ((itr_v_b_replicate_found + 1) == itr_v_s_peptide_data->v_p_replicate_data.end()) {
					fout_file_peptide_summary << str_fout_ws;
				}
			}
			fout_file_peptide_summary << "   " << "Member of " << itr_v_s_peptide_data->v_p_peptideassociation_distinct.size() << " gene families - ";
			std::vector<std::pair<fpf_data::s_multinomial_element_data*, double>> con_itr_v_p_peptideassociation_distinct = itr_v_s_peptide_data->v_p_peptideassociation_distinct;
			for (auto itr_v_p_peptideassociation_distinct = con_itr_v_p_peptideassociation_distinct.begin(); itr_v_p_peptideassociation_distinct != con_itr_v_p_peptideassociation_distinct.end(); ++itr_v_p_peptideassociation_distinct) {
				fout_file_peptide_summary << std::get<0>(*itr_v_p_peptideassociation_distinct)->str_multinomial_element_name;
				fout_file_peptide_summary << "(" << std::fixed << std::setprecision(2) << std::get<1>(*itr_v_p_peptideassociation_distinct) << ")";
				if ((itr_v_p_peptideassociation_distinct + 1) != con_itr_v_p_peptideassociation_distinct.end()) {
					fout_file_peptide_summary << ", ";
				}
			}
		}
	}

	std::vector<string_type> create_v_str_peptide_total_observed(std::vector<fpf_filesystem::s_filesystem> par_v_s_filesystem) {
		std::vector<string_type> con_v_str_peptide_total_observed;
		for (auto itr_v_s_filesystem = par_v_s_filesystem.begin(); itr_v_s_filesystem != par_v_s_filesystem.end(); ++itr_v_s_filesystem) {
			auto con_v_s_peptide_data = itr_v_s_filesystem->v_s_peptide_data;
			for (auto itr_v_s_peptide_data = con_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
				if (std::find(con_v_str_peptide_total_observed.begin(), con_v_str_peptide_total_observed.end(), itr_v_s_peptide_data->str_peptide)
					== con_v_str_peptide_total_observed.end()) {
					con_v_str_peptide_total_observed.push_back(itr_v_s_peptide_data->str_peptide);
				}
			}
		}
		std::sort(con_v_str_peptide_total_observed.begin(), con_v_str_peptide_total_observed.end());
		return con_v_str_peptide_total_observed;
	}

	std::vector<string_type> create_v_str_peptide_filtered_distinct_total_observed(std::vector<fpf_filesystem::s_filesystem> par_v_s_filesystem) {
		std::vector<string_type> con_v_str_peptide_filtered_distinct_total_observed;
		for (auto itr_v_s_filesystem = par_v_s_filesystem.begin(); itr_v_s_filesystem != par_v_s_filesystem.end(); ++itr_v_s_filesystem) {
			auto con_v_s_peptide_data = itr_v_s_filesystem->v_s_peptide_data_filtered_distinct;
			for (auto itr_v_s_peptide_data = con_v_s_peptide_data.begin(); itr_v_s_peptide_data != con_v_s_peptide_data.end(); ++itr_v_s_peptide_data) {
				if (std::find(con_v_str_peptide_filtered_distinct_total_observed.begin(), con_v_str_peptide_filtered_distinct_total_observed.end(), itr_v_s_peptide_data->str_peptide)
					== con_v_str_peptide_filtered_distinct_total_observed.end()) {
					con_v_str_peptide_filtered_distinct_total_observed.push_back(itr_v_s_peptide_data->str_peptide);
				}
			}
		}
		std::sort(con_v_str_peptide_filtered_distinct_total_observed.begin(), con_v_str_peptide_filtered_distinct_total_observed.end());
		return con_v_str_peptide_filtered_distinct_total_observed;
	}

	s_filesystem_analysis create_s_filesystem_analysis(fpf_filesystem::s_filesystem& par_s_filesystem, std::vector<string_type>& par_v_str_peptide_total_observed, std::vector<string_type>& par_v_str_peptide_filtered_distinct_total_observed) {
		s_filesystem_analysis con_s_filesystem_analysis;
		static fpf_data::peptide_data_type static_default_s_peptide_data = fpf_data::peptide_data_type();
		for (auto itr_v_str_peptide_total_observed = par_v_str_peptide_total_observed.begin(); itr_v_str_peptide_total_observed != par_v_str_peptide_total_observed.end(); ++itr_v_str_peptide_total_observed) {
			auto find_v_s_peptide_data = std::find_if(par_s_filesystem.v_s_peptide_data.begin(), par_s_filesystem.v_s_peptide_data.end(), 
				[itr_v_str_peptide_total_observed](fpf_data::peptide_data_type& par_s_peptide_data) {
				return par_s_peptide_data.str_peptide == *itr_v_str_peptide_total_observed; });
			if (find_v_s_peptide_data != par_s_filesystem.v_s_peptide_data.end()) {
				con_s_filesystem_analysis.v_s_peptide_data_analysis.push_back(*find_v_s_peptide_data);
			}
			else {
				con_s_filesystem_analysis.v_s_peptide_data_analysis.push_back(static_default_s_peptide_data);
			}
		}
		static fpf_data::peptide_data_type static_default_s_peptide_data_filtered_distinct = fpf_data::peptide_data_type();
		for (auto itr_v_str_peptide_filtered_distinct_total_observed = par_v_str_peptide_filtered_distinct_total_observed.begin(); itr_v_str_peptide_filtered_distinct_total_observed != par_v_str_peptide_filtered_distinct_total_observed.end(); ++itr_v_str_peptide_filtered_distinct_total_observed) {
			auto find_v_s_peptide_data_filtered_distinct = std::find_if(par_s_filesystem.v_s_peptide_data_filtered_distinct.begin(), par_s_filesystem.v_s_peptide_data_filtered_distinct.end(),
				[itr_v_str_peptide_filtered_distinct_total_observed](fpf_data::peptide_data_type& par_s_peptide_data_filtered_distinct) {
				return par_s_peptide_data_filtered_distinct.str_peptide == *itr_v_str_peptide_filtered_distinct_total_observed; });
			if (find_v_s_peptide_data_filtered_distinct != par_s_filesystem.v_s_peptide_data_filtered_distinct.end()) {
				con_s_filesystem_analysis.v_s_peptide_data_filtered_distinct_analysis.push_back(*find_v_s_peptide_data_filtered_distinct);
			}
			else {
				con_s_filesystem_analysis.v_s_peptide_data_filtered_distinct_analysis.push_back(static_default_s_peptide_data_filtered_distinct);
			}
		}
		con_s_filesystem_analysis.s_filesystem = &par_s_filesystem;
		return con_s_filesystem_analysis;
	}

	//inline bool predicate_v_s_filesystem_peptide_data_by_st_spectralcount(const std::pair<string_type, std::vector<s_filesystem_peptide_data>>& i, const std::pair<string_type, std::vector<s_filesystem_peptide_data>>& j) {
	//	size_type sum_i_st_spectralcount = size_type();
	//	size_type sum_j_st_spectralcount = size_type();
	//	for (auto itr_i = std::get<1>(i).begin(); itr_i != std::get<1>(i).end(); ++itr_i) {
	//		sum_i_st_spectralcount += itr_i->st_spectralcount;
	//	}
	//	for (auto itr_j = std::get<1>(j).begin(); itr_j != std::get<1>(j).end(); ++itr_j) {
	//		sum_j_st_spectralcount += itr_j->st_spectralcount;
	//	}
	//	return (sum_i_st_spectralcount > sum_j_st_spectralcount);
	//}

	//inline void sort_v_s_peptide_data_st_spectralcount(std::vector<std::pair<string_type, std::vector<s_filesystem_peptide_data>>>& par_v_p_filesystem_peptide_data) {
	//	std::sort(par_v_p_filesystem_peptide_data.begin(), par_v_p_filesystem_peptide_data.end(), predicate_v_s_filesystem_peptide_data_by_st_spectralcount);
	//}

	void fout_filesystem_peptide_data_summary(std::vector<s_filesystem_analysis>& par_v_s_filesystem_analysis, std::vector<string_type> par_v_str_peptide_total_observed, string_type par_str_output) {
		std::string output_v_s_filesystem_analysis_peptide_summary = "summary_data\\" + par_str_output;
		std::ofstream fout_v_s_filesystem_analysis_peptide_summary;
		fout_v_s_filesystem_analysis_peptide_summary.open(output_v_s_filesystem_analysis_peptide_summary);
		
		fout_v_s_filesystem_analysis_peptide_summary << ",";
		for (auto itr_v_s_filesystem_analysis = par_v_s_filesystem_analysis.begin(); itr_v_s_filesystem_analysis != par_v_s_filesystem_analysis.end(); ++itr_v_s_filesystem_analysis) {
			fout_v_s_filesystem_analysis_peptide_summary << itr_v_s_filesystem_analysis->s_filesystem->str_filename << ",";			
		}
		fout_v_s_filesystem_analysis_peptide_summary << std::endl;
		for (auto itr_v_s_filesystem_analysis = par_v_s_filesystem_analysis.begin(); itr_v_s_filesystem_analysis != par_v_s_filesystem_analysis.end(); ++itr_v_s_filesystem_analysis) {
			fout_v_s_filesystem_analysis_peptide_summary << itr_v_s_filesystem_analysis->s_filesystem->str_patientstatus << ",";
		}
		fout_v_s_filesystem_analysis_peptide_summary << std::endl;
		for (auto itr_v_s_filesystem_analysis = par_v_s_filesystem_analysis.begin(); itr_v_s_filesystem_analysis != par_v_s_filesystem_analysis.end(); ++itr_v_s_filesystem_analysis) {
			fout_v_s_filesystem_analysis_peptide_summary << itr_v_s_filesystem_analysis->s_filesystem->st_replicate_count << ",";
		}
		fout_v_s_filesystem_analysis_peptide_summary << std::endl;
		for (auto itr_v_s_peptide_data_analysis = size_type(); itr_v_s_peptide_data_analysis < par_v_s_filesystem_analysis.begin()->v_s_peptide_data_analysis.size(); ++itr_v_s_peptide_data_analysis) {
			fout_v_s_filesystem_analysis_peptide_summary << par_v_str_peptide_total_observed[itr_v_s_peptide_data_analysis];
			fout_v_s_filesystem_analysis_peptide_summary << ",";
			for (auto itr_v_s_filesystem_analysis = par_v_s_filesystem_analysis.begin(); itr_v_s_filesystem_analysis != par_v_s_filesystem_analysis.end(); ++itr_v_s_filesystem_analysis) {				
				fout_v_s_filesystem_analysis_peptide_summary << itr_v_s_filesystem_analysis->v_s_peptide_data_analysis[itr_v_s_peptide_data_analysis].st_spectralcount;
				fout_v_s_filesystem_analysis_peptide_summary << ",";
			}
			fout_v_s_filesystem_analysis_peptide_summary << std::endl;
		}
	}

	void fout_filesystem_summary(std::vector<s_filesystem_analysis>& par_v_s_filesystem_analysis) {
		std::string output_v_s_filesystem_summary = "summary_data\\filesystem_summary.txt";
		std::ofstream fout_v_s_filesystem_summary;
		fout_v_s_filesystem_summary.open(output_v_s_filesystem_summary);

		for (auto itr_v_s_filesystem_analysis = par_v_s_filesystem_analysis.begin(); itr_v_s_filesystem_analysis != par_v_s_filesystem_analysis.end(); ++itr_v_s_filesystem_analysis) {
			fout_v_s_filesystem_summary << "ID: " << std::get<0>(itr_v_s_filesystem_analysis->s_filesystem->p_filesystemid) << "," << std::get<1>(itr_v_s_filesystem_analysis->s_filesystem->p_filesystemid) << ";\n";
			fout_v_s_filesystem_summary << "FILE: " << itr_v_s_filesystem_analysis->s_filesystem->str_filename << ";\n";
			fout_v_s_filesystem_summary << "VERSION: " << IgFamily::version << ";\n";
			fout_v_s_filesystem_summary << "REPLICATES: ";
			for (std::vector<std::pair<string_type, string_type>>::iterator itr_v_p_replicates = itr_v_s_filesystem_analysis->s_filesystem->v_p_replicates.begin(); itr_v_p_replicates != itr_v_s_filesystem_analysis->s_filesystem->v_p_replicates.end(); ++itr_v_p_replicates) {
				fout_v_s_filesystem_summary << std::get<0>(*itr_v_p_replicates) << "," << std::get<1>(*itr_v_p_replicates);
				if ((itr_v_p_replicates + 1) != itr_v_s_filesystem_analysis->s_filesystem->v_p_replicates.end()) {
					fout_v_s_filesystem_summary << ",";
				}
			}
			fout_v_s_filesystem_summary << ";\n";
			fout_v_s_filesystem_summary << "STATUS: " << itr_v_s_filesystem_analysis->s_filesystem->str_patientstatus << ";";
			fout_v_s_filesystem_summary << std::endl << std::endl;
		}
	}

}

#endif