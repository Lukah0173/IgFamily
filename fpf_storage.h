//std::ofstream fout_peptidecomparison;
//fout_peptidecomparison.open(itr_v_s_filesystem.directory + "peptidecomparison.txt");
//fout_peptidecomparison << "\n\n\n";
//for (auto itr_1 : main_v_c_parse_csv_proteinpeptides_data) {
//	fout_peptidecomparison << "\n" << itr_1.return_str_parse_proteinpeptides_csv_peptide();
//}
//fout_peptidecomparison << "\n\n\n";
//for (auto itr_2 : main_v_c_parse_csv_denovopeptides_data) {
//	fout_peptidecomparison << "\n" << itr_2.str_parse_PEAKS_denovopeptides_csv_peptide;
//}
//fout_peptidecomparison << "\n\n\n";
//for (auto itr_1 : main_v_c_parse_csv_proteinpeptides_data) {
//	bool itr_1_find = bool();
//	for (auto itr_2 = main_v_c_parse_csv_denovopeptides_data.begin(); itr_2 != main_v_c_parse_csv_denovopeptides_data.end(); ++itr_2) {
//		if (itr_1.return_str_parse_proteinpeptides_csv_peptide() == itr_2->str_parse_PEAKS_denovopeptides_csv_peptide) {
//			itr_1_find = true;
//		}
//	}
//	if (!itr_1_find) {
//		fout_peptidecomparison << "\n" << itr_1.return_str_parse_proteinpeptides_csv_peptide();
//	}
//}
//fout_peptidecomparison << "\n\n\n";
//for (auto itr_2 = main_v_c_parse_csv_denovopeptides_data.begin(); itr_2 != main_v_c_parse_csv_denovopeptides_data.end(); ++itr_2) {
//	bool itr_2_find = bool();
//	for (auto itr_1 = main_v_c_parse_csv_proteinpeptides_data.begin(); itr_1 != main_v_c_parse_csv_proteinpeptides_data.end(); ++itr_1) {
//		if (itr_1->return_str_parse_proteinpeptides_csv_peptide() == itr_2->str_parse_PEAKS_denovopeptides_csv_peptide) {
//			itr_2_find = true;
//		}
//	}
//	if (!itr_2_find) {
//		fout_peptidecomparison << "\n" << itr_2->str_parse_PEAKS_denovopeptides_csv_peptide;
//	}
//}


//size_t count_class_1 = size_t();
//size_t count_class_2 = size_t();
//for (auto itr_count_class = map_main_v_c_analysis_distict.begin(); itr_count_class != map_main_v_c_analysis_distict.end(); ++itr_count_class) {
//	if ((*itr_count_class)->return_str_genefamily_class() == "IGHV") {
//		++count_class_1;
//	}
//	if ((*itr_count_class)->return_str_genefamily_class() == "IGLV") {
//		++count_class_2;
//	}
//}
//
//count_class_1 = size_t();
//count_class_2 = size_t();
//for (auto itr_count_class = map_main_v_c_analysis.begin(); itr_count_class != map_main_v_c_analysis.end(); ++itr_count_class) {
//	if ((*itr_count_class)->return_str_genefamily_class() == "IGHV") {
//		++count_class_1;
//	}
//	if ((*itr_count_class)->return_str_genefamily_class() == "IGLV") {
//		++count_class_2;
//	}
//}

//if (IgFamily::FILESYSTEM_MODE) {
//	vector<fpf_filesystem::filesystem> v_s_filesystem_replicate_combined;
//	for (vector<fpf_filesystem::filesystem>::iterator itr_v_s_filesystem = v_s_filesystem.begin(); itr_v_s_filesystem != v_s_filesystem.end(); ++itr_v_s_filesystem) {
//		if (itr_v_s_filesystem->v_filesystem_replicates.empty()) {
//			v_s_filesystem_replicate_combined.push_back(*itr_v_s_filesystem);
//		}
//		else {
//			bool b_replicate_found = bool();
//			for (auto& itr_v_s_filesystem_replicate_combined : v_s_filesystem_replicate_combined) {
//				for (vector<pair<string, string>>::iterator itr_v_p_replicates = itr_v_s_filesystem->v_filesystem_replicates.begin(); itr_v_p_replicates != itr_v_s_filesystem->v_filesystem_replicates.end(); ++itr_v_p_replicates) {
//					if (itr_v_s_filesystem_replicate_combined.filesystem_id == *itr_v_p_replicates) {
//						b_replicate_found = true;
//						++itr_v_s_filesystem_replicate_combined.filesystem_replicate_count;
//						for (size_t itr_v_c_peptide_data_count = 0; itr_v_c_peptide_data_count < itr_v_s_filesystem->v_c_peptide_data.size(); ++itr_v_c_peptide_data_count) {
//							fpf_data::c_peptide_data con_c_peptide_data = itr_v_s_filesystem->v_c_peptide_data.at(itr_v_c_peptide_data_count);
//							con_c_peptide_data.filesystem_sample_replicate_count = itr_v_s_filesystem_replicate_combined.filesystem_replicate_count;
//							itr_v_s_filesystem_replicate_combined.v_c_peptide_data.push_back(con_c_peptide_data);
//						}
//						if (itr_v_s_filesystem_replicate_combined.filesystem_replicate_count == itr_v_s_filesystem->v_filesystem_replicates.size()) {
//							std::cout << "\n\n\n\n ...combining replicate data";
//							for (vector<pair<string, string>>::iterator itr_v_p_replicates_2 = itr_v_s_filesystem->v_filesystem_replicates.begin(); itr_v_p_replicates_2 != itr_v_s_filesystem->v_filesystem_replicates.end(); ++itr_v_p_replicates_2) {
//								std::cout << "\n\n * " << std::get<0>(*itr_v_p_replicates_2) << "   " << std::get<1>(*itr_v_p_replicates_2);
//							}
//							for (vector<fpf_data::c_genefamily_data>::iterator itr_v_c_analysis = itr_v_s_filesystem_replicate_combined.v_multinomial_category.begin(); itr_v_c_analysis != itr_v_s_filesystem_replicate_combined.v_multinomial_category.end(); ++itr_v_c_analysis) {
//								itr_v_c_analysis->ref_v_c_analysis_polyassociation().clear();
//								itr_v_c_analysis->set_d_score(1);
//							}
//							std::cout << "\n\n ...determining peptide replicate co-occurence";
//							fpf_filesystem::create_v_p_replicate_data(itr_v_s_filesystem_replicate_combined);
//
//							vector<fpf_data::c_peptide_data> filesystem_v_c_peptide_data_filtered = fpf_data::create_v_c_peptide_data_filtered(itr_v_s_filesystem_replicate_combined.v_c_peptide_data);
//							vector<fpf_data::c_peptide_data> filesystem_v_c_peptide_data_distinct = fpf_data::create_v_c_peptide_data_distinct(itr_v_s_filesystem_replicate_combined.v_c_peptide_data);
//							vector<fpf_data::c_peptide_data> filesystem_v_c_peptide_data_filtered_distinct = fpf_data::create_v_c_peptide_data_filtered_distinct(filesystem_v_c_peptide_data_filtered);
//
//							itr_v_s_filesystem_replicate_combined.v_c_multinomial_catagory_distinct = fpf_data::create_v_c_analysis_distinct(itr_v_s_filesystem_replicate_combined.v_multinomial_category);
//							create_global_score_mean(itr_v_s_filesystem_replicate_combined.v_multinomial_category);
//
//							if (SIMPLE_SCORE) {
//								std::cout << "\n\n ...assigning peptides to gene families";
//								fpf_data::create_v_c_analysis_v_c_peptide_data(itr_v_s_filesystem_replicate_combined.v_multinomial_category, itr_v_s_filesystem_replicate_combined.v_c_peptide_data);
//								fpf_data::create_v_c_analysis_v_c_peptide_data_distinct_filtered(itr_v_s_filesystem_replicate_combined.v_multinomial_category, filesystem_v_c_peptide_data_filtered_distinct);
//								fpf_data::create_v_c_analysis_str_alignment(itr_v_s_filesystem_replicate_combined.v_multinomial_category);
//								std::cout << "\n\n ...determining sequence coverage and total spectral count";
//								fpf_data::create_v_c_analysis_st_totalspectralcount(itr_v_s_filesystem_replicate_combined.v_multinomial_category);
//								fpf_data::create_v_c_analysis_d_coverage(itr_v_s_filesystem_replicate_combined.v_multinomial_category);
//								std::cout << "\n\n ...determining peptide association co-occurence";
//								fpf_data::create_v_c_peptide_v_p_peptideassociation(itr_v_s_filesystem_replicate_combined.v_multinomial_category, itr_v_s_filesystem_replicate_combined.v_c_peptide_data);
//								fpf_data::create_v_c_peptide_v_str_peptideassociation_distinct(itr_v_s_filesystem_replicate_combined.v_c_peptide_data);
//								fpf_data::create_v_c_peptide_v_p_peptideassociation_distinct(itr_v_s_filesystem_replicate_combined.v_c_multinomial_catagory_distinct, itr_v_s_filesystem_replicate_combined.v_c_peptide_data);
//								fpf_data::create_v_c_analysis_v_c_peptide_data(itr_v_s_filesystem_replicate_combined.v_multinomial_category, itr_v_s_filesystem_replicate_combined.v_c_peptide_data);
//								std::cout << "\n\n ...calculating score";
//								fpf_data::create_v_c_analysis_d_score(itr_v_s_filesystem_replicate_combined.v_c_peptide_data, itr_v_s_filesystem_replicate_combined.v_multinomial_category, itr_v_s_filesystem_replicate_combined.v_c_multinomial_catagory_distinct);
//								fpf_data::create_v_c_analysis_v_c_peptide_data(itr_v_s_filesystem_replicate_combined.v_multinomial_category, itr_v_s_filesystem_replicate_combined.v_c_peptide_data);
//								std::cout << "\n\n ...training score";
//								fpf_data::train_v_c_analysis_d_score(itr_v_s_filesystem_replicate_combined.v_c_peptide_data, itr_v_s_filesystem_replicate_combined.v_multinomial_category, itr_v_s_filesystem_replicate_combined.v_c_multinomial_catagory_distinct);
//								std::cout << "\n\n ...formatting output";
//								fpf_data::sort_v_c_peptide_data_str_peptide(itr_v_s_filesystem_replicate_combined.v_c_peptide_data);
//								fpf_data::sort_v_c_peptide_data_str_peptide(filesystem_v_c_peptide_data_filtered_distinct);
//								fpf_data::create_v_c_analysis_v_c_peptide_data(itr_v_s_filesystem_replicate_combined.v_multinomial_category, itr_v_s_filesystem_replicate_combined.v_c_peptide_data);
//								fpf_data::update_v_c_analysis_distinct(itr_v_s_filesystem_replicate_combined.v_c_multinomial_catagory_distinct);
//								std::cout << "\n\n ...streaming output";
//								fpf_filesystem_analysis::fout_file_genefamily_summary(itr_v_s_filesystem_replicate_combined);
//								std::cout << "\n\n  * file " << itr_v_s_filesystem_replicate_combined.filename + "_genefamily_summary.txt" << " output";
//								fpf_filesystem_analysis::fout_file_peptide_summary(itr_v_s_filesystem_replicate_combined);
//								std::cout << "\n\n  * file " << itr_v_s_filesystem_replicate_combined.filename + "_peptide_summary.txt" << " output";
//							}
//
//							IgFamily::GLOBAL_ITERATOR = size_t();
//						}
//					}
//				}
//			}
//			if (!b_replicate_found) {
//				v_s_filesystem_replicate_combined.push_back(*itr_v_s_filesystem);
//			}
//		}
//
//if (IgFamily::FILESYSTEM_MODE == 1) {
//	vector<string> v_str_peptide_total_observed = fpf_filesystem_analysis::create_v_peptide_data_total_observed(v_s_filesystem);
//	for (auto itr_v_s_filesystem = v_s_filesystem.begin(); itr_v_s_filesystem != v_s_filesystem.end(); ++itr_v_s_filesystem) {
//		vector<fpf_data::c_peptide_data> filesystem_v_c_peptide_data_filtered = fpf_data::create_v_c_peptide_data_filtered(itr_v_s_filesystem->v_c_peptide_data);
//		vector<fpf_data::c_peptide_data> filesystem_v_c_peptide_data_filtered_distinct = fpf_data::create_v_c_peptide_data_filtered_distinct(filesystem_v_c_peptide_data_filtered);
//		itr_v_s_filesystem->v_c_peptide_data_filtered_distinct = filesystem_v_c_peptide_data_filtered_distinct;
//	}
//	vector<string> v_str_peptide_filtered_distinct_total_observed = fpf_filesystem_analysis::create_v_str_peptide_filtered_distinct_total_observed(v_s_filesystem);
//
//	vector<fpf_filesystem_analysis::filesystem_analysis> v_s_filesystem_analysis;
//	for (auto itr_v_s_filesystem = v_s_filesystem.begin(); itr_v_s_filesystem != v_s_filesystem.end(); ++itr_v_s_filesystem) {
//		v_s_filesystem_analysis.push_back(fpf_filesystem_analysis::create_filesystem_analysis(*itr_v_s_filesystem, v_str_peptide_total_observed, v_str_peptide_filtered_distinct_total_observed));
//	}
//	fpf_filesystem_analysis::fout_filesystem_peptide_data_summary(v_s_filesystem_analysis, v_str_peptide_total_observed, "peptide_summary.csv");
//
//	vector<fpf_filesystem_analysis::filesystem_analysis> v_s_filesystem_analysis_replicate_combined;
//	for (auto itr_v_s_filesystem_replicate_combined = v_s_filesystem_replicate_combined.begin(); itr_v_s_filesystem_replicate_combined != v_s_filesystem_replicate_combined.end(); ++itr_v_s_filesystem_replicate_combined) {
//		v_s_filesystem_analysis_replicate_combined.push_back(fpf_filesystem_analysis::create_filesystem_analysis(*itr_v_s_filesystem_replicate_combined, v_str_peptide_total_observed, v_str_peptide_filtered_distinct_total_observed));
//	}
//	fpf_filesystem_analysis::fout_filesystem_peptide_data_summary(v_s_filesystem_analysis_replicate_combined, v_str_peptide_total_observed, "peptide_summary_replicate_combined.csv");
//}
//	}