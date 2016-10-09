// * * fpf_core.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_CORE
#define	FPF_CORE

#include "IgFamily.h"
#include "fpf_convert.h"
#include "fpf_data_analysis.h"
#include "fpf_denovo.h"
#include "fpf_filesystem.h"
#include "fpf_homology_analysis.h"
#include "fpf_multinomial.h"
#include "fpf_report.h"


namespace fpf_core {

	using std::string;
	using std::vector;

	using fpf_filesystem::filesystem;
	using fpf_filesystem::sample_analysis;
	using fpf_parse::csv_data;

	inline void core_perform_wiff_fileconversion(filesystem& par_filesystem) {
		fpf_convert::perform_fileconversion(par_filesystem);
	}

	inline void core_perform_novor_denovo(filesystem& par_filesystem) {
		fpf_denovo::perform_novor_denovo(par_filesystem);
	}

	void core_parse_data(filesystem& par_filesystem, vector<string>& par_v_select_peptide_assignment_method) {
		if ((IgFamily::FILESYSTEM_UPDATE_ALL) || (par_filesystem.fileversion != IgFamily::version)) {
			std::cout << "\n\n\n\n parsing data for file ";
			std::cout << par_filesystem.filename;
			std::cout << "...\n";
			for (const auto& itr_v_select_peptide_assignment_method : par_v_select_peptide_assignment_method) {
				sample_analysis temp_sample_analysis{};
				if (itr_v_select_peptide_assignment_method == "PEAKS database match") {
					vector<csv_data> main_v_csv_PEAKS_database_peptides = fpf_filesystem::parse_filesystem_PEAKS_database_peptides(fpf_filesystem::read_filesystem_PEAKS_database_peptides(par_filesystem.directory));					
					temp_sample_analysis.file_found = fpf_parse::check_csv_PEAKS_database_peptides_empty(main_v_csv_PEAKS_database_peptides);
					temp_sample_analysis.peptide_assignment_method = "PEAKS_database";
					std::cout << "\n\n\n * parsing " << par_filesystem.filename << " PEAKS database matched peptides...";
					temp_sample_analysis.v_csv_data = std::move(main_v_csv_PEAKS_database_peptides);
				}
				if (itr_v_select_peptide_assignment_method == "PEAKS de novo") {
					vector<csv_data> main_v_csv_PEAKS_denovo_peptides = fpf_filesystem::parse_filesystem_PEAKS_denovo_peptides(fpf_filesystem::read_filesystem_PEAKS_denovo_peptides(par_filesystem.directory));
					temp_sample_analysis.file_found = fpf_parse::check_csv_PEAKS_denovo_peptides_empty(main_v_csv_PEAKS_denovo_peptides);
					temp_sample_analysis.peptide_assignment_method = "PEAKS_denono";
					std::cout << "\n\n\n * parsing " << par_filesystem.filename << " PEAKS de novo peptides...";
					temp_sample_analysis.v_csv_data = std::move(main_v_csv_PEAKS_denovo_peptides);
					if (!IgFamily::FILESYSTEM_MODE) {
						par_filesystem.filename = temp_sample_analysis.v_csv_data.begin()->csv_sourcefile;
						std::cout << temp_sample_analysis.v_csv_data.begin()->csv_sourcefile;
					}
				}
				if (itr_v_select_peptide_assignment_method == "NOVOR de novo") {
					vector<csv_data> main_v_csv_NOVOR_denovo_peptides = fpf_filesystem::parse_filesystem_NOVOR_denovo_peptides(fpf_filesystem::read_filesystem_NOVOR_denovo_peptides(par_filesystem.directory));
					temp_sample_analysis.file_found = fpf_parse::check_csv_NOVOR_denovo_peptides_empty(main_v_csv_NOVOR_denovo_peptides);
					temp_sample_analysis.peptide_assignment_method = "NOVOR_denono";
					std::cout << "\n\n\n * parsing " << par_filesystem.filename << " NOVOR de novo matched peptides...";
					temp_sample_analysis.v_csv_data = std::move(main_v_csv_NOVOR_denovo_peptides);
				}
				if (temp_sample_analysis.file_found) {
					std::cout << "\n\n peptides parsed - " << temp_sample_analysis.v_csv_data.size();
					par_filesystem.v_sample_analysis.push_back(std::move(temp_sample_analysis));
				}
				else {
					std::cout << "\n ...no data file found";
				}
			}
		}
	}

	void core_create_data_structures(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::cout << "\n\n\n\n creating data structures for file ";
		std::cout << par_filesystem.filename;
		std::cout << "...";
		par_sample_analysis.v_peptide_data = fpf_data::create_v_peptide_data(par_sample_analysis.v_csv_data);
		par_sample_analysis.v_peptide_data_map = fpf_data::create_v_peptide_data_map(par_sample_analysis.v_peptide_data);
		par_sample_analysis.v_peptide_analysis = fpf_data_analysis::create_v_peptide_analysis(par_sample_analysis.v_peptide_data);
		par_sample_analysis.v_peptide_analysis_map = fpf_data_analysis::create_v_peptide_analysis_map(par_sample_analysis.v_peptide_analysis);
		par_sample_analysis.v_protein_data = fpf_data::create_v_protein_data(par_filesystem.v_FASTA_data);
		par_sample_analysis.v_protein_data_map = fpf_data::create_v_protein_data_map(par_sample_analysis.v_protein_data);
		std::cout << "\n ...data structures assigned";
	}

	bool core_create_FASTA_data_structures(filesystem& par_filesystem, const string& par_select_FASTA) {		
		fpf_parse::custom_FASTA_output(par_select_FASTA);
		std::cout << "\n\n\n * parsing FASTA file... \n";
		par_filesystem.v_FASTA_data = fpf_parse::parse_FASTA(par_select_FASTA);
		if (check_FASTA_file_exists(par_filesystem.v_FASTA_data)) { return false; }
		fpf_parse::output_v_FASTA_data(par_filesystem.v_FASTA_data);
		return true;
	}

	void core_homology_analysis(filesystem& par_filesystem, sample_analysis& par_sample_analysis, bool par_refined) {
		if (!par_refined) {
			std::cout << "\n\n\n analysing homology for file ";
			std::cout << par_filesystem.filename;
			std::cout << "...\n";
			fpf_homology_analysis::create_blastp_input(par_filesystem, par_sample_analysis);
			fpf_homology_analysis::create_blastp_database(par_filesystem, par_sample_analysis);
		}
		else {
			std::cout << "\n\n determining most-probable germline representation...\n";
			fpf_data_analysis::select_protein_analysis_by_score(par_sample_analysis);
			std::cout << "\n\n reanalysing homology for file ";
			std::cout << par_filesystem.filename;
			std::cout << "...\n";
			fpf_homology_analysis::create_blastp_database_refined(par_filesystem, par_sample_analysis);
		}
		std::cout << "\n * * * Calling blastp.exe * * *";
		fpf_homology_analysis::systemcall_blastp(par_filesystem, par_sample_analysis);
		std::cout << "\n\n\n * * * Closing blastp.exe * * *";
		std::cout << "\n\n ...homology analysis complete";
		std::cout << "\n\n\n creating homology data structures for file ";
		std::cout << par_filesystem.filename;
		std::cout << "...\n";
		fpf_homology_analysis::parse_homology_data(par_filesystem, par_sample_analysis);
		size_t count_blastp_fail{};
		while (par_sample_analysis.v_homology_data.empty() && (count_blastp_fail <= size_t(3))) {
			++count_blastp_fail;
			std::cout << "\n\n did blastp have a runtime error? ";
			std::cout << "\n\n trying again... \n\n";
			std::cout << "\n * * * Calling blastp.exe * * *";
			fpf_homology_analysis::systemcall_blastp(par_filesystem, par_sample_analysis);
			std::cout << "\n\n\n * * * Closing blastp.exe * * *";
			std::cout << "\n\n ...homology analysis complete";
			std::cout << "\n\n\n creating homology data structures for file ";
			std::cout << par_filesystem.filename;
			std::cout << "...\n";
			fpf_homology_analysis::parse_homology_data(par_filesystem, par_sample_analysis);
		}
		fpf_homology_analysis::modify_homology_data(par_sample_analysis);
		if (!par_refined) {
			fpf_homology_analysis::associate_homology_data_to_protein_data(par_sample_analysis);
		}
		else {
			fpf_homology_analysis::associate_homology_data_to_protein_analysis_refined(par_sample_analysis);
		}
		fpf_homology_analysis::associate_homology_data_to_v_peptide_data(par_sample_analysis);
		fpf_homology_analysis::create_blastp_query_alignment(par_sample_analysis);
		fpf_homology_analysis::transform_homology_data(par_sample_analysis);
		fpf_homology_analysis::determine_homology_data_parameters(par_sample_analysis, false);
		std::cout << " ...homology data structures assigned";
		std::cout << "\n\n outputting homology summary for file ";
		std::cout << par_filesystem.filename;;
		std::cout << "...\n";
		fpf_report::fout_v_homology_data(par_filesystem, par_sample_analysis);
		std::cout << " ...homology file output\n";
	}

	void core_data_analysis(filesystem& par_filesystem, sample_analysis& par_sample_analysis, bool par_refined) {
		std::cout << "\n\n scoring proteins...\n";
		fpf_data_analysis::create_v_protein_analysis(par_sample_analysis, 0, par_refined);
		std::cout << " ...proteins scored\n";
		if (!par_refined) {
			fpf_data_analysis::train_homology_analysis_parameter_score(par_filesystem, par_sample_analysis, IgFamily::SELECT_N_MANY_GENE_FAMILIES_INITIAL_TRAIN, par_refined);
		}
		else {
			fpf_data_analysis::train_homology_analysis_parameter_score(par_filesystem, par_sample_analysis, IgFamily::SELECT_N_MANY_GENE_FAMILIES, par_refined);
			fpf_data_analysis::determine_protein_score_density(par_sample_analysis);
			fpf_data_analysis::create_proteinconstruct_from_denovo(par_sample_analysis);
			fpf_data_analysis::determine_sequence_coverage(par_sample_analysis);
		}
		for (auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			//fpf_data_analysis::sort_v_homology_data_combined_by_protein(itr_v_protein_analysis.v_homology_data_combined_by_protein);
		}
		fpf_data_analysis::sort_v_protein_analysis(par_sample_analysis.v_protein_analysis);
	}

	void core_multinomial(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::cout << "\n\n creating multinomial data frames...\n";
		fpf_multinomial::create_multinomial_data(par_sample_analysis);
	}

	void core_report(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::cout << "\n\n\n producing summary reports...\n";
		std::cout << " ...outputting data reports for " << par_filesystem.filename;
		std::cout << "\n";
		fpf_report::fout_v_peptide_data(par_filesystem, par_sample_analysis);
		fpf_report::fout_v_protein_data(par_filesystem, par_sample_analysis);
		fpf_report::fout_v_peptide_analysis(par_filesystem, par_sample_analysis);
		fpf_report::fout_v_protein_analysis(par_filesystem, par_sample_analysis);
		std::cout << " ...outputting multinomial reports for " << par_filesystem.filename;
		std::cout << "\n ...outputting multinomial data frame for " << par_filesystem.filename;
		fpf_report::fout_multinomial(par_filesystem, par_sample_analysis);
		std::cout << "\n ...outputting multinomial peptide list for " << par_filesystem.filename;
		fpf_report::fout_multinomial_element(par_filesystem, par_sample_analysis);
		std::cout << "\n ...outputting filtered multinomial peptide list for " << par_filesystem.filename;
		fpf_report::fout_multinomial_element_nomatch(par_filesystem, par_sample_analysis);
		std::cout << "\n ...outputting contaminants report for " << par_filesystem.filename;
		fpf_report::fout_multinomial_contaminants_report(par_filesystem, par_sample_analysis);
		std::cout << "\n ...outputting contaminants list for " << par_filesystem.filename;
		fpf_report::fout_multinomial_contaminants_list(par_filesystem, par_sample_analysis);
		std::cout << "\n ...outputting protein score comparison for " << par_filesystem.filename;
		fpf_report::fout_multinomial_protein_score(par_filesystem, par_sample_analysis);
		fpf_report::fout_multinomial_protein_density(par_filesystem, par_sample_analysis);
		fpf_report::fout_protein_pseudoabundance_score(par_filesystem, par_sample_analysis, 0, "allpeptides_");
		fpf_report::fout_protein_pseudoabundance_score(par_filesystem, par_sample_analysis, 0.5, "0.5peptides_");
		fpf_report::fout_protein_pseudoabundance_score(par_filesystem, par_sample_analysis, 0.99, "uniquepeptides_");
		std::cout << "\n ...outputting html report for " << par_filesystem.filename;
		std::cout << "\n";
		fpf_report::fout_html_report(par_filesystem, par_sample_analysis, true, true);
		fpf_report::fout_html_report(par_filesystem, par_sample_analysis, false, true);
		fpf_report::fout_html_report(par_filesystem, par_sample_analysis, true, false);
		fpf_report::fout_html_report(par_filesystem, par_sample_analysis, false, false);
		std::cout << "\n";
		fpf_filesystem::fout_filesystem(par_filesystem);
		IgFamily::PROTEIN_SCORE_THRESHOLD = IgFamily::DEFAULT_PROTEIN_SCORE_THRESHOLD;
		IgFamily::LOGISTIC_CONJUGATION_RANGE = IgFamily::DEFAULT_LOGISTIC_CONJUGATION_RANGE;
		IgFamily::LOGISTIC_CONJUGATION_MIDPOINT = IgFamily::DEFAULT_LOGISTIC_CONJUGATION_MIDPOINT;
	}
}
#endif

