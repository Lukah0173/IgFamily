// * * fpf_core.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_CORE
#define	FPF_CORE

#include <windows.h>

#include "IgFamily.h"
#include "fpf_convert.h"
#include "fpf_data_analysis.h"
#include "fpf_denovo.h"
#include "fpf_filesystem.h"
#include "fpf_homology_analysis.h"
#include "fpf_multinomial.h"
#include "fpf_parameters.h"
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
			std::cout << "\n\n\n\n parsing parameters file...\n";
			fpf_parameters::read_parameters_file();
			std::cout << "\n ...parameters assigned";
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
					if (!IgFamily::FILESYSTEM_MODE) {
						CreateDirectory((temp_sample_analysis.v_csv_data.begin()->csv_sourcefile + "_" + temp_sample_analysis.peptide_assignment_method).c_str(), NULL);
						par_filesystem.directory = temp_sample_analysis.v_csv_data.begin()->csv_sourcefile + "_" + temp_sample_analysis.peptide_assignment_method + "\\";
						par_filesystem.filename = temp_sample_analysis.v_csv_data.begin()->csv_sourcefile;
						std::cout << temp_sample_analysis.v_csv_data.begin()->csv_sourcefile;
					}
				}
				if (itr_v_select_peptide_assignment_method == "PEAKS de novo") {
					vector<csv_data> main_v_csv_PEAKS_denovo_peptides = fpf_filesystem::parse_filesystem_PEAKS_denovo_peptides(fpf_filesystem::read_filesystem_PEAKS_denovo_peptides(par_filesystem.directory));
					temp_sample_analysis.file_found = fpf_parse::check_csv_PEAKS_denovo_peptides_empty(main_v_csv_PEAKS_denovo_peptides);
					temp_sample_analysis.peptide_assignment_method = "PEAKS_denovo";
					std::cout << "\n\n\n * parsing " << par_filesystem.filename << " PEAKS de novo peptides...";
					temp_sample_analysis.v_csv_data = std::move(main_v_csv_PEAKS_denovo_peptides);
					if (!IgFamily::FILESYSTEM_MODE) {
						CreateDirectory((temp_sample_analysis.v_csv_data.begin()->csv_sourcefile + "_" + temp_sample_analysis.peptide_assignment_method).c_str(), NULL);
						par_filesystem.directory = temp_sample_analysis.v_csv_data.begin()->csv_sourcefile + "_" + temp_sample_analysis.peptide_assignment_method + "\\";
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
					if (!IgFamily::FILESYSTEM_MODE) {
						CreateDirectory((temp_sample_analysis.v_csv_data.begin()->csv_sourcefile + "_" + temp_sample_analysis.peptide_assignment_method).c_str(), NULL);
						par_filesystem.directory = temp_sample_analysis.v_csv_data.begin()->csv_sourcefile + "_" + temp_sample_analysis.peptide_assignment_method + "\\";
						par_filesystem.filename = temp_sample_analysis.v_csv_data.begin()->csv_sourcefile;
						std::cout << temp_sample_analysis.v_csv_data.begin()->csv_sourcefile;
					}
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
		std::cout << "\n\n\n * parsing FASTA file... ";
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
		while (par_sample_analysis.v_homology_data.empty() && (count_blastp_fail <= size_t(5))) {
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
		fpf_homology_analysis::determine_homology_data_uniquenesss(par_sample_analysis);
		fpf_homology_analysis::create_blastp_query_alignment(par_sample_analysis);
		fpf_homology_analysis::transform_homology_data(par_sample_analysis);
		fpf_homology_analysis::determine_HomologyDataParameters(par_sample_analysis, false);
	}

	void core_data_analysis(filesystem& par_filesystem, sample_analysis& par_sample_analysis, bool par_refined) {
		std::cout << "\n\n scoring proteins...\n";
		fpf_data_analysis::create_v_protein_analysis(par_sample_analysis, 0, par_refined, false);
		std::cout << " ...proteins scored\n";
		if (!par_refined) {			
			//fpf_data_analysis::conjugate_homology(par_filesystem, par_sample_analysis, IgFamily::SELECT_N_MANY_INITIAL_TRAIN_GENE_FAMILIES, par_refined);
		}
		else {
			fpf_data_analysis::determine_peptide_sequence_identity(par_sample_analysis);
			fpf_data_analysis::conjugate_homology(par_filesystem, par_sample_analysis, IgFamily::SELECT_N_MANY_GENE_FAMILIES, par_refined);
			//fpf_homology_analysis::determine_HomologyDataParameters(par_sample_analysis, true);
			//fpf_homology_analysis::normalise_v_HomologyData(par_sample_analysis);
			fpf_data_analysis::determine_ProteinScoreDensity(par_sample_analysis);
			fpf_data_analysis::create_ProteinConstruct(par_sample_analysis);
			fpf_data_analysis::determine_SequenceCoverage(par_sample_analysis);
		}
		fpf_data_analysis::sort_v_protein_analysis(par_sample_analysis.v_protein_analysis);
	}

	void core_multinomial(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::cout << "\n\n creating multinomial data frames...\n";
		fpf_multinomial::create_MultinomialData(par_sample_analysis);
	}

	void core_report(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::cout << "\n\n\n producing summary reports...";
		std::cout << "\n ...outputting data reports for ";
		std::cout << par_filesystem.filename;
		fpf_report::fout_v_PeptideData(par_filesystem, par_sample_analysis);
		fpf_report::fout_v_ProteinData(par_filesystem, par_sample_analysis);
		fpf_report::fout_v_PeptideAnalysis(par_filesystem, par_sample_analysis);
		fpf_report::fout_PeptideReport(par_filesystem, par_sample_analysis, "peptide_report");
		fpf_report::fout_v_ProteinAnalysis(par_filesystem, par_sample_analysis);
		fpf_report::fout_v_HomologyData(par_filesystem, par_sample_analysis);
		std::cout << "\n ...outputting homology data report for ";
		std::cout << par_filesystem.filename;
		std::cout << "\n ...outputting multinomial reports for ";
		std::cout << par_filesystem.filename;
		std::cout << "\n ...outputting multinomial data frame for ";
		std::cout << par_filesystem.filename;
		fpf_report::fout_Multinomial(par_filesystem, par_sample_analysis);
		//std::cout << "\n ...outputting multinomial peptide list for ";
		//std::cout << par_filesystem.filename;
		//fpf_report::fout_MultinomialElement(par_filesystem, par_sample_analysis);
		//std::cout << "\n ...outputting filtered multinomial peptide list for ";
		//std::cout << par_filesystem.filename;
		//fpf_report::fout_MultinomialElementNoMatch(par_filesystem, par_sample_analysis);
		//std::cout << "\n ...outputting contaminants list for ";
		//std::cout << par_filesystem.filename;
		//fpf_report::fout_MultinomialContaminantsList(par_filesystem, par_sample_analysis);
		std::cout << "\n ...outputting protein score comparison for ";
		std::cout << par_filesystem.filename;
		fpf_report::fout_MultinomialProteinScore(par_filesystem, par_sample_analysis);
		fpf_report::fout_MultinomialProteinDensity(par_filesystem, par_sample_analysis);
		//fpf_report::fout_ProteinPseudoabundanceScore(par_filesystem, par_sample_analysis, 1, 0, "allpeptides_");
		//fpf_report::fout_ProteinPseudoabundanceScore(par_filesystem, par_sample_analysis, 20, 0.5, "0.5peptides_");
		//fpf_report::fout_ProteinPseudoabundanceScore(par_filesystem, par_sample_analysis, 20, 0.99, "uniquepeptides_");
		std::cout << "\n ...outputting peptide sequence identity report for ";
		std::cout << par_filesystem.filename;
		fpf_report::fout_v_peptide_by_sequence_identity(par_filesystem, par_sample_analysis);
		std::cout << "\n ...outputting html report for ";
		std::cout << par_filesystem.filename;
		std::cout << "\n";
		fpf_homology_analysis::aggregate_v_homology_data_by_homology_distribution(par_sample_analysis);
		fpf_report::fout_HTMLReport(par_filesystem, par_sample_analysis, false, true, 0, "report_allIG_summary");
		fpf_report::fout_HTMLReport(par_filesystem, par_sample_analysis, false, false, 0, "report_allIG");
		fpf_report::fout_HTMLReport(par_filesystem, par_sample_analysis, false, true, IgFamily::REPORT_PROTEIN_SCORE_THRESHOLD, "report_topIG_summary");
		fpf_report::fout_HTMLReport(par_filesystem, par_sample_analysis, false, false, IgFamily::REPORT_PROTEIN_SCORE_THRESHOLD, "report_topIG");
		std::cout << "\n";
		fpf_filesystem::fout_Filesystem(par_filesystem);
		IgFamily::PROTEIN_SCORE_THRESHOLD = IgFamily::DEFAULT_PROTEIN_SCORE_THRESHOLD;
		IgFamily::PARAMETER_LOGISTIC_CONJUGATION_FACTOR = IgFamily::PARAMETER_DEFAULT_LOGISTIC_CONJUGATION_FACTOR;
	}
}
#endif

