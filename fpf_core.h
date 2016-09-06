// * * fpf_core.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_CORE
#define	FPF_CORE

#include "IgFamily.h"
#include "fpf_data_analysis.h"
#include "fpf_filesystem.h"
#include "fpf_homology_analysis.h"
#include "fpf_multinomial.h"
#include "fpf_report.h"


namespace fpf_core {

	using std::string;
	using std::vector;

	typedef fpf_filesystem::filesystem filesystem;
	typedef fpf_filesystem::sample_analysis sample_analysis;

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
		fpf_homology_analysis::create_v_homology_data(par_filesystem, par_sample_analysis);
		fpf_homology_analysis::modify_filesystem_homology_data(par_sample_analysis);
		fpf_homology_analysis::associate_homology_data_to_v_protein_data(par_sample_analysis);
		fpf_homology_analysis::associate_homology_data_to_v_peptide_data(par_sample_analysis);
		if (par_refined) {
			fpf_homology_analysis::create_protein_from_protein_analysis(par_sample_analysis);
		}
		fpf_homology_analysis::create_query_alignment(par_sample_analysis);
		fpf_homology_analysis::transform_homology_data(par_sample_analysis);
		fpf_homology_analysis::determine_blastp_parameter_density(par_sample_analysis);
		std::cout << " ...homology data structures assigned";
		std::cout << "\n\n outputting homology summary for file ";
		std::cout << par_filesystem.filename;;
		std::cout << "...\n";
		fpf_report::fout_blastp_summary(par_filesystem, par_sample_analysis);
		std::cout << " ...homology file output\n";
	}

	void core_data_analysis(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::cout << "\n\n scoring proteins...\n";
		fpf_data_analysis::create_protein_analysis(par_sample_analysis);
		std::cout << " ...proteins scored\n";
		std::cout << "\n training protein scores...\n";
		for (auto i = 0; i < IgFamily::MULTINOMIAL_CONJUGATION_ITERATION; ++i) {
			fpf_data_analysis::train_homology_analysis_parameter_score(par_filesystem, par_sample_analysis);
			std::cout << " ...iteration ";
			std::cout << (i + 1);
			std::cout << "\n";
		}				
		fpf_data_analysis::create_proteinconstruct_from_denovo(par_sample_analysis);
		fpf_data_analysis::determine_sequence_coverage(par_sample_analysis);
		for (auto& itr_v_protein_analysis : par_sample_analysis.v_protein_analysis) {
			fpf_data_analysis::sort_v_homology_data_with_spectralcount(itr_v_protein_analysis.v_homology_data_combined_by_protein);
		}
		fpf_data_analysis::sort_v_protein_analysis(par_sample_analysis.v_protein_analysis);
	}

	void core_multinomial(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::cout << "\n\n creating multinomial data frames...\n";
		fpf_multinomial::create_filesystem_multinomial_data(par_sample_analysis);
		fpf_report::fout_multinomial(par_filesystem, par_sample_analysis);
		fpf_report::fout_multinomial_element(par_filesystem, par_sample_analysis);
		fpf_report::fout_multinomial_element_nomatch(par_filesystem, par_sample_analysis);
	}

	void core_report(filesystem& par_filesystem, sample_analysis& par_sample_analysis) {
		std::cout << "\n\n\n producing summary reports...\n";
		std::cout << " ...generating peptide report for " << par_filesystem.filename;
		std::cout << "\n";
		fpf_report::fout_v_peptide_data(par_filesystem, par_sample_analysis);
		fpf_report::fout_v_peptide_analysis(par_filesystem, par_sample_analysis);
		fpf_report::fout_v_protein_analysis(par_filesystem, par_sample_analysis);
		std::cout << " ...generating multinomial report for " << par_filesystem.filename;
		std::cout << "\n";
		fpf_report::fout_multinomial_comparison(par_filesystem, par_sample_analysis);
		std::cout << " ...generating html report for " << par_filesystem.filename;
		std::cout << "\n";
		fpf_report::fout_html_report(par_filesystem, par_sample_analysis, true, true);
		fpf_report::fout_html_report(par_filesystem, par_sample_analysis, false, true);
		fpf_report::fout_html_report(par_filesystem, par_sample_analysis, true, false);
		fpf_report::fout_html_report(par_filesystem, par_sample_analysis, false, false);
		std::cout << "\n";
	}
}
#endif

