// * * IgFamily.h * *
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef IgFamily
#define	IgFamily

#include <math.h>
#include <string>


namespace IgFamily {

	using std::string;

	const string version{ "v0.12.0e" };

	const bool FILESYSTEM_MODE{ 1 };
	const bool FILESYSTEM_UPDATE_ALL{ 1 };
	const bool MAP_FOUT_BY_SCORE{ 1 };
	const bool MAP_FOUT_BY_DISTINCT{ 1 };
	const bool MAP_FOUT_PEPTIDE_SUMMARY_BY_SPECTRALCOUNT{ 1 };
	const bool OUTPUT_FASTA{ 1 };
	const bool BLASTP_BY_SELECTED_PEPTIDE{ 1 };

	const string IGFAMILY_ROOT_DIR{ "filesystem_directory\\Thilini-Mendis.txt" };
	const string DEFAULT_IGFAMILY_DIRECTORY{ "" };
	const string DEFAULT_FASTA_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "FASTA\\" };
	const string DEFAULT_FASTA_MODULE_DIRECTORY{ DEFAULT_FASTA_DIRECTORY + "FASTA_modules\\" };
	const string DEFAULT_BLASTP_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "blast_directory\\" };
	const string DEFAULT_MSCONVERT_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "msconvert\\" };
	const string DEFAULT_NOVOR_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "novor\\win\\" };
	const string DEFAULT_TRANSCRIPT_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "transcript_data\\" };
	const string DEFAULT_GENOME_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "genome_data\\" };
	const string DEFAULT_INPUT_FASTA{ "IGHV_IGLV_IGKV_CONT_UNIPROT_20160827.fasta" };
	const string DEFAULT_PEPTIDE_ASSIGNMENT_METHOD{ "PEAKS de novo" };

	const int OUTPUT_FASTA_ACCESSION_WIDTH{ 60 };
	const double DENOVO_PEPTIDE_CONFIDENCE_THRESHOLD{ 50 };
	const double DENOVO_LOCAL_CONFIDENCE_MOVING_AVERAGE_THRESHOLD{ 85 };
	const double DENOVO_PEPTIDE_SIZE_THRESHOLD{ 5 };
	const string SELECT_TYPE_GENE_FAMILIES{ "IGV" };
	const size_t SELECT_N_MANY_GENE_FAMILIES_INITIAL_TRAIN{ 500 };
	const size_t SELECT_N_MANY_GENE_FAMILIES{ 10 };
	const double MULTINOMIAL_ELEMENT_OUTPUT_THRESHOLD{ 0.1 };	
	const double REPORT_QUERY_ALIGNMENT_TOTALSCORE_OUTPUT_THRESHOLD{ 0.1 };
	const double PROTEIN_CONSTRUCT_PARAMETER_DENSITY_CONJUGATED_THRESHOLD{ 0.1 };
	const double REPORT_QUERY_PARAMETER_DENSITY_CONJUGATED_THRESHOLD{ 0.1 };
	const double REPORT_PROTEIN_DENSITY_THRESHOLD{ 0.001 };
	const double REPORT_QUERY_PARAMETER_SCORE_DENSITY_THRESHOLD{ 0.001 };
	const double REPORT_V_HOMOLOGY_DATA_AGGREGATED_BY_PROTEIN_COJUGATED_DENSITY_THRESHOLD{ 0.1 };

	const double PARAMETER_HOMOLOGY_WEIGHT{ 2.6 }; // >= 1
	const double PARAMETER_SCORE_CONJUGATION_WEIGHT{ 1.3 };
	const double PARAMETER_SCORE_MISMATCH_WEIGHT{ 0.90 };
	const double PARAMETER_SCORE_COVERAGE_DELTA_WEIGHT{ 0.90 };
	const double DEFAULT_PROTEIN_SCORE_THRESHOLD{ 3 };
	double PROTEIN_SCORE_THRESHOLD{ DEFAULT_PROTEIN_SCORE_THRESHOLD };
	const double HOMOLOGY_SCORE_THRESHOLD_FACTOR{ 50 };
	const double HOMOLOGY_QUERY_ALIGNMENT_COVERAGE_THRESHOLD{ 0 };

	const double DEFAULT_LOGISTIC_CONJUGATION_RANGE{ 0.9995 };
	const double DEFAULT_LOGISTIC_CONJUGATION_MIDPOINT{ (double(1) - DEFAULT_LOGISTIC_CONJUGATION_RANGE) * double(2) };
	double LOGISTIC_CONJUGATION_RANGE{ DEFAULT_LOGISTIC_CONJUGATION_RANGE };
	double LOGISTIC_CONJUGATION_MIDPOINT{ DEFAULT_LOGISTIC_CONJUGATION_MIDPOINT };
	const double LOGISTIC_ITERATION_FACTOR{ DEFAULT_LOGISTIC_CONJUGATION_RANGE / double(100000) };

	const double REPORT_PROTEIN_SCORE_OUTPUT_SCALE{ std::pow(50, PARAMETER_HOMOLOGY_WEIGHT) };
	const double REPORT_PROTEIN_SCORE_THRESHOLD { 1 };
}

#endif
