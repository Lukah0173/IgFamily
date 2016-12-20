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

	const string version{ "v0.12.5b" };

	const bool FILESYSTEM_MODE{ 0 };
	const bool FILESYSTEM_UPDATE_ALL{ 1 };
	const bool MAP_FOUT_BY_SCORE{ 1 };
	const bool OUTPUT_FASTA{ 1 };
	const bool BLASTP_BY_SELECTED_PEPTIDE{ 1 };

	const string IGFAMILY_ROOT_DIR{ "filesystem_directory\\Will-Murray-Brown_20160927.txt" };
	const string DEFAULT_IGFAMILY_DIRECTORY{ "" };
	const string DEFAULT_PARAMETERS_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "parameters\\" };
	const string DEFAULT_FASTA_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "FASTA\\" };
	const string DEFAULT_FASTA_MODULE_DIRECTORY{ DEFAULT_FASTA_DIRECTORY + "FASTA_modules\\" };
	const string DEFAULT_BLASTP_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "blast_directory\\" };
	const string DEFAULT_MSCONVERT_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "msconvert\\" };
	const string DEFAULT_NOVOR_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "novor\\win\\" };
	const string DEFAULT_TRANSCRIPT_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "transcript_data\\" };
	const string DEFAULT_GENOME_DIRECTORY{ DEFAULT_IGFAMILY_DIRECTORY + "genome_data\\" };
	const string DEFAULT_INPUT_FASTA{ "IGHV_IGLV_IGKV_CONT_UNIPROT_20160827.fasta" };
	const string DEFAULT_PEPTIDE_ASSIGNMENT_METHOD{ "PEAKS de novo" };

	int OUTPUT_FASTA_ACCESSION_WIDTH{ 60 };
	double DENOVO_PEPTIDE_CONFIDENCE_THRESHOLD{ 50 };
	double DENOVO_LOCAL_CONFIDENCE_MOVING_AVERAGE_THRESHOLD{ 84 };
	double DENOVO_PEPTIDE_SIZE_THRESHOLD{ 5 };
	string SELECT_TYPE_GENE_FAMILIES{ "IGV" };
	size_t SELECT_N_MANY_INITIAL_TRAIN_GENE_FAMILIES{ 500 };
	size_t SELECT_N_MANY_GENE_FAMILIES{ 7 };
	double MULTINOMIAL_ELEMENT_OUTPUT_THRESHOLD{ 0.1 };	
	double REPORT_QUERY_ALIGNMENT_TOTALSCORE_OUTPUT_THRESHOLD{ 0.1 };
	double PROTEIN_CONSTRUCT_PARAMETER_DENSITY_CONJUGATED_THRESHOLD{ 0.1 };
	double PROTEIN_CONSTRUCT_SCORE_THRESHOLD{ 1 };
	double REPORT_QUERY_PARAMETER_DENSITY_CONJUGATED_THRESHOLD{ 0.001 };
	double REPORT_PROTEIN_SCORE_THRESHOLD{ 10 };
	double REPORT_PROTEIN_DENSITY_THRESHOLD{ 0.001 };
	double REPORT_QUERY_SCORE_THRESHOLD{ 0.1 };
	double REPORT_QUERY_PARAMETER_SCORE_DENSITY_THRESHOLD{ 0.001 };
	double REPORT_V_HOMOLOGY_DATA_AGGREGATED_BY_PROTEIN_CONJUGATED_DENSITY_THRESHOLD{ 0.1 };

	double PARAMETER_HOMOLOGY_SCALING_FACTOR{ 100 }; // Default: 100
	double PARAMETER_HOMOLOGY_WEIGHT{ 3.5 }; // >= 1 ~ Default: 3.5
	double PARAMETER_HOMOLOGY_MISMATCH_WEIGHT{ 0.30 }; // Default: 0.30
	double PARAMETER_HOMOLOGY_DELTA_ALIGNMENT_WEIGHT{ 0.95 }; // Default: 0.95
	double PARAMETER_SCORE_CONJUGATION_WEIGHT{ 1.0 }; // Default: 1.0
	double PARAMETER_PRIOR_DISTRIBUTION_WEIGHT{ 0.005 }; // Default: 0.005
	double PARAMETER_DEFAULT_LOGISTIC_CONJUGATION_FACTOR{ 1.0 }; // Default: 1.0
	double PARAMETER_LOGISTIC_CONJUGATION_FACTOR{ PARAMETER_DEFAULT_LOGISTIC_CONJUGATION_FACTOR };
	double PARAMETER_LOGISTIC_ITERATION_FACTOR{ 0.001 };  // Default: 0.001

	double DEFAULT_PROTEIN_SCORE_THRESHOLD{ 3 };
	double PROTEIN_SCORE_THRESHOLD{ DEFAULT_PROTEIN_SCORE_THRESHOLD };
	double HOMOLOGY_DENSITY_THRESHOLD { 0.05 };
	double HOMOLOGY_DENSITY_CONJUGATED_THRESHOLD{ 0.05 };
	double CLUSTER_PROPORTION_THRESHOLD{ 0.80 };
}

#endif
