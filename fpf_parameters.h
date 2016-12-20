// * * fpf_core.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_PARAMETERS
#define	FPF_PARAMETERS

#include <cstdlib>
#include <iostream>
#include <istream>
#include <string>

#include "IgFamily.h"


namespace fpf_parameters {

	using std::string;

	void read_parameters_file()
	{
		string stream_parameters = IgFamily::DEFAULT_PARAMETERS_DIRECTORY + "parameters.txt";
		std::ifstream fin_input_parameters(stream_parameters);
		char read_parameters{};
		string parse_parameters{};
		while (fin_input_parameters.get(read_parameters))
		{
			parse_parameters += read_parameters;
			if (read_parameters == '\n')
			{
				parse_parameters.clear();
			}
			if (parse_parameters == "PARAMETER_HOMOLOGY_SCALING_FACTOR")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{					
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::PARAMETER_HOMOLOGY_SCALING_FACTOR = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "PARAMETER_HOMOLOGY_WEIGHT")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::PARAMETER_HOMOLOGY_WEIGHT = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "PARAMETER_HOMOLOGY_MISMATCH_WEIGHT")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::PARAMETER_HOMOLOGY_MISMATCH_WEIGHT = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "PARAMETER_HOMOLOGY_DELTA_ALIGNMENT_WEIGHT")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::PARAMETER_HOMOLOGY_DELTA_ALIGNMENT_WEIGHT = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "PARAMETER_SCORE_CONJUGATION_WEIGHT")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::PARAMETER_SCORE_CONJUGATION_WEIGHT = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "PARAMETER_PRIOR_DISTRIBUTION_WEIGHT")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::PARAMETER_PRIOR_DISTRIBUTION_WEIGHT = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "PARAMETER_DEFAULT_LOGISTIC_CONJUGATION_FACTOR")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::PARAMETER_DEFAULT_LOGISTIC_CONJUGATION_FACTOR = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "PARAMETER_LOGISTIC_ITERATION_FACTOR")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::PARAMETER_LOGISTIC_ITERATION_FACTOR = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "OUTPUT_FASTA_ACCESSION_WIDTH")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::OUTPUT_FASTA_ACCESSION_WIDTH = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "DENOVO_PEPTIDE_CONFIDENCE_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::DENOVO_PEPTIDE_CONFIDENCE_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "DENOVO_LOCAL_CONFIDENCE_MOVING_AVERAGE_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::DENOVO_LOCAL_CONFIDENCE_MOVING_AVERAGE_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "DENOVO_PEPTIDE_SIZE_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::DENOVO_PEPTIDE_SIZE_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "SELECT_TYPE_GENE_FAMILIES")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::SELECT_TYPE_GENE_FAMILIES = parse_parameter_value;
			}
			if (parse_parameters == "SELECT_N_MANY_INITIAL_TRAIN_GENE_FAMILIES")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::SELECT_N_MANY_INITIAL_TRAIN_GENE_FAMILIES = std::stoi(parse_parameter_value);
			}
			if (parse_parameters == "SELECT_N_MANY_GENE_FAMILIES")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::SELECT_N_MANY_GENE_FAMILIES = std::stoi(parse_parameter_value);
			}
			if (parse_parameters == "MULTINOMIAL_ELEMENT_OUTPUT_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::MULTINOMIAL_ELEMENT_OUTPUT_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "REPORT_QUERY_ALIGNMENT_TOTALSCORE_OUTPUT_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::REPORT_QUERY_ALIGNMENT_TOTALSCORE_OUTPUT_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "PROTEIN_CONSTRUCT_PARAMETER_DENSITY_CONJUGATED_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::PROTEIN_CONSTRUCT_PARAMETER_DENSITY_CONJUGATED_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "PROTEIN_CONSTRUCT_SCORE_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::PROTEIN_CONSTRUCT_SCORE_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "REPORT_QUERY_PARAMETER_DENSITY_CONJUGATED_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::REPORT_QUERY_PARAMETER_DENSITY_CONJUGATED_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "REPORT_PROTEIN_SCORE_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::REPORT_PROTEIN_SCORE_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "REPORT_PROTEIN_DENSITY_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::REPORT_PROTEIN_DENSITY_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "REPORT_QUERY_SCORE_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::REPORT_QUERY_SCORE_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "REPORT_QUERY_PARAMETER_SCORE_DENSITY_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::REPORT_QUERY_PARAMETER_SCORE_DENSITY_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "REPORT_V_HOMOLOGY_DATA_AGGREGATED_BY_PROTEIN_CONJUGATED_DENSITY_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::REPORT_V_HOMOLOGY_DATA_AGGREGATED_BY_PROTEIN_CONJUGATED_DENSITY_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "DEFAULT_PROTEIN_SCORE_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::DEFAULT_PROTEIN_SCORE_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "HOMOLOGY_DENSITY_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::HOMOLOGY_DENSITY_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "HOMOLOGY_DENSITY_CONJUGATED_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != ',')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::HOMOLOGY_DENSITY_CONJUGATED_THRESHOLD = std::stod(parse_parameter_value);
			}
			if (parse_parameters == "CLUSTER_PROPORTION_THRESHOLD")
			{
				string parse_parameter_value{};
				while (read_parameters != '\n')
				{
					fin_input_parameters.get(read_parameters);
					if ((read_parameters != ':') && (read_parameters != ' ') && (read_parameters != ','))
					{
						parse_parameter_value += read_parameters;
					}
				}
				IgFamily::CLUSTER_PROPORTION_THRESHOLD = std::stod(parse_parameter_value);
			}
		}
	}
}
#endif