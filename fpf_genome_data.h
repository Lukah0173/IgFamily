// * * fpf_genome_data.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_GENOME_DATA
#define	FPF_GENOME_DATA

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "IgFamily.h"


namespace fpf_genome_data {

	using std::string;
	using std::vector;

	struct genome_data {
	public:
		string genome_ID;
		string genome_Vregion_genefamily;
		string genome_translation_sequence;
	};

	struct genome_analysis {
	public:
		genome_data* p_genome_data;
		size_t count_genome_data_replicate;
	};

	vector<genome_data> create_v_genome_data() {
		string genomic_data_input = IgFamily::DEFAULT_GENOME_DIRECTORY + "Sjogrens_syndrome_patient_data\\SS1_AA.csv";
		std::ifstream fin_input_genomic_data(genomic_data_input);
		char genomic_data_read{};
		vector<genome_data> temp_v_genomic_data{};
		genome_data temp_genomic_data{};
		string temp_parse_genomic_data{};
		size_t genome_data_count_delimit{};
		size_t genome_data_count_delimit_width{};
		bool header_parsed{};
		bool skip_this{};
		while (fin_input_genomic_data.get(genomic_data_read)) {
			if (header_parsed) {
				if (genomic_data_read == ',') {
					if (genome_data_count_delimit % genome_data_count_delimit_width == 1) {
						temp_genomic_data.genome_ID = temp_parse_genomic_data;
					}
					if (genome_data_count_delimit % genome_data_count_delimit_width == 3) {
						if (temp_parse_genomic_data.find("or") != std::string::npos) {
							skip_this = true;
						}
						else {
							temp_genomic_data.genome_Vregion_genefamily.clear();
							bool read_temp_parse_genomic_data{};
							for (const auto& itr_temp_parse_genomic_data : temp_parse_genomic_data) {
								if (itr_temp_parse_genomic_data == ' ') {
									if (read_temp_parse_genomic_data) {
										read_temp_parse_genomic_data = false;
									}
									else {
										read_temp_parse_genomic_data = true;
									}
								}
								if (read_temp_parse_genomic_data) {
									temp_genomic_data.genome_Vregion_genefamily += itr_temp_parse_genomic_data;
								}
							}
						}
					}
					if (genome_data_count_delimit % genome_data_count_delimit_width == 6) {
						temp_genomic_data.genome_translation_sequence = temp_parse_genomic_data;
					}
					temp_parse_genomic_data.clear();
					++genome_data_count_delimit;
				}
				if ((genomic_data_read != ',') && (genomic_data_read != '\n')) {
					temp_parse_genomic_data += genomic_data_read;
				}
				if (genomic_data_read == '\n') {
					if (!skip_this) {
						temp_v_genomic_data.push_back(temp_genomic_data);
					}
					temp_parse_genomic_data.clear();
					skip_this = false;
				}
			}
			if (!header_parsed) {
				if (genomic_data_read == ',') {
					++genome_data_count_delimit;
				}
				if (genomic_data_read == '\n') {
					header_parsed = true;
					genome_data_count_delimit_width = genome_data_count_delimit;
				}
			}
		}
		return temp_v_genomic_data;
	}

	vector<genome_analysis> create_v_genome_analysis(vector<genome_data>& par_v_genome_data) {
		vector<genome_analysis> temp_v_genome_analysis{};
		for (auto& itr_v_genome_data : par_v_genome_data) {	
			genome_analysis temp_genome_analysis{};
			auto& find_genome_analysis = std::find_if(temp_v_genome_analysis.begin(), temp_v_genome_analysis.end(),
				[itr_v_genome_data](const genome_analysis& par_genome_analysis) {
				return par_genome_analysis.p_genome_data->genome_translation_sequence == itr_v_genome_data.genome_translation_sequence;
			});
			if (find_genome_analysis == temp_v_genome_analysis.end()) {
				temp_genome_analysis.p_genome_data = &itr_v_genome_data;
				++temp_genome_analysis.count_genome_data_replicate;
				temp_v_genome_analysis.push_back(temp_genome_analysis);
			}
			else {
				++find_genome_analysis->count_genome_data_replicate;
			}
		}
		return temp_v_genome_analysis;
	}

	void fout_v_genome_data(const vector<genome_data>& par_v_genome_data) {
		std::string output_v_genome_data = IgFamily::DEFAULT_GENOME_DIRECTORY + "genome_data.csv";
		std::ofstream fout_v_genome_data;
		fout_v_genome_data.open(output_v_genome_data);
		fout_v_genome_data << "genome_ID,genome_Vregion_genefamily,genome_translation_sequence,\n";
		for (const auto& itr_v_genome_data : par_v_genome_data) {
			fout_v_genome_data << itr_v_genome_data.genome_ID << ",";
			fout_v_genome_data << itr_v_genome_data.genome_Vregion_genefamily << ",";
			fout_v_genome_data << itr_v_genome_data.genome_translation_sequence << ",";
			fout_v_genome_data << ",\n";
		}
	}

	void fout_v_genome_analysis(const vector<genome_analysis>& par_v_genome_analysis) {
		std::string output_v_genome_analysis = IgFamily::DEFAULT_GENOME_DIRECTORY + "genome_analysis.csv";
		std::ofstream fout_v_genome_analysis;
		fout_v_genome_analysis.open(output_v_genome_analysis);
		fout_v_genome_analysis << "count_genome_data_replicate,genome_Vregion_genefamily,genome_translation_sequence,\n";
		for (const auto& itr_v_genome_analysis : par_v_genome_analysis) {
			fout_v_genome_analysis << itr_v_genome_analysis.count_genome_data_replicate << ",";
			fout_v_genome_analysis << itr_v_genome_analysis.p_genome_data->genome_Vregion_genefamily << ",";
			fout_v_genome_analysis << itr_v_genome_analysis.p_genome_data->genome_translation_sequence << ",";
			fout_v_genome_analysis << ",\n";
		}
	}
}

#endif
