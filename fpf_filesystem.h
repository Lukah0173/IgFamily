// * * fpf_filesystem.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_FILESYSTEM
#define	FPF_FILESYSTEM

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <map>
#include <utility>
#include <vector>

#include "IgFamily.h"
#include "fpf_data.h"
#include "fpf_parse.h"


namespace fpf_filesystem {

	using std::pair;
	using std::string;
	using std::map;
	using std::multimap;
	using std::vector;

	typedef fpf_data::homology_data homology_data;
	typedef fpf_data::protein_data protein_data;
	typedef fpf_data::peptide_analysis peptide_analysis;
	typedef fpf_data::peptide_data peptide_data;
	typedef fpf_data::protein_analysis protein_analysis;
	typedef fpf_data::multinomial multinomial;
	typedef fpf_parse::csv_data csv_data;
	typedef fpf_parse::FASTA_data FASTA_data;

	struct filesystem;
	struct sample_analysis;
	
	struct sample_analysis {
		string peptide_assignment_method;
		bool file_found;
		vector<csv_data> main_v_csv_peptides;
		vector<protein_data> v_protein_data;
		map<string, protein_data*> v_protein_data_map;
		vector<peptide_data> v_peptide_data;
		multimap<string, peptide_data*> v_peptide_data_map;
		vector<peptide_analysis> v_peptide_analysis;
		map<string, peptide_analysis*> v_peptide_analysis_map;
		vector<homology_data> v_homology_data;
		vector<protein_analysis> v_protein_analysis;
		map<string, protein_analysis*> v_protein_analysis_map;
		vector<protein_analysis> v_protein_analysis_selected_by_polymorphism;
		multinomial multinomial_data;
		double protein_analysis_score_mean;
	};

	struct filesystem {
	public:

		/* filesystem */

		pair<string, string> filesystem_id;
		string directory;
		string filename;
		string fileversion;
		vector<pair<string, string>> v_filesystem_replicates;
		size_t filesystem_replicate_count;

		/* sample factors */

		string patientstatus;
		string enzyme;	

		/* mgf conversion parameters */

		bool perform_wiff_fileconversion;

		/* de novo parameters */

		bool perform_novor_denovo;
		string denono_deltamass;

		/* sample analysis */

		std::vector<sample_analysis> v_sample_analysis;
	};

	vector<string> read_root_dir(string par_root_directory) {
		std::cout << "\n\n\n reading root directory...\n";
		vector<string> temp_v_IgFamily_root{};
		string fin_IgFamily_root{};
		if (IgFamily::FILESYSTEM_MODE) {
			std::ifstream fin_input_csv(par_root_directory);
			char stream_IgFamily_root{};
			while (fin_input_csv.std::istream::get(stream_IgFamily_root)) {
				if ((stream_IgFamily_root != '\n') && (stream_IgFamily_root != ',')) {
					fin_IgFamily_root += stream_IgFamily_root;
				}
				if (stream_IgFamily_root == ',') {
					temp_v_IgFamily_root.push_back(fin_IgFamily_root);
					std::cout << "\n * " << fin_IgFamily_root;
					fin_IgFamily_root.clear();
				}
			}
		}
		return temp_v_IgFamily_root;
	}

	vector<filesystem> read_filesystem(vector<string> par_root_directory) {
		filesystem temp_filesystem{};
		vector<filesystem> temp_v_filesystem{};
		if (IgFamily::FILESYSTEM_MODE) {
			std::cout << "\n";
			for (const auto& itr_root_directory : par_root_directory) {
				string stream_filesystem = itr_root_directory + "filesystem.data";
				std::ifstream fin_input_filesystem(stream_filesystem);
				char read_fin_filesystem{};
				string stream_fin_filesystem{};
				size_t switch_fin_filesystem{};
				pair<string, string> temp_filesystem_id{};
				string temp_filesystem_date{};
				string temp_filesystem_filename{};
				string temp_filesystem_version{};
				string temp_filesystem_status{};
				string temp_filesystem_enzyme{};
				string temp_filesystem_denovo_deltamass{};
				string temp_filesystem_replicatedate{};
				size_t switch__filesystem_replicatepair{};
				vector<pair<string, string>> v_filesystem_replicates{};
				while (fin_input_filesystem.std::istream::get(read_fin_filesystem)) {
					if ((read_fin_filesystem != ';') && (read_fin_filesystem != '\n') && (read_fin_filesystem != ',')) {
						stream_fin_filesystem += read_fin_filesystem;
					}
					if (switch__filesystem_replicatepair == 2) {
						switch__filesystem_replicatepair = 0;
					}
					if ((switch_fin_filesystem == 1) && (read_fin_filesystem == ',')) {
						temp_filesystem_date = stream_fin_filesystem;
						stream_fin_filesystem.clear();
					}
					if (stream_fin_filesystem == "ID: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 1;
					}
					if (stream_fin_filesystem == "FILE: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 2;
					}
					if (stream_fin_filesystem == "VERSION: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 3;
					}
					if ((switch_fin_filesystem == 4) && (read_fin_filesystem == ',') && (switch__filesystem_replicatepair == 1)) {
						v_filesystem_replicates.push_back(std::make_pair(temp_filesystem_replicatedate, stream_fin_filesystem));
						switch__filesystem_replicatepair = 2;
						temp_filesystem_replicatedate.clear();
						stream_fin_filesystem.clear();
					}
					if ((switch_fin_filesystem == 4) && (read_fin_filesystem == ';') && (switch__filesystem_replicatepair == 0)) {
						temp_filesystem_replicatedate = stream_fin_filesystem;
						switch__filesystem_replicatepair = 1;
						stream_fin_filesystem.clear();
					}
					if (stream_fin_filesystem == "REPLICATES: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 4;
					}
					if (stream_fin_filesystem == "STATUS: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 5;
					}
					if (stream_fin_filesystem == "ENZYME: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 6;
					}
					if (stream_fin_filesystem == "DENOVO_DELTAMASS: ") {
						stream_fin_filesystem.clear();
						switch_fin_filesystem = 7;
					}
					if (read_fin_filesystem == ';') {
						if (switch_fin_filesystem == 1) {
							temp_filesystem_id = make_pair(temp_filesystem_date, stream_fin_filesystem);
							temp_filesystem_date.clear();
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 2) {
							temp_filesystem_filename = stream_fin_filesystem;
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 3) {
							temp_filesystem_version = stream_fin_filesystem;
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 4) {
							if (stream_fin_filesystem != "") {
								v_filesystem_replicates.push_back(std::make_pair(temp_filesystem_replicatedate, stream_fin_filesystem));
							}
							switch__filesystem_replicatepair = size_t();
							temp_filesystem_replicatedate.clear();
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 5) {
							temp_filesystem_status = stream_fin_filesystem;
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 6) {
							temp_filesystem_enzyme = stream_fin_filesystem;
							stream_fin_filesystem.clear();
						}
						if (switch_fin_filesystem == 7) {
							temp_filesystem_denovo_deltamass = stream_fin_filesystem;
							stream_fin_filesystem.clear();
							temp_filesystem.directory = itr_root_directory;
							temp_filesystem.filesystem_id = temp_filesystem_id;
							temp_filesystem.filename = temp_filesystem_filename;
							temp_filesystem.fileversion = temp_filesystem_version;
							temp_filesystem.patientstatus = temp_filesystem_status;
							temp_filesystem.v_filesystem_replicates = v_filesystem_replicates;
							temp_filesystem.enzyme = temp_filesystem_enzyme;
							temp_filesystem.denono_deltamass = temp_filesystem_denovo_deltamass;
							temp_filesystem.filesystem_replicate_count = size_t{ 1 };
							v_filesystem_replicates.clear();
							temp_v_filesystem.push_back(temp_filesystem);
							std::cout << "\n * reading: " << temp_filesystem_filename << "   version - " << temp_filesystem_version;
							if (temp_filesystem_version == IgFamily::version) {
								std::cout << "   ! up to date !";
							}
						}
					}
				}
				fin_input_filesystem.clear();
				fin_input_filesystem.seekg(0, std::ios::beg);
			}
		}
		else {
			std::cout << "\n * reading: Local file";
			temp_filesystem.filename = "local_file";
			temp_filesystem.fileversion = IgFamily::version;
			temp_v_filesystem.push_back(temp_filesystem);
		}
		return temp_v_filesystem;
	}

	string read_filesystem_PEAKS_database_peptides(string par_root_directory) {
		string temp_root_PEAKS_database_peptides = par_root_directory + "database_peptides_PEAKS.csv";
		return temp_root_PEAKS_database_peptides;
	}

	string read_filesystem_PEAKS_denovo_peptides(string par_root_directory) {
		string temp_root_PEAKS_denovo_peptides{};
		temp_root_PEAKS_denovo_peptides = par_root_directory + "denovo_peptides_PEAKS.csv";
		return temp_root_PEAKS_denovo_peptides;
	}

	string read_filesystem_NOVOR_denovo_peptides(string par_root_directory) {
		string temp_root_NOVOR_database_peptides{};
		temp_root_NOVOR_database_peptides = par_root_directory + "denovo_peptides_NOVOR.csv";
		return temp_root_NOVOR_database_peptides;
	}

	vector<csv_data> parse_filesystem_PEAKS_database_peptides(string par_fin_root_directory) {
		vector<csv_data> temp_v_csv_proteinpeptides{};
		std::ifstream fin_input_csv(par_fin_root_directory);
		temp_v_csv_proteinpeptides = fpf_parse::parse_proteinpeptides(fin_input_csv, par_fin_root_directory);
		return temp_v_csv_proteinpeptides;
	}

	vector<csv_data> parse_filesystem_PEAKS_denovo_peptides(string par_fin_root_directory) {
		vector<csv_data> temp_csv_PEAKS_denovo_peptides{};
		std::ifstream fin_input_csv(par_fin_root_directory);
		temp_csv_PEAKS_denovo_peptides = fpf_parse::parse_csv_PEAKS_denovopeptides(fin_input_csv, par_fin_root_directory);
		return temp_csv_PEAKS_denovo_peptides;
	}

	vector<csv_data> parse_filesystem_NOVOR_denovo_peptides(string par_fin_root_directory) {
		vector<csv_data> temp_csv_NOVOR_denovo_peptides{};
		std::ifstream fin_input_csv(par_fin_root_directory);
		temp_csv_NOVOR_denovo_peptides = fpf_parse::parse_csv_NOVOR_denovopeptides(fin_input_csv, par_fin_root_directory);
		return temp_csv_NOVOR_denovo_peptides;
	}

	void fout_filesystem(filesystem& par_filesystem) {
		if (IgFamily::FILESYSTEM_MODE) {
			string output_filesystem = par_filesystem.directory + "filesystem.data";
			std::ofstream fout_filesystem;
			fout_filesystem.open(output_filesystem);
			fout_filesystem << "ID: " << std::get<0>(par_filesystem.filesystem_id) << "," << std::get<1>(par_filesystem.filesystem_id) << ";\n";
			fout_filesystem << "FILE: " << par_filesystem.filename << ";\n";
			fout_filesystem << "VERSION: " << IgFamily::version << ";\n";
			par_filesystem.fileversion = IgFamily::version;
			fout_filesystem << "REPLICATES: ";
			for (const auto& itr_v_p_replicates : par_filesystem.v_filesystem_replicates) {
				fout_filesystem << std::get<0>(itr_v_p_replicates) << ",";
				fout_filesystem << std::get<1>(itr_v_p_replicates) << ";";
			}
			fout_filesystem << "\n";
			fout_filesystem << "STATUS: " << par_filesystem.patientstatus << ";\n";
			fout_filesystem << "ENZYME: " << par_filesystem.enzyme << ";\n";
			fout_filesystem << "DENOVO_DELTAMASS: " << par_filesystem.denono_deltamass << ";\n";
		}
	}
}

#endif