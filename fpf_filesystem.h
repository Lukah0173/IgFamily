// * * fpf_filesystem.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_FILESYSTEM
#define	FPF_FILESYSTEM

#include <cstdlib> // provides - size_t
#include <vector> // provides - std::vector
#include <utility> // provides - std::pair
#include <iostream> // provides - std::istream
#include <algorithm> // provides - std::find
#include <tuple> // provides - std::tuple

#include "IgFamily.h"
#include "fpf_parse.h"
#include "fpf_data.h"
#include "fpf_convert.h"


namespace fpf_filesystem {

	using std::string;
	using std::vector;
	using std::pair;

	typedef fpf_convert::fileconversion_parameters fileconversion_parameters;
	typedef fpf_data::peptide_data peptide_data;
	typedef fpf_data::peptide_analysis peptide_analysis;
	typedef fpf_data::FASTA_category FASTA_category;
	typedef fpf_data::blastp_data blastp_data;
	typedef fpf_data::multinomial multinomial;
	typedef fpf_data::category_analysis category_analysis;
	typedef fpf_parse::csv_data csv_data;
	typedef fpf_parse::FASTA_data FASTA_data;

	struct filesystem;

	struct filesystem {
	public:
		bool proteinpeptides_exist;
		bool denovopeptides_exist;

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

		bool fileconversion = true;
		fileconversion_parameters fileconversion_parameters;

		/* de novo parameters */

		string denono_deltamass;

		/* sample data */

		vector<FASTA_category> v_FASTA_category;
		vector<peptide_data> v_peptide_data;
		vector<peptide_analysis> v_peptide_analysis;

		/* blastp data */

		vector<blastp_data> v_blastp_data;

		/* category analysis */

		vector<category_analysis> v_category_analysis;
		vector<category_analysis> v_category_analysis_selected_by_polymorphism;

		/* multinomial data */
		
		multinomial multinomial_data;
	};

	void display_default() {
		std::cout << " Current settings:";
		std::cout << "\n\n\n";
		std::cout << " FASTA file -           ";
		std::cout << IgFamily::DEFAULT_INPUT_FASTA;
		std::cout << "\n\n\n";
		std::cout << " Workflow ~ ";
		std::cout << "\n\n\n";
		std::cout << "  [FILE].wiff + [FILE].wiff.scan";
		std::cout << "\n\n";
		std::cout << "     |\n";
		std::cout << "     V";
		std::cout << "\n\n";
		std::cout << "  [FILE].mgf ";
		std::cout << "\n\n";
		std::cout << "     |\n";
		std::cout << "     V";
		std::cout << "\n\n";
		std::cout << "  [FILE]_denovopeptides.csv    ";
		std::cout << "\n\n";
		std::cout << "     |\n";
		std::cout << "     V";
		std::cout << "\n\n";
		std::cout << "  [FILE]_report.html        \n";
		std::cout << "  [FILE]_report_IG.html      ";
	}

	string display_menu() {
		std::cout << "\n\n";
		std::cout << " FASTA utilities:     [F] ";
		std::cout << "\n Continue:            [X]";
		string menu_selection{};
		while ((menu_selection != ("F")) && (menu_selection != ("X"))) {
			menu_selection.clear();
			std::cout << "\n\n Input selection: \n\n > ";
			std::cin >> menu_selection;
		}
		return menu_selection;
	}

	int perform_menu_selection(string par_menu_selection) {
		if (par_menu_selection == "F") {
			std::cout << "\n Read FASTA format:              [R] ";
			std::cout << "\n Output custom FASTA format:     [C] ";
			string menu_selection{};
			while ((menu_selection != ("R")) && (menu_selection != ("C"))) {
				menu_selection.clear();
				std::cout << "\n\n Input selection: \n\n > ";
				std::cin >> menu_selection;
			}
			if (menu_selection == "R") {
				fpf_parse::check_FASTA_format(IgFamily::DEFAULT_INPUT_FASTA);
			}
			if (menu_selection == "C") {
				fpf_parse::output_custom_FASTA_format(IgFamily::DEFAULT_INPUT_FASTA);
			}
			return 0;
		}
		if (par_menu_selection == "X") {
			return 0;
		}
		return 1;
	}

	vector<string> read_root_dir(string par_root_directory) {
		std::ifstream fin_input_csv(par_root_directory);
		vector<string> temp_v_IgFamily_root{};
		string fin_IgFamily_root{};
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
		return temp_v_IgFamily_root;
	}

	vector<filesystem> read_filesystem(vector<string> par_root_directory) {
		filesystem temp_filesystem{};
		vector<filesystem> temp_v_filesystem{};
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
		return temp_v_filesystem;
	}

	void perform_fileconversion(filesystem& par_filesystem) {
		par_filesystem.fileconversion_parameters = fpf_convert::create_fileconversion_parameters(fpf_convert::prompt_defaultconversion());
		string fileconversion_command{};
		fileconversion_command += "msconvert.exe ";
		fileconversion_command += "\"Z:\\Lukah_Dykes\\IgFamily\\";
		fileconversion_command += par_filesystem.directory;
		fileconversion_command += par_filesystem.filename;
		fileconversion_command += ".wiff\"";
		////fileconversion_command += " --64";
		////fileconversion_command += " --mz64";
		//fileconversion_command += " -v";
		//fileconversion_command += " --mgf";
		//if (par_filesystem.fileconversion_parameters.peakpicking.peakpicking) {
		//	//fileconversion_command += " --filter \"peakPicking cwt ";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion.peakpicking.peakpicking_mslevel_from);
		//	//fileconversion_command += "-";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion.peakpicking.peakpicking_mslevel_to);
		//	//fileconversion_command += "\"";
		//}
		//if (par_filesystem.fileconversion_parameters.threshold.intensity_threshold) {
		//	//fileconversion_command += " --filter \"threshold absolute ";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.threshold.threshold);
		//	//fileconversion_command += " most-intense\"";
		//}
		//if (par_filesystem.fileconversion_parameters.ms2denoise.ms2denoise) {
		//	//fileconversion_command += " --filter \"MS2Denoise ";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.ms2denoise.ms2denoise_peaksinwindow);
		//	//fileconversion_command += " ";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.ms2denoise.ms2denoise_windowwidth);
		//	//fileconversion_command += " true\"";
		//}
		//if (par_filesystem.fileconversion_parameters.ms2deisotope.ms2deisotope) {
		//	//fileconversion_command += " --filter MS2Deisotope";
		//}
		//if (par_filesystem.fileconversion_parameters.chargestatepredictor.chargestatepredictor) {
		//	fileconversion_command += " --filter \"chargeStatePredictor true ";
		//	fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.chargestatepredictor.chargestatepredictor_maxcharge);
		//	fileconversion_command += " ";
		//	fileconversion_command += std::to_string(par_filesystem.fileconversion_parameters.chargestatepredictor.chargestatepredictor_mincharge);
		//	fileconversion_command += " ";
		//	fileconversion_command += "0.9";
		//	fileconversion_command += "\"";
		//}
		//fileconversion_command += " outdir=Z:\\Lukah_Dykes\\IgFamily\\";
		//fileconversion_command += par_filesystem.directory;
		fileconversion_command += " --mgf --filter \"chargeStatePredictor true 3 2 0.9\"";
		std::cout << "\n\n" << fileconversion_command;
		fpf_convert::sys_msconvert(fileconversion_command, par_filesystem.directory);
	}

	string read_filesystem_proteinpeptides(string par_root_directory) {
		string temp_root_proteinpeptides = par_root_directory + "protein_peptides.csv";
		return temp_root_proteinpeptides;
	}

	string read_filesystem_denovopeptides(string par_root_directory) {
		string temp_root_denovopeptides{};
		if (IgFamily::NOVOR_DENOVO) {
			temp_root_denovopeptides = par_root_directory + "denovo_peptides_NOVOR.csv";
		}
		else {
			temp_root_denovopeptides = par_root_directory + "denovo_peptides.csv";
		}
		return temp_root_denovopeptides;
	}

	vector<csv_data> parse_filesystem_proteinpeptides(string par_fin_root_directory) {
		vector<csv_data> temp_v_csv_proteinpeptides{};
		std::ifstream fin_input_csv(par_fin_root_directory);
		temp_v_csv_proteinpeptides = fpf_parse::parse_proteinpeptides(fin_input_csv, par_fin_root_directory);
		return temp_v_csv_proteinpeptides;
	}

	vector<csv_data> parse_filesystem_denovopeptides(string par_fin_root_directory) {
		vector<csv_data> temp_csv_denovopeptides{};
		std::ifstream fin_input_csv(par_fin_root_directory);
		if (IgFamily::NOVOR_DENOVO) {
			temp_csv_denovopeptides = fpf_parse::parse_csv_NOVOR_denovopeptides(fin_input_csv, par_fin_root_directory);
		}
		else {
			temp_csv_denovopeptides = fpf_parse::parse_csv_PEAKS_denovopeptides(fin_input_csv, par_fin_root_directory);
		}
		return temp_csv_denovopeptides;
	}

	void fout_filesystem(filesystem par_filesystem) {
		string output_filesystem = par_filesystem.directory + "filesystem.data";
		std::ofstream fout_filesystem;
		fout_filesystem.open(output_filesystem);
		fout_filesystem << "ID: " << std::get<0>(par_filesystem.filesystem_id) << "," << std::get<1>(par_filesystem.filesystem_id) << ";\n";
		fout_filesystem << "FILE: " << par_filesystem.filename << ";\n";
		fout_filesystem << "VERSION: " << IgFamily::version << ";\n";
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

#endif