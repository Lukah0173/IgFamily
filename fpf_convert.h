// * * fpf_convert.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_CONVERT
#define	FPF_CONVERT

#include <cstdlib>	
#include <iostream>
#include <string>
#include <vector>

#include "IgFamily.h"
#include "fpf_filesystem.h"


namespace fpf_convert {

	using std::string;
	using std::vector;

	typedef fpf_filesystem::filesystem filesystem;

	struct conversion_absolute_threshold;
	struct conversion_peakpicking;
	struct conversion_ms2deisotope;
	struct conversion_ms2denoise;
	struct conversion_chargestatepredictor;
	struct file_conversion;

	struct conversion_absolute_threshold {
		bool intensity_threshold;
		size_t threshold;
	};

	struct conversion_peakpicking {
		bool peakpicking;
		size_t peakpicking_mslevel_from;
		size_t peakpicking_mslevel_to;
	};

	struct conversion_ms2deisotope {
		bool ms2deisotope;
	};

	struct conversion_ms2denoise {
		bool ms2denoise;
		size_t ms2denoise_windowwidth;
		size_t ms2denoise_peaksinwindow;
	};

	struct conversion_chargestatepredictor {
	public:
		bool chargestatepredictor;
		bool chargestateoverride;
		size_t chargestatepredictor_maxcharge;
		size_t chargestatepredictor_mincharge;
		double chargestatepredictor_chargefraction;
	};

	struct file_conversion {
		conversion_absolute_threshold threshold;
		conversion_peakpicking peakpicking;
		conversion_ms2deisotope ms2deisotope;
		conversion_ms2denoise ms2denoise;
		conversion_chargestatepredictor chargestatepredictor;
	};

	bool prompt_defaultconversion() {
		std::cout << "\n\n\n\n default .wiff conversion parameters - ";
		std::cout << "\n\n   ~ absolute intensity threshold";
		std::cout << "\n\n     - intensity: 3";
		std::cout << "\n\n   ~ peak picking";
		std::cout << "\n\n     - MS levels: 1-2";
		std::cout << "\n\n   ~ deisotope";
		std::cout << "\n\n   ~ denoise ms2";
		std::cout << "\n\n     - window width: 30 Da";
		std::cout << "\n\n     - most intense peaks in window: 6";
		std::cout << "\n\n   ~ charge state predictor";
		std::cout << "\n\n     - minimum charge: 2";
		std::cout << "\n\n     - maximum charge: 3";
		std::cout << "\n\n     - charge fraction: 0.9";
		std::cout << "\n\n\n use default .wiff conversion parameters?";
		std::cout << "\n\n (y = yes, any other value = no)";
		
		string str_prompt;
		while (str_prompt != "y") {
			std::cout << "\n\n -> ";
			std::cin >> str_prompt;
		}
		if (str_prompt == "y") {
			return true;
		}
		return true;
	}

	file_conversion create_fileconversion_parameters(bool par_defaultconversion) {
		file_conversion fileconversion = file_conversion();
		if (par_defaultconversion) {
			fileconversion.threshold.intensity_threshold = true;
			fileconversion.threshold.threshold = 3;
			fileconversion.peakpicking.peakpicking = true;
			fileconversion.peakpicking.peakpicking_mslevel_from = 1;
			fileconversion.peakpicking.peakpicking_mslevel_to = 2;
			fileconversion.ms2deisotope.ms2deisotope = true;
			fileconversion.ms2denoise.ms2denoise = true;
			fileconversion.ms2denoise.ms2denoise_windowwidth = 30;
			fileconversion.ms2denoise.ms2denoise_peaksinwindow = 6;
			fileconversion.chargestatepredictor.chargestatepredictor = true;
			fileconversion.chargestatepredictor.chargestatepredictor = true;
			fileconversion.chargestatepredictor.chargestatepredictor_maxcharge = 3;
			fileconversion.chargestatepredictor.chargestatepredictor_mincharge = 2;
			fileconversion.chargestatepredictor.chargestatepredictor_chargefraction = 0.9;
		}
		return fileconversion;
	}

	void sys_msconvert(string par_str_msconvert_command, string par_str_filesystem_directory) {
		std::cout << "\n\n";
		string string_system = "CD " + IgFamily::DEFAULT_MSCONVERT_DIRECTORY;
		string_system += " && ";
		string_system += par_str_msconvert_command;
		system(string_system.c_str());
	}

	void perform_fileconversion(filesystem& par_filesystem) {
		file_conversion temp_file_conversion = fpf_convert::create_fileconversion_parameters(fpf_convert::prompt_defaultconversion());
		string fileconversion_command{};
		fileconversion_command += "msconvert.exe ";
		fileconversion_command += "\"" + DEFAULT_IGFAMILY_DIRECTORY;
		fileconversion_command += par_filesystem.directory;
		fileconversion_command += par_filesystem.filename;
		fileconversion_command += ".wiff\"";
		////fileconversion_command += " --64";
		////fileconversion_command += " --mz64";
		//fileconversion_command += " -v";
		//fileconversion_command += " --mgf";
		//if (temp_file_conversion.peakpicking.peakpicking) {
		//	//fileconversion_command += " --filter \"peakPicking cwt ";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion.peakpicking.peakpicking_mslevel_from);
		//	//fileconversion_command += "-";
		//	//fileconversion_command += std::to_string(par_filesystem.fileconversion.peakpicking.peakpicking_mslevel_to);
		//	//fileconversion_command += "\"";
		//}
		//if (temp_file_conversion.threshold.intensity_threshold) {
		//	//fileconversion_command += " --filter \"threshold absolute ";
		//	//fileconversion_command += std::to_string(par_filesystem.file_conversion.threshold.threshold);
		//	//fileconversion_command += " most-intense\"";
		//}
		//if (temp_file_conversion.ms2denoise.ms2denoise) {
		//	//fileconversion_command += " --filter \"MS2Denoise ";
		//	//fileconversion_command += std::to_string(par_filesystem.file_conversion.ms2denoise.ms2denoise_peaksinwindow);
		//	//fileconversion_command += " ";
		//	//fileconversion_command += std::to_string(par_filesystem.file_conversion.ms2denoise.ms2denoise_windowwidth);
		//	//fileconversion_command += " true\"";
		//}
		//if (temp_file_conversion.ms2deisotope.ms2deisotope) {
		//	//fileconversion_command += " --filter MS2Deisotope";
		//}
		//if (temp_file_conversion.chargestatepredictor.chargestatepredictor) {
		//	fileconversion_command += " --filter \"chargeStatePredictor true ";
		//	fileconversion_command += std::to_string(par_filesystem.file_conversion.chargestatepredictor.chargestatepredictor_maxcharge);
		//	fileconversion_command += " ";
		//	fileconversion_command += std::to_string(par_filesystem.file_conversion.chargestatepredictor.chargestatepredictor_mincharge);
		//	fileconversion_command += " ";
		//	fileconversion_command += "0.9";
		//	fileconversion_command += "\"";
		//}
		//fileconversion_command += " outdir=Z:\\Lukah_Dykes\\IgFamily\\";
		//fileconversion_command += par_filesystem.directory;
		fileconversion_command += " --mgf --filter \"chargeStatePredictor true 3 2 0.9\"";
		std::cout << "\n\n" << fileconversion_command;
		sys_msconvert(fileconversion_command, par_filesystem.directory);
	}
}

#endif