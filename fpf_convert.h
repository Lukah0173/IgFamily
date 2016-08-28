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


namespace fpf_convert {

	using std::string;
	using std::vector;

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
}

#endif