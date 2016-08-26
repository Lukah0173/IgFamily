// * * fpf_convert.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_CONVERT
#define	FPF_CONVERT
#include <cstdlib>						// provides - size_t
#include <string>						// provides - std::string
#include <vector>						// provides - std::vector



namespace fpf_convert {

	using std::string;
	using std::vector;

	struct threshold_type;
	struct peakpicking_type;
	struct ms2deisotope_type;
	struct ms2denoise_type;
	struct chargestatepredictor_type;
	struct fileconversion_parameters;

	struct threshold_type {
		bool intensity_threshold;
		size_t threshold;
	};

	struct peakpicking_type {
		bool peakpicking;
		size_t peakpicking_mslevel_from;
		size_t peakpicking_mslevel_to;
	};

	struct ms2deisotope_type {
		bool ms2deisotope;
	};

	struct ms2denoise_type {
		bool ms2denoise;
		size_t ms2denoise_windowwidth;
		size_t ms2denoise_peaksinwindow;
	};

	struct chargestatepredictor_type {
	public:
		bool chargestatepredictor;
		bool chargestateoverride;
		size_t chargestatepredictor_maxcharge;
		size_t chargestatepredictor_mincharge;
		double chargestatepredictor_chargefraction;
	};

	struct fileconversion_parameters {
		threshold_type threshold;
		peakpicking_type peakpicking;
		ms2deisotope_type ms2deisotope;
		ms2denoise_type ms2denoise;
		chargestatepredictor_type chargestatepredictor;
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

	fileconversion_parameters create_fileconversion_parameters(bool par_defaultconversion) {
		fileconversion_parameters fileconversion = fileconversion_parameters();
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
		string string_system = "CD Z:\\Lukah_Dykes\\IgFamily\\proteowizard\\";
		string_system += " && ";
		string_system += par_str_msconvert_command;
		system(string_system.c_str());
	}
}

#endif