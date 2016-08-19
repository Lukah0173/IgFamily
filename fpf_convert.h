// * * fpf_convert.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_CONVERT
#define	FPF_CONVERT
#include <cstdlib>						// provides - size_t
#include <string>						// provides - string
#include <vector>						// provides - vector



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
		bool b_threshold;
		size_t st_threshold;
	};

	struct peakpicking_type {
		bool b_peakpicking;
		size_t st_peakpicking_mslevel_from;
		size_t st_peakpicking_mslevel_to;
	};

	struct ms2deisotope_type {
		bool b_ms2deisotope;
	};

	struct ms2denoise_type {
		bool b_ms2denoise;
		size_t st_ms2denoise_windowwidth;
		size_t st_ms2denoise_peaksinwindow;
	};

	struct chargestatepredictor_type {
	public:
		bool b_chargestatepredictor;
		bool b_chargestateoverride;
		size_t st_chargestatepredictor_maxcharge;
		size_t st_chargestatepredictor_mincharge;
		double d_chargestatepredictor_chargefraction;
	};

	struct fileconversion_parameters {
		threshold_type s_threshold;
		peakpicking_type s_peakpicking;
		ms2deisotope_type s_ms2deisotope;
		ms2denoise_type s_ms2denoise;
		chargestatepredictor_type s_chargestatepredictor;
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
			fileconversion.s_threshold.b_threshold = true;
			fileconversion.s_threshold.st_threshold = 3;
			fileconversion.s_peakpicking.b_peakpicking = true;
			fileconversion.s_peakpicking.st_peakpicking_mslevel_from = 1;
			fileconversion.s_peakpicking.st_peakpicking_mslevel_to = 2;
			fileconversion.s_ms2deisotope.b_ms2deisotope = true;
			fileconversion.s_ms2denoise.b_ms2denoise = true;
			fileconversion.s_ms2denoise.st_ms2denoise_windowwidth = 30;
			fileconversion.s_ms2denoise.st_ms2denoise_peaksinwindow = 6;
			fileconversion.s_chargestatepredictor.b_chargestatepredictor = true;
			fileconversion.s_chargestatepredictor.b_chargestatepredictor = true;
			fileconversion.s_chargestatepredictor.st_chargestatepredictor_mincharge = 2;
			fileconversion.s_chargestatepredictor.st_chargestatepredictor_maxcharge = 3;
			fileconversion.s_chargestatepredictor.d_chargestatepredictor_chargefraction = 0.9;
		}
		return fileconversion;
	}

	void sys_msconvert(string par_str_msconvert_command, string par_str_filesystem_directory) {
		std::cout << "\n\n";
		string string_system = "CD Z:\\Lukah_Dykes\\IgFamily\\ProteoWizard\\";
		string_system += " && ";
		string_system += par_str_msconvert_command;
		system(string_system.c_str());
	}
}

#endif