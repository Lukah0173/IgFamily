// * * fpf_dirichlet_mixture_model.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_DIRICHLET_MIXTURE_MODEL
#define	FPF_DIRICHLET_MIXTURE_MODEL
#include <cstdlib>						// provides - size_t
#include <string>						// provides - std::string
#include <vector>						// provides - std::vector
#include "fpf_data.h"
#include "fpf_filesystem.h"
#include "fpf_filesystem_analysis.h"
#include "fpf_blastp_analysis.h"

/*Eigen*/


namespace fpf_dirichlet_mixture_model {

	typedef size_t size_type;
	typedef std::string string_type;
	typedef fpf_data::s_multinomial_element_data s_multinomial_element_data;

	struct def_s_model_data;
	struct def_s_model_parameters;

	struct def_s_model_data {
		def_s_model_data() {
			st_N = size_type();
			st_S = size_type();
			v2_d_mnom_val = std::vector<std::vector<double>>();
			v_str_mnom_xlabel = std::vector<string_type>();
			v_str_mnom_ylabel = std::vector<string_type>();
		};

	public:
		size_type st_N;									// number of samples
		size_type st_S;									// number of multinomial dimensions
		std::vector<string_type> v_str_mnom_xlabel;		// sample filename
		std::vector<string_type> v_str_mnom_ylabel;		// gene family
		std::vector<std::vector<double>> v2_d_mnom_val; // multinomial values
	};

	struct def_s_model_parameters {
		def_s_model_parameters() {
			st_mix_comp = size_type();
			str_inputfile = string_type();
		};

	public:
		size_type st_mix_comp;		// mixture components
		string_type str_inputfile;	// input file
	};

	void initialise() {
		std::cout << "\n\n\n\n" << " ...initialising Dirichlet mixture model";
		//std::cout << "\n\n" << " Dirichlet mixture components? (Recommended: 5)";
		std::cout << "\n\n";
		//std::cin >> s_model_parameters.st_mix_comp;
	}

	//def_s_model_data create_s_model_data(std::vector<s_filesystem_blastp> par_v_blastp_filesystem_data, 
	//									 std::vector<s_multinomial_element_data> par_v_s_multinomial_element_data) {
	//	def_s_model_data con_s_model_data;
	//	con_s_model_data.st_N = par_v_blastp_filesystem_data.size();
	//	con_s_model_data.st_S = par_v_s_multinomial_element_data.size();
	//	std::vector<std::vector<double>> con_v2_d_mnom_val(con_s_model_data.st_S, std::vector<double>(con_s_model_data.st_N));
	//	for (auto itr_v_s_multinomial_element_data : par_v_s_multinomial_element_data) {
	//		con_s_model_data.v_str_mnom_ylabel.push_back(itr_v_s_multinomial_element_data.return_str_genefamily());
	//	}
	//	for (auto itr_v_blastp_filesystem_data : par_v_blastp_filesystem_data) {
	//		con_s_model_data.v_str_mnom_xlabel.push_back(itr_v_blastp_filesystem_data.str_blastp_file);
	//		for (auto itr_v_s_filesystem_mnom : itr_v_blastp_filesystem_data.v_s_filesystem_mnom) {
	//			auto find_v_str_mnom_xlabel = std::find_if(con_s_model_data.v_str_mnom_xlabel.begin(),
	//													   con_s_model_data.v_str_mnom_xlabel.end(),
	//													   [itr_v_s_filesystem_mnom](string_type par_str_mnom_xlabel) {
	//				return par_str_mnom_xlabel == itr_v_s_filesystem_mnom.str_mnom_comp;
	//			});
	//			auto index_mnom_xlabel = std::distance(con_s_model_data.v_str_mnom_xlabel.begin(), find_v_str_mnom_xlabel);
	//			con_v2_d_mnom_val[index_mnom_xlabel].push_back(itr_v_s_filesystem_mnom.d_mnom_value);
	//		}
	//	}
	//	return con_s_model_data;
	//}
}

#endif