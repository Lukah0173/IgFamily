// * * fpf_dirichlet_mixture_model.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_DIRICHLET_MIXTURE_MODEL
#define	FPF_DIRICHLET_MIXTURE_MODEL
#include <cstdlib>						// provides - size_t
#include <string>						// provides - string
#include <vector>						// provides - vector
#include "fpf_data.h"
#include "fpf_filesystem.h"
#include "fpf_filesystem_analysis.h"
#include "fpf_blastp_analysis.h"

/*Eigen*/



namespace fpf_dirichlet_mixture_model {

	using std::string;
	using std::vector;

	typedef fpf_data::FASTA_category FASTA_category;

	struct def_s_model_data;
	struct def_s_model_parameters;

	struct def_s_model_data {
		def_s_model_data() {
			st_N = size_t();
			st_S = size_t();
			v2_d_mnom_val = vector<vector<double>>();
			v_str_mnom_xlabel = vector<string>();
			v_str_mnom_ylabel = vector<string>();
		};

	public:
		size_t st_N;									// number of samples
		size_t st_S;									// number of multinomial dimensions
		vector<string> v_str_mnom_xlabel;		// sample filename
		vector<string> v_str_mnom_ylabel;		// gene family
		vector<vector<double>> v2_d_mnom_val; // multinomial values
	};

	struct def_s_model_parameters {
		def_s_model_parameters() {
			st_mix_comp = size_t();
			str_inputfile = string();
		};

	public:
		size_t st_mix_comp;		// mixture components
		string str_inputfile;	// input file
	};

	void initialise() {
		std::cout << "\n\n\n\n" << " ...initialising Dirichlet mixture model";
		//std::cout << "\n\n" << " Dirichlet mixture components? (Recommended: 5)";
		std::cout << "\n\n";
		//std::cin >> s_model_parameters.st_mix_comp;
	}

	//def_s_model_data create_s_model_data(vector<s_filesystem_blastp> par_v_blastp_filesystem_data, 
	//									 vector<FASTA_category> par_v_s_multinomial_element_data) {
	//	def_s_model_data con_s_model_data;
	//	con_s_model_data.st_N = par_v_blastp_filesystem_data.size();
	//	con_s_model_data.st_S = par_v_s_multinomial_element_data.size();
	//	vector<vector<double>> con_v2_d_mnom_val(con_s_model_data.st_S, vector<double>(con_s_model_data.st_N));
	//	for (auto itr_v_s_multinomial_element_data : par_v_s_multinomial_element_data) {
	//		con_s_model_data.v_str_mnom_ylabel.push_back(itr_v_s_multinomial_element_data.return_str_genefamily());
	//	}
	//	for (auto itr_v_blastp_filesystem_data : par_v_blastp_filesystem_data) {
	//		con_s_model_data.v_str_mnom_xlabel.push_back(itr_v_blastp_filesystem_data.str_blastp_file);
	//		for (auto itr_v_filesystem_mnom : itr_v_blastp_filesystem_data.v_s_filesystem_mnom) {
	//			auto find_v_str_mnom_xlabel = std::find_if(con_s_model_data.v_str_mnom_xlabel.begin(),
	//													   con_s_model_data.v_str_mnom_xlabel.end(),
	//													   [itr_v_filesystem_mnom](string par_str_mnom_xlabel) {
	//				return par_str_mnom_xlabel == itr_v_filesystem_mnom.str_mnom_comp;
	//			});
	//			auto index_mnom_xlabel = std::distance(con_s_model_data.v_str_mnom_xlabel.begin(), find_v_str_mnom_xlabel);
	//			con_v2_d_mnom_val[index_mnom_xlabel].push_back(itr_v_filesystem_mnom.d_mnom_value);
	//		}
	//	}
	//	return con_s_model_data;
	//}
}

#endif