// * * fpf_denovo.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_DENOVO
#define	FPF_DENOVO

#include <cstdlib>	
#include <iostream>
#include <string>
#include <vector>

#include "IgFamily.h"
#include "fpf_filesystem.h"


namespace fpf_denovo {

	using std::string;
	using std::vector;

	using fpf_filesystem::filesystem;

	void perform_novor_denovo (filesystem& par_filesystem) {
		std::cout << "\n\n\n * * * calling novor.bat * * * \n\n\n";
		string system_command_novor{};
		system_command_novor += IgFamily::DEFAULT_NOVOR_DIRECTORY;
		system_command_novor += "novor.bat -p \"";
		system_command_novor += IgFamily::DEFAULT_NOVOR_DIRECTORY;
		system_command_novor += "parameters.txt\" \"";
		system_command_novor += IgFamily::DEFAULT_IGFAMILY_DIRECTORY;
		system_command_novor += par_filesystem.directory;
		system_command_novor += par_filesystem.filename;
		system_command_novor += ".mgf\" -f";
		system(system_command_novor.c_str());
		system_command_novor = "";
		system_command_novor += "CD \"";
		system_command_novor += IgFamily::DEFAULT_IGFAMILY_DIRECTORY;
		system_command_novor += par_filesystem.directory;
		system_command_novor += "\" && del denovo_peptides_NOVOR.csv && rename ";
		system_command_novor += par_filesystem.filename;
		system_command_novor += ".mgf.csv ";
		system_command_novor += "denovo_peptides_NOVOR.csv";
		system(system_command_novor.c_str());
		std::cout << "\n\n\n * * * closing novor.bat * * * ";
	}
}

#endif