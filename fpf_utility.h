
// * * fpf_utility.h * * 
// 
// Lukah Dykes - Flinders Proteomics Facility - 2016
// 
// * * * * *



#ifndef FPF_UTILITY
#define	FPF_UTILITY

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


namespace fpf_utility {
	
	using std::string;
	using std::vector;

	struct transcript {
	public: 
		size_t transcript_key;
		string transcript_sequence;
		string translation_sequence_frame_1;
		string translation_sequence_frame_2;
		string translation_sequence_frame_3;
	};

	string parse_transcript_data() {
		string input_transcript_data = "transcript.csv";
		std::ifstream fin_transcript_data(input_transcript_data);
		char get_transcript_data{};
		string transcript_sequence{};
		while (fin_transcript_data.get(get_transcript_data)) {
			if (get_transcript_data != '\n') {
				transcript_sequence += get_transcript_data;
			}
		}
		return transcript_sequence;
	}

	void translate_read_codon(const string& par_codon, string& par_translation_sequence_frame) {
		if (par_codon == "TTT") {
			par_translation_sequence_frame += "F";
		}
		if (par_codon == "TTC") {
			par_translation_sequence_frame += "F";
		}
		if (par_codon == "TTA") {
			par_translation_sequence_frame += "L";
		}
		if (par_codon == "TTG") {
			par_translation_sequence_frame += "L";
		}
		if (par_codon == "TCT") {
			par_translation_sequence_frame += "S";
		}
		if (par_codon == "TCC") {
			par_translation_sequence_frame += "S";
		}
		if (par_codon == "TCA") {
			par_translation_sequence_frame += "S";
		}
		if (par_codon == "TCG") {
			par_translation_sequence_frame += "S";
		}
		if (par_codon == "TAT") {
			par_translation_sequence_frame += "Y";
		}
		if (par_codon == "TAC") {
			par_translation_sequence_frame += "Y";
		}
		if (par_codon == "TAA") {
			par_translation_sequence_frame += "*";
		}
		if (par_codon == "TAG") {
			par_translation_sequence_frame += "*";
		}
		if (par_codon == "TGT") {
			par_translation_sequence_frame += "C";
		}
		if (par_codon == "TGC") {
			par_translation_sequence_frame += "C";
		}
		if (par_codon == "TGA") {
			par_translation_sequence_frame += "*";
		}
		if (par_codon == "TGG") {
			par_translation_sequence_frame += "W";
		}
		if (par_codon == "CTT") {
			par_translation_sequence_frame += "L";
		}
		if (par_codon == "CTC") {
			par_translation_sequence_frame += "L";
		}
		if (par_codon == "CTA") {
			par_translation_sequence_frame += "L";
		}
		if (par_codon == "CTG") {
			par_translation_sequence_frame += "L";
		}
		if (par_codon == "CCT") {
			par_translation_sequence_frame += "P";
		}
		if (par_codon == "CCC") {
			par_translation_sequence_frame += "P";
		}
		if (par_codon == "CCA") {
			par_translation_sequence_frame += "P";
		}
		if (par_codon == "CCG") {
			par_translation_sequence_frame += "P";
		}
		if (par_codon == "CAT") {
			par_translation_sequence_frame += "H";
		}
		if (par_codon == "CAC") {
			par_translation_sequence_frame += "H";
		}
		if (par_codon == "CAA") {
			par_translation_sequence_frame += "Q";
		}
		if (par_codon == "CAG") {
			par_translation_sequence_frame += "Q";
		}
		if (par_codon == "CGT") {
			par_translation_sequence_frame += "R";
		}
		if (par_codon == "CGC") {
			par_translation_sequence_frame += "R";
		}
		if (par_codon == "CGA") {
			par_translation_sequence_frame += "R";
		}
		if (par_codon == "CGG") {
			par_translation_sequence_frame += "R";
		}
		if (par_codon == "ATT") {
			par_translation_sequence_frame += "I";
		}
		if (par_codon == "ATC") {
			par_translation_sequence_frame += "I";
		}
		if (par_codon == "ATA") {
			par_translation_sequence_frame += "I";
		}
		if (par_codon == "ATG") {
			par_translation_sequence_frame += "M";
		}
		if (par_codon == "ACT") {
			par_translation_sequence_frame += "T";
		}
		if (par_codon == "ACC") {
			par_translation_sequence_frame += "T";
		}
		if (par_codon == "ACA") {
			par_translation_sequence_frame += "T";
		}
		if (par_codon == "ACG") {
			par_translation_sequence_frame += "T";
		}
		if (par_codon == "AAT") {
			par_translation_sequence_frame += "N";
		}
		if (par_codon == "AAC") {
			par_translation_sequence_frame += "N";
		}
		if (par_codon == "AAA") {
			par_translation_sequence_frame += "K";
		}
		if (par_codon == "AAG") {
			par_translation_sequence_frame += "K";
		}
		if (par_codon == "AGT") {
			par_translation_sequence_frame += "S";
		}
		if (par_codon == "AGC") {
			par_translation_sequence_frame += "S";
		}
		if (par_codon == "AGA") {
			par_translation_sequence_frame += "R";
		}
		if (par_codon == "AGG") {
			par_translation_sequence_frame += "R";
		}
		if (par_codon == "GTT") {
			par_translation_sequence_frame += "V";
		}
		if (par_codon == "GTC") {
			par_translation_sequence_frame += "V";
		}
		if (par_codon == "GTA") {
			par_translation_sequence_frame += "V";
		}
		if (par_codon == "GTG") {
			par_translation_sequence_frame += "V";
		}
		if (par_codon == "GCT") {
			par_translation_sequence_frame += "A";
		}
		if (par_codon == "GCC") {
			par_translation_sequence_frame += "A";
		}
		if (par_codon == "GCA") {
			par_translation_sequence_frame += "A";
		}
		if (par_codon == "GCG") {
			par_translation_sequence_frame += "A";
		}
		if (par_codon == "GAT") {
			par_translation_sequence_frame += "D";
		}
		if (par_codon == "GAC") {
			par_translation_sequence_frame += "D";
		}
		if (par_codon == "GAA") {
			par_translation_sequence_frame += "E";
		}
		if (par_codon == "GAG") {
			par_translation_sequence_frame += "E";
		}
		if (par_codon == "GGT") {
			par_translation_sequence_frame += "G";
		}
		if (par_codon == "GGC") {
			par_translation_sequence_frame += "G";
		}
		if (par_codon == "GGA") {
			par_translation_sequence_frame += "G";
		}
		if (par_codon == "GGG") {
			par_translation_sequence_frame += "G";
		}
	}

	transcript translate_transcript() {
		transcript temp_transcript{};
		temp_transcript.transcript_sequence = parse_transcript_data();;
		size_t count_transcript_sequence{};
		string transcript_codon_frame_1{};
		string transcript_codon_frame_2{};
		string transcript_codon_frame_3{};
		for(const auto& itr_transcript_sequence : temp_transcript.transcript_sequence) {
			++count_transcript_sequence;
			transcript_codon_frame_1 += itr_transcript_sequence;
			transcript_codon_frame_2 += itr_transcript_sequence;
			transcript_codon_frame_3 += itr_transcript_sequence;
			if(count_transcript_sequence % 3 == 0) { 
				translate_read_codon(transcript_codon_frame_1, temp_transcript.translation_sequence_frame_1);
				transcript_codon_frame_1.clear();
			}
			if (count_transcript_sequence % 3 == 1) {
				translate_read_codon(transcript_codon_frame_2, temp_transcript.translation_sequence_frame_2);
				transcript_codon_frame_2.clear();
			}
			if (count_transcript_sequence % 3 == 2) {
				translate_read_codon(transcript_codon_frame_3, temp_transcript.translation_sequence_frame_3);
				transcript_codon_frame_3.clear();
			}
		}
		return temp_transcript;
	}

	void fout_transcript_and_translation (const transcript& par_transcript) {
		std::string output_transcript = "transcript_and_translation.csv";
		std::ofstream fout_transcript;
		fout_transcript.open(output_transcript);
		//fout_transcript << par_transcript.transcript_key;
		fout_transcript << par_transcript.transcript_sequence << ",";
		fout_transcript << par_transcript.translation_sequence_frame_1 << ",";
		fout_transcript << par_transcript.translation_sequence_frame_2 << ",";
		fout_transcript << par_transcript.translation_sequence_frame_3 << ",";
	}
}

#endif

