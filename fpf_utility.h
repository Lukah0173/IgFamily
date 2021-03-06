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
		size_t key_transcript;
		string transcript_sequence;
		string translation_sequence_5to3_frame_1;
		string translation_sequence_5to3_frame_2;
		string translation_sequence_5to3_frame_3;
		string translation_sequence_3to5_frame_1;
		string translation_sequence_3to5_frame_2;
		string translation_sequence_3to5_frame_3;
	};

	struct sample_transcript_and_translation {
	public:
		size_t key_sample;
		vector<transcript> sample_v_transcript;
	};

	inline double log_base(double d, double base) {
		return (log(d) / log(base));
	}

	vector<string> parse_transcript_directory() {
		vector<string> temp_v_transcript_directory{};
		string temp_transcript_directory{};
		string input_transcript_directory = IgFamily::DEFAULT_TRANSCRIPT_DIRECTORY + "transcript_directory.txt";
		std::ifstream fin_transcript_directory(input_transcript_directory);
		char get_transcript_directory{};
		while (fin_transcript_directory.get(get_transcript_directory)) {
			if ((get_transcript_directory != ',') && (get_transcript_directory != '\n')) {
				temp_transcript_directory += get_transcript_directory;
			}
			if (get_transcript_directory == ',') {
				temp_transcript_directory = IgFamily::DEFAULT_TRANSCRIPT_DIRECTORY + temp_transcript_directory;
				temp_v_transcript_directory.push_back(temp_transcript_directory);
				temp_transcript_directory.clear();
			}
		}
		return temp_v_transcript_directory;
	}

	vector<sample_transcript_and_translation> parse_transcript_data() {
		vector<sample_transcript_and_translation> temp_v_transcript_and_translation{};
		sample_transcript_and_translation temp_transcript_and_translation{};
		vector<string> v_transcript_directory = parse_transcript_directory();
		size_t count_transcript{};
		for (const auto& itr_v_transcript_directory : v_transcript_directory) {
			vector<transcript> temp_v_transcript{};
			transcript temp_transcript{};
			std::ifstream fin_transcript_data(itr_v_transcript_directory);
			char get_transcript_data{};
			string transcript_sequence{};
			bool ignore_header{};
			bool read_transcript{};
			while (fin_transcript_data.get(get_transcript_data)) {
				if (get_transcript_data == '\n') {
					if (!ignore_header) {
						ignore_header = true;
					}
					read_transcript = false;
				}
				if (ignore_header && !read_transcript) {
					if (get_transcript_data == ' ') {
						temp_transcript.transcript_sequence = transcript_sequence;
						temp_v_transcript.push_back(temp_transcript);
						read_transcript = true;
						transcript_sequence.clear();
						++count_transcript;
						if (count_transcript % 100 == 0) {
							std::cout << "\n transcripts parsed: " << count_transcript;
						}
					}
					if (get_transcript_data != '\n') {
						transcript_sequence += get_transcript_data;
					}
				}
			}
			temp_transcript_and_translation.sample_v_transcript = temp_v_transcript;
			temp_v_transcript_and_translation.push_back(temp_transcript_and_translation);
		}
		return temp_v_transcript_and_translation;
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
			par_translation_sequence_frame += "-";
		}
		if (par_codon == "TAG") {
			par_translation_sequence_frame += "-";
		}
		if (par_codon == "TGT") {
			par_translation_sequence_frame += "C";
		}
		if (par_codon == "TGC") {
			par_translation_sequence_frame += "C";
		}
		if (par_codon == "TGA") {
			par_translation_sequence_frame += "-";
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

	void translate_v_transcript(vector<sample_transcript_and_translation>& par_v_sample_transcript_and_translation) {
		size_t count_translation{};
		for (auto& itr_v_sample_transcript_and_translation : par_v_sample_transcript_and_translation) {
			for (auto& itr_v_transcript : itr_v_sample_transcript_and_translation.sample_v_transcript) {
				size_t count_transcript_sequence{};
				string transcript_codon_frame_1{};
				string transcript_codon_frame_2{};
				string transcript_codon_frame_3{};
				string transcript_codon_reverseframe_1{};
				string transcript_codon_reverseframe_2{};
				string transcript_codon_reverseframe_3{};
				for (const auto& itr_transcript_sequence : itr_v_transcript.transcript_sequence) {
					++count_transcript_sequence;
					transcript_codon_frame_1 += itr_transcript_sequence;
					transcript_codon_frame_2 += itr_transcript_sequence;
					transcript_codon_frame_3 += itr_transcript_sequence;
					if (count_transcript_sequence % 3 == 0) {
						translate_read_codon(transcript_codon_frame_1, itr_v_transcript.translation_sequence_5to3_frame_1);
						transcript_codon_frame_1.clear();
					}
					if (count_transcript_sequence % 3 == 1) {
						translate_read_codon(transcript_codon_frame_2, itr_v_transcript.translation_sequence_5to3_frame_2);
						transcript_codon_frame_2.clear();
					}
					if (count_transcript_sequence % 3 == 2) {
						translate_read_codon(transcript_codon_frame_3, itr_v_transcript.translation_sequence_5to3_frame_3);
						transcript_codon_frame_3.clear();
					}
				}
				for (string::reverse_iterator& itr_transcript_sequence = itr_v_transcript.transcript_sequence.rbegin(); itr_transcript_sequence != itr_v_transcript.transcript_sequence.rend(); ++itr_transcript_sequence) {
					++count_transcript_sequence;
					if (*itr_transcript_sequence == 'A') {
						transcript_codon_reverseframe_1 += 'T';
						transcript_codon_reverseframe_2 += 'T';
						transcript_codon_reverseframe_3 += 'T';
					}
					if (*itr_transcript_sequence == 'T') {
						transcript_codon_reverseframe_1 += 'A';
						transcript_codon_reverseframe_2 += 'A';
						transcript_codon_reverseframe_3 += 'A';
					}
					if (*itr_transcript_sequence == 'G') {
						transcript_codon_reverseframe_1 += 'C';
						transcript_codon_reverseframe_2 += 'C';
						transcript_codon_reverseframe_3 += 'C';
					}
					if (*itr_transcript_sequence == 'C') {
						transcript_codon_reverseframe_1 += 'G';
						transcript_codon_reverseframe_2 += 'G';
						transcript_codon_reverseframe_3 += 'G';
					}
					if ((count_transcript_sequence + itr_v_transcript.transcript_sequence.size() + 2) % 3 == 0) {
						translate_read_codon(transcript_codon_reverseframe_1, itr_v_transcript.translation_sequence_3to5_frame_1);
						transcript_codon_reverseframe_1.clear();
					}
					if ((count_transcript_sequence + itr_v_transcript.transcript_sequence.size() + 2) % 3 == 1) {
						translate_read_codon(transcript_codon_reverseframe_2, itr_v_transcript.translation_sequence_3to5_frame_2);
						transcript_codon_reverseframe_2.clear();
					}
					if ((count_transcript_sequence + itr_v_transcript.transcript_sequence.size() + 2) % 3 == 2) {
						translate_read_codon(transcript_codon_reverseframe_3, itr_v_transcript.translation_sequence_3to5_frame_3);
						transcript_codon_reverseframe_3.clear();
					}
				}
				++count_translation;
				if (count_translation % 100 == 0) {
					std::cout << "\n transcripts translated: " << count_translation;
				}
			}
		}
	}

	void fout_transcript_and_translation(const vector<sample_transcript_and_translation>& par_v_sample_transcript_and_translation) {
		for (auto& itr_v_sample_transcript_and_translation : par_v_sample_transcript_and_translation) {
			std::string output_transcript = "transcript_and_translation.txt";
			std::ofstream fout_transcript;
			fout_transcript.open(output_transcript);
			for (const auto& itr_v_transcript : itr_v_sample_transcript_and_translation.sample_v_transcript) {
				//fout_transcript << par_transcript.transcript_key;
				fout_transcript << itr_v_transcript.transcript_sequence << "\n\n";
				fout_transcript << itr_v_transcript.translation_sequence_5to3_frame_1 << "\n";
				fout_transcript << itr_v_transcript.translation_sequence_5to3_frame_2 << "\n";
				fout_transcript << itr_v_transcript.translation_sequence_5to3_frame_3 << "\n";
				fout_transcript << itr_v_transcript.translation_sequence_3to5_frame_1 << "\n";
				fout_transcript << itr_v_transcript.translation_sequence_3to5_frame_2 << "\n";
				fout_transcript << itr_v_transcript.translation_sequence_3to5_frame_3 << "\n";
				fout_transcript << "\n\n";
			}
		}
	}
}

#endif

