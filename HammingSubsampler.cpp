#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <mutex>
#include <omp.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/stat.h>
#include <tmmintrin.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include <filesystem>
#include <cmath>
#include <bits/stdc++.h>
#include "include/Eigen/Dense"
#include "HammingSubsampler.h"
#include "utils.h"
#include "include/zstr.hpp"
 
void Hammer::generateControlMatrix(Eigen::RowVectorXd vec, int i){
    if (i == r) {
		Eigen::RowVectorXd curr_row = parity_m.row(cpt);
		if(vec.sum() != 0 and vec.sum() != 1){
			for(int j = 0; j < vec.size(); ++j){
				curr_row(j) = vec(j);
			}
			parity_m.row(cpt) = curr_row;
			++cpt;
		}
        return;
    }
    vec(i) = 0;
    generateControlMatrix(vec, i + 1);
    vec(i) = 1;
    generateControlMatrix(vec, i + 1);
}

void Hammer::create_matrix(){
    Eigen::RowVectorXd vec(r);
	uint64_t size = pow(2, r);
	uint64_t size_tmp = pow(2, r) - (r+1);
	uint64_t sum(0);
    generateControlMatrix(vec, 0);
	for(auto row: parity_m.rowwise()){
		sum = row.sum();
		if(sum%2 == 0){
			row(r) = 1;
		}
	}
    Eigen::MatrixXd id_m = Eigen::MatrixXd::Identity(r+1, r+1);
	Eigen::MatrixXd tmp(size_tmp, r+1);
	tmp = parity_m;
	parity_m.resize(size, r+1);
	parity_m << tmp, id_m;

}

uint64_t Hammer::findInMat(Eigen::RowVectorXd& res){
	for(int i = 0; i < parity_m.rows(); ++i){
  		if(parity_m.row(i) == res){
			return i;
		}
	}
	return -1;
}

string Hammer::extract_name(const string& str){
    string result;
    uint64_t begin(0);
    for(uint i(0);i<str.size();++i){
        if(str[i]=='/'){
            begin=i+1;
        }
    }
    for(uint i(begin);i<str.size();++i){
        if(str[i]=='.'){
            return result;
        }else{
            result.push_back(str[i]);
        }
    }
    return result;
}

void Hammer::parse_fasta(const string& input_file, const string& output_prefix, const char* file_type) {
	cout << file_type << endl;
    uint64_t total_kmer_number(0), selected_kmer_number(0);
    uint64_t read_kmer(0);
	string tmp;
	zstr::ifstream* input_stream = openFile(input_file);
    if(input_stream==NULL){
        cout<<"Can't open file: "<<input_file<<endl;
        return;
    }
    if(not input_stream->good()){
        cout<<"Can't open file: "<<input_file<<endl;
        return;
    }
    string clean_input_file=extract_name(input_file);
    hammed_file=output_prefix +clean_input_file+".gz";
	zstr::ofstream* out_file = (new zstr::ofstream(hammed_file,21,9));
	string ref, useless, curr_kmer;
	unordered_map<string,uint64_t> sketch;
	Eigen::RowVectorXd kmer_vect(k*2), res(r+1), hamming(r+1);
	create_matrix();
	// cout << parity_m << endl;
	map<vector<int>, pair<uint64_t, uint64_t>> combinations;
	// map<string, set<pair<char, char>>> unique_kmer_per_HW;
	// int po = 0;
	for (uint64_t i = 0; i < cpt+r ; ++i){
		for (uint64_t j = i+1; j < cpt+r+1; ++j){
				std::vector<int> res;
			for (uint64_t l = 0; l < parity_m.row(i).size(); ++l){
				res.push_back((int)(parity_m.row(i)(l)+ parity_m.row(j)(l))%2);
			}
			// po++;
			if (combinations.find(res) == combinations.end()){
				pair<uint64_t, uint64_t> pos;
				pos.first = i;
				pos.second = j;
				combinations[res] = pos;
			}
		}
	}
	cout << "Matrix creation done. Starting computations..." << endl;
//cout << parity_m << endl;
// cout << po << endl;
// cout << combinations.size() << endl;
	while (not input_stream->eof()) {
		ref = "";
		Biogetline(input_stream,ref,*file_type,k);
		if (ref.size() < k) {
			ref = "";
		} else {
			read_kmer += ref.size() - k + 1;
		}

		// FOREACH sequence
		if (not ref.empty()) {
			bool is_rev, old_rev;
			uint64_t last_position(0);
			uint64_t i(0);
			// FOREACH KMER
			for (; i + k < ref.size(); ++i) {
				curr_kmer = canonize(ref.substr(i, k));
				kmer_vect = str2vect(curr_kmer, k);
				nb_kmer_seen++;
				res = (kmer_vect * parity_m).unaryExpr([](int x){return (double)(x%2);});
				if(res == hamming){
					//STORE K-MER & INCREMENTE COMPTEUR
					// cout << "mot de hamming" << endl;
					if (sketch.find(curr_kmer) == sketch.end()) {
						sketch[curr_kmer] = 1;
						nb_hamming++;
					}else{	
						sketch[curr_kmer]++;
					}
					// pair<char, char> positions;
					// positions.first = -1;
					// positions.second = -1;
					// unique_kmer_per_HW[curr_kmer].insert(positions);
				}else{
					// FINDINMAT RETOURNE LA POS DANS LA MATRICE OU RENVOIE -1 SINON
					uint64_t pos = findInMat(res);
					if(pos != -1){
						// SWITCH LE BIT A LA POS DE L'ERREUR
						// STORE K-MER & INCREMENTE COMPTEUR
						kmer_vect(pos) = (double)(((int)kmer_vect(pos)+1)%2);
						nb_1_error++;
						string k_mer = vect2strv(kmer_vect);
						// pair<char, char> positions;
						// positions.first = pos;
						// positions.second = -1;
						// unique_kmer_per_HW[curr_kmer].insert(positions);
						if (sketch.find(k_mer) == sketch.end()) {
							sketch[k_mer] = 1;
							nb_hamming++;
						}else{
							sketch[k_mer]++;
						}			
					}
					else{
						// CHECK SI 2 ERREURS
						std::vector<int> result;
                        for (uint64_t l = 0; l < res.size(); ++l){
                            int vect = res(l);
                            result.push_back(vect);
                        }
                        if (combinations.find(result) != combinations.end()){
							nb_2_error++;
                            pair<uint64_t, uint64_t> pos = combinations[result];
                            kmer_vect(pos.first) = (double)(((int)kmer_vect(pos.first)+1)%2);
                            kmer_vect(pos.second) = (double)(((int)kmer_vect(pos.second)+1)%2);
							string k_mer = vect2strv(kmer_vect);
							// unique_kmer_per_HW[curr_kmer].insert(pos);
                            // if (nb_kmer_saved.find(unrevhash(str2num(curr_kmer))) == nb_kmer_saved.end()) {
  							// 	nb_kmer_saved.insert(unrevhash(str2num(curr_kmer)));
							// }
                            if (sketch.find(k_mer) == sketch.end()) {
								sketch[k_mer] = 1;
								nb_hamming++;
							}else{
								sketch[k_mer]++;
							}
                        }
					}
				}
			}
		}
	}
	cout << "Writing output..." << endl;
	cout << "I have seen " << intToString(read_kmer) << endl;
	cout << "nb hamming word??? " << intToString(sketch.size()) << endl;
	string line_1 = "";
	// map<string, pair<uint64_t, uint64_t>>tmp_res;
	// pair<uint64_t, uint64_t> tmp_pair;
	// for(auto const& [h_word, pos] : unique_kmer_per_HW){
	// 	tmp_pair.first = pos.size();
	// 	tmp_res[h_word] = tmp_pair;
	// }
	for(auto const& [h_word, nb]: sketch){
		line_1 = ">" + intToString(nb) + "\n";
		out_file->write(line_1.c_str(), line_1.size());
		out_file->write(h_word.c_str(), h_word.size());
		out_file->write("\n", 1);
	}
	// for(auto const& [h_word, paire]: tmp_res){
	// 	line_1 = ">" + intToString(paire.second) + " | " = intToString(paire.first) + "\n";
	// 	out_file->write(line_1.c_str(), line_1.size());
	// 	out_file->write(h_word.c_str(), h_word.size());
	// 	out_file->write("\n", 1);
	// }
	delete input_stream;
	delete out_file;
}

int main(int argc, char** argv) {
char ch, *file_type = new char('A');
	string input, inputfof, query, output("hammed_");
	uint64_t k(31), r(5);
	uint c(8);
    bool verbose=true;

	while ((ch = getopt(argc, argv, "hdag:q:k:r:m:n:s:t:b:e:f:i:o:v:t:")) != -1) {
		switch (ch) {
			case 'i': input = optarg; break;
			case 'f': inputfof = optarg; break;
			case 'r': r = stoi(optarg); break;
			case 'k': k = stoi(optarg); break;
			case 'c': c = stoi(optarg); break;
			case 'o': output = optarg; break;
            case 'v': verbose = stoi(optarg); break;
			case 't': file_type = optarg; break;
		}
	}
	if ((input == "" && inputfof == "")) {
		cout << "Core arguments:" << endl
		     << "	-i Input file" << endl
			 << "	-f Input file of file" << endl
			 << "	-o Output prefix (hammed)" << endl
			 << "	-t File type (A for fasta, Q for fastq)" << endl
		     << "	-k Kmer size used  (31) " << endl
             << "	-t Threads used  (8) " << endl
             << "	-v Verbose level (1) " << endl
             ;
		return 0;
	}else{
        cout<<" I use r="<<r<<" which means k=" << (pow(2,r))/2<<endl;	
		if(input != ""){
			Hammer h = Hammer(r);
			h.parse_fasta(input, output, file_type);
			if(verbose){
				cout << "I have seen " << intToString(h.nb_kmer_seen) << " k-mers"/*, among which " << intToString(h.diff_kmer_seen.size()) << " unique k-mers and I saved " << intToString(h.nb_kmer_saved.size()) <<*/ " under " << intToString(h.nb_hamming) << " Hamming words." << endl;
                // cout << "This means " << (double)h.nb_kmer_saved.size()/h.nb_hamming << " k-mer per Hamming words in average" << endl;
                // cout << "This means a subsampling rate of " /*<< (double)h.diff_kmer_seen.size()/h.nb_kmer_saved.size() << " Or "*/ << (double)h.nb_kmer_saved.size()/h.nb_hamming << " k-mers per hamming word." << endl;
				cout << "Output file is " << intToString(std::filesystem::file_size(h.hammed_file)/1000) << "KB" << endl;
				cout << "Input file is " << intToString(std::filesystem::file_size(input)/1000) << "KB" << endl;
				cout << "There are " << intToString(h.nb_1_error) << " k-mers at 1 error distance from their Hamming word and " << intToString(h.nb_2_error) << " k-mers at 2 error distance." << endl;
            }
		}
        
    }
	    
}
