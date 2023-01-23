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
#include "include/Eigen/Dense"
#include "HammingSubsampler.h"
#include "utils.h"
#include "include/zstr.hpp"

Eigen::MatrixXd create_matrix(uint64_t& k){
	Eigen::MatrixXd right_m(k, k);
	Eigen::MatrixXd left_m = Eigen::MatrixXd::Identity(k, k);
	for(int i = 0; i < k; ++i){
		Eigen::VectorXd curr_col = right_m.col(i);
		int pos(0), j(pow(2,i)-1);
		while(pos<k){
			pos = j;
			for(; j < pos+pow(2, i) and j < k; ++j){
				curr_col(j) = 1;
			}
			j += pow(2,i);
		}
		right_m.col(i) = curr_col;
	}
	cout << right_m << endl;
	cin.get();
	cout << left_m << endl;
	cin.get();
	Eigen::MatrixXd res(left_m.rows()+right_m.rows(), left_m.cols());
	res << left_m, right_m;
	cout << res << endl;
	cin.get();
	return res;
}

// MAT *m_fill( MAT *A, double x)
// /* MAT *m_random_fill( MAT *A ) */
// {
//    int i, j;
//    for ( i = 0; i < A->m; i++ ) for ( j = 0; j < A->n; j++ )
//       { A->me[i][j] = x ; }
//       /* { A->me[i][j] = m_random ; } */
//    return A;
// }

string extract_name(const string& str){
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

void parse_fasta(const string& input_file, const string& output_prefix, uint64_t& k) {
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
    string subsampled_file=output_prefix +clean_input_file+".gz";
	zstr::ofstream* out_file_skmer = (new zstr::ofstream(subsampled_file,21,9));
	string ref, useless;
	Eigen::RowVectorXd kmer_vect(k*2), res(k*2);
	Eigen::MatrixXd parity_m = create_matrix(k);
	while (not input_stream->eof()) {
		ref = "";
		Biogetline(input_stream,ref,'A',k);
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
				kmer_vect = str2vect(ref.substr(i,k), k);
				res = kmer_vect * parity_m;
				//res = (parity_m.colwise() * kmer_vect).colwise().redux(logical_xor());
				cout << res << endl;
				cin.get();
			}
		}
	}
}

int main(int argc, char** argv) {
	char ch;
	string input, inputfof, query, output("hammed_");
	uint64_t k(31);
	uint c(8);
    bool verbose=true;

	while ((ch = getopt(argc, argv, "hdag:q:k:m:n:s:t:b:e:f:i:o:v:")) != -1) {
		switch (ch) {
			case 'i': input = optarg; break;
			case 'f': inputfof = optarg; break;
			case 'k': k = stoi(optarg); break;
			case 't': c = stoi(optarg); break;
			case 'o': output = optarg; break;
            case 'v': verbose = stoi(optarg); break;
		}
	}
	if ((input == "" && inputfof == "")) {
		cout << "Core arguments:" << endl
		     << "	-i Input file" << endl
			 << "	-f Input file of file" << endl
			 << "	-o Output prefix (hammed)" << endl
		     << "	-k Kmer size used  (31) " << endl
             << "	-t Threads used  (8) " << endl
             << "	-v Verbose level (1) " << endl
             ;
		return 0;
	}else{
        cout<<" I use k="<<k<<endl;	
		if(input != ""){
			parse_fasta(input, output, k);
			if(verbose){
                cout << "VERBOSE" << endl;
            }
		}
        
    }
	    
}
