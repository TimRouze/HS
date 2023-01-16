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
#include "HammingSubsampler.h"
#include "utils.h"

#include "include/zstr.hpp"


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
	{
		string ref, useless;
		while (not input_stream->eof()) {
			ref = "";
			{
                Biogetline(input_stream,ref,'A',k);
				if (ref.size() < k) {
					ref = "";
				} else {
					read_kmer += ref.size() - k + 1;
				}
			}
			// FOREACH sequence
			if (not ref.empty()) {
				bool is_rev, old_rev;
				uint64_t last_position(0);
				kmer seq(str2num(ref.substr(0, k)));
				uint64_t i(0);
				// FOREACH KMER
				for (; i + k < ref.size(); ++i) {
					updateK(seq, ref[i + k], k);
                    cout << num2str(seq, k) << endl;
                    cin.get();
				}
			}
		}
	}
}

int main(int argc, char** argv) {
	char ch;
	string input, inputfof, query, output("hammed_");
	uint64_t k(32);
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
		     << "	-k Kmer size used  (32) " << endl
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
