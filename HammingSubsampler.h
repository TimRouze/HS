#ifndef HAM
#define HAM



#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <filesystem>
#include "include/Eigen/Dense"
#include "utils.h"


using namespace std;
using namespace Eigen;


class Hammer {
  public:
  //CONSTANTS
  	uint64_t k, r, cpt, nb_kmer_seen, nb_hamming, nb_1_error, nb_2_error;
    //VARIABLES
    // unordered_set<string> diff_kmer_seen;
    unordered_set<uint64_t> nb_kmer_saved;
    string hammed_file;
    Eigen::MatrixXd parity_m;
    Hammer(uint64_t ir){
        r=ir;
        nb_kmer_seen = 0;
        nb_hamming = 0;
        nb_1_error = 0;
        nb_2_error = 0;
        uint64_t size = pow(2,r)-(r+1);
        k = pow(2, r)/2;
        Eigen::MatrixXd tmp(size, r+1);
        parity_m = tmp;
        hammed_file = "";
        cpt = 0;
    }

    uint64_t findInMat(Eigen::RowVectorXd& res);
    void create_matrix();
    void generateControlMatrix(Eigen::RowVectorXd vec, int i);
    void parse_fasta(const string& input_file, const string& output_prefix, const char* file_type);
    string extract_name(const string& str);
};

#endif