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
  	uint64_t k, r, cpt, nb_kmer_seen, nb_kmer_saved, nb_hamming;
    //VARIABLES
    Eigen::MatrixXd parity_m;
    Hammer(uint64_t ir){
        r=ir;
        nb_kmer_seen = 0;
        nb_kmer_saved = 0;
        nb_hamming = 0;
        uint64_t size = pow(2,r)-(r+1);
        k = pow(2, r)/2;
        Eigen::MatrixXd tmp(size, r+1);
        parity_m = tmp;
        cpt = 0;
    }

    uint64_t findInMat(Eigen::RowVectorXd& res);
    void create_matrix();
    void generateControlMatrix(Eigen::RowVectorXd vec, int i);
    void parse_fasta(const string& input_file, const string& output_prefix);
    string extract_name(const string& str);
};

#endif