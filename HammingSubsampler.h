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


Eigen::MatrixXd create_matrix(uint64_t& k);
void parse_fasta(const string& input_file, const string& output_prefix, uint64_t& k);
string extract_name(const string& str);

#endif