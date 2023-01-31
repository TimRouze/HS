#ifndef UTIL
#define UTIL
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <chrono>
#include <omp.h>
#include <tmmintrin.h>
#include <math.h>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include "include/zstr.hpp"
#include "include/Eigen/Dense"

#define kmer uint64_t


using namespace std;
using Eigen::RowVectorXd;

uint64_t canonize(uint64_t x, uint64_t n);
uint64_t rcbc(uint64_t in, uint64_t n);
void updateK(kmer & min, char nuc, uint64_t& k);
kmer nuc2int(char c);
kmer nuc2intrc(char c);
string intToString(uint64_t n);
kmer min_k(const kmer& k1, const kmer& k2);
kmer str2num(const string& str);
string num2str(uint64_t num,uint k);
uint64_t revhash(uint64_t x);
uint64_t unrevhash(uint64_t x);
char revCompChar(char c);
string revComp(const string& s);
vector<bool> str2boolv(const string& str);
RowVectorXd str2vect(const string& str, uint64_t& k);
string vect2strv(const RowVectorXd& v);
string bool2strv(const vector<bool>& v);
zstr::ifstream* openFile(const string& input_file);
void Biogetline(zstr::ifstream* in,string& result,char type,uint K);

#endif