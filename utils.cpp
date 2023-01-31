#include "utils.h"
#include <memory>
#include <map>
#include <algorithm>
#include <string>
#include <sys/stat.h>


kmer nuc2int(char c) {
	return (c / 2) % 4;
}

kmer nuc2intrc(char c) {
	return ((c / 2) % 4) ^ 2;
}

string intToString(uint64_t n) {
	if (n < 1000) {
		return to_string(n);
	}
	string end(to_string(n % 1000));
	if (end.size() == 3) {
		return intToString(n / 1000) + "," + end;
	}
	if (end.size() == 2) {
		return intToString(n / 1000) + ",0" + end;
	}
	return intToString(n / 1000) + ",00" + end;
}

char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}

string revComp(const string& s) {
	string rc(s.size(), 0);
	for (int i((int)s.length() - 1); i >= 0; i--) {
		rc[s.size() - 1 - i] = revCompChar(s[i]);
	}
	return rc;
}

kmer str2num(const string& str) {
	kmer res(0);
	for (uint64_t i(0); i < str.size(); i++) {
		res <<= 2;
		res += (str[i] / 2) % 4;
	}
	return res;
}


string num2str(uint64_t num,uint k){
	string str;
	int nuc;
	for(uint i(0);i<k;++i){
		nuc=num%4;
		switch (nuc){
			case 0:str.push_back('A');break;
			case 1:str.push_back('C');break;
			case 2:str.push_back('T');break;
			case 3:str.push_back('G');break;
		}
		num>>=2;
	}
	reverse( str.begin(), str.end());
	return str;
}

static kmer hashtest( kmer u ){
  kmer v = u * 3935559000370003845 + 2691343689449507681;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >>  4;
  v *= 4768777513237032717;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v <<  5;
  return v;
  }


uint64_t unrevhash(uint64_t x) {
	//~ return murmur64(x);
    return hashtest(x);
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x);
	return x;
}

RowVectorXd str2vect(const string& str, uint64_t& k) {
	Eigen::RowVectorXd res(k*2);
	uint64_t j(0);
	for (uint64_t i(0); i < str.size(); ++i) {
		if (str[i] == 'G' or str[i] == 'T') {
			res(j) = 1;
		} else {
			res(j) = 0;
		}
		if (str[i] == 'C' or str[i] == 'G') {
			res(j+1) = 1;
		} else {
			res(j+1) = 0;
		}
		j += 2;
	}
	return res;
}

string vect2strv(const RowVectorXd& v) {
	string res;
	for (uint64_t i(0); i < v.size(); i += 2) {
		if (v(i)) {
			if (v(i + 1)) {
				res += 'G';
			} else {
				res += 'T';
			}
		} else {
			if (v(i + 1)) {
				res += 'C';
			} else {
				res += 'A';
			}
		}
	}
	return res;
}

vector<bool> str2boolv(const string& str) {
	vector<bool> res;
	for (uint64_t i(0); i < str.size(); ++i) {
		if (str[i] == 'G' or str[i] == 'T') {
			res.push_back(true);
		} else {
			res.push_back(false);
		}
		if (str[i] == 'C' or str[i] == 'G') {
			res.push_back(true);
		} else {
			res.push_back(false);
		}
	}
	return res;
}


string bool2strv(const vector<bool>& v) {
	string res;
	for (uint64_t i(0); i < v.size(); i += 2) {
		if (v[i]) {
			if (v[i + 1]) {
				res += 'G';
			} else {
				res += 'T';
			}
		} else {
			if (v[i + 1]) {
				res += 'C';
			} else {
				res += 'A';
			}
		}
	}
	return res;
}

zstr::ifstream* openFile(const string& input_file){
	zstr::ifstream* input_stream = new zstr::ifstream(input_file);
	if(not input_stream-> good()){
		cout << "Problem with file opening" << endl;
        return NULL;
	}
	return input_stream;
}

kmer min_k(const kmer& k1, const kmer& k2) {
	if (k1 <= k2) {
		return k1;
	}
	return k2;
}

void updateK(kmer& min, char nuc, uint64_t& k) {
	min <<= 2;
	min += nuc2int(nuc);
	min &= ((kmer)1<<(2*k))-1;
}

uint64_t canonize(uint64_t x, uint64_t n) {
	return min(x, rcbc(x, n));
}

void Biogetline(zstr::ifstream* in,string& result,char type,uint K) {
  string discard;
  result.clear();
  switch(type){
    case 'Q':
      getline(*in,discard);
      getline(*in,result);
      getline(*in,discard);
      getline(*in,discard);
      break;
    case 'A':
      getline(*in,discard);
      char c=in->peek();
      while(c!='>' and c!=EOF){
        getline(*in,discard);
        transform(discard.begin(),discard.end(),discard.begin(),::toupper);
        result+=discard;
        c=in->peek();
      }
      break;
  }
  if(result.size()< K){
    result.clear();
  }
}