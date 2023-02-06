# HS
Subsampling strategy based on Hamming codes. This is a proof of concept.

## Install
```sh
git clone https://github.com/TimRouze/HS.git
```
```sh
cd HS
```
## Usage
```sh
make
```
```sh
./HS -i my_fasta_file.fa(.gz) [-r 5] [-t A] [-o "Hammed"]
```
## Options
-i Input file (MANDATORY)  
-o Output prefix (Base value: "Hammed")  
-r R value, k = (2^r)/2 (Base value: 5 (k = 16))  
-t File type, A for fasta. Q for fastq. (Base value is A).  
      Reads gz files too!!

## Authors
Gautier Couture, Timothé Rouzé, Mikael Salson.
