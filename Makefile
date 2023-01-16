CXX=g++
CC=gcc



all: HS



#Other compilation flags
CFLAGS_BLIGHT+= -Ofast -flto -march=native -mtune=native -std=c++17 -pipe -lz -fopenmp -msse4

#Needed object files
BLO=utils.o


HS: HS.o $(BLO)
	$(CXX) -o $@ $^ $(CFLAGS_BLIGHT)

HS.o: HammingSubsampler.cpp HammingSubsampler.h
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

utils.o: utils.cpp $(INC) utils.h
	$(CXX) -o $@ -c $< $(CFLAGS_BLIGHT)

clean:
	rm -rf *.o
	


rebuild: clean all
