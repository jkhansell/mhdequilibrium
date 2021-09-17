CC=g++
COMMON_FLAGS=-O3 -std=c++0x -o bs-solctra-seq -Wall
SOURCE=solctra_sequential.h solctra_sequential.cpp main_sequential.cpp utils.h utils.cpp
RPT_FLAGS=
FP_FLAGS=


gcc: $(SOURCE)
	g++ -O3 -std=c++11 -o bs-solctra-seq $(SOURCE)
	#g++ -O0 -g -std=c++11 -o bs-solctra-seq $(SOURCE) 
	cp bs-solctra-seq resultados/;
	rm bs-solctra-seq;

debug: $(SOURCE)
	g++ -O3 -std=c++11 -pg -o bs-solctra-seq $(SOURCE)
	cp bs-solctra-seq resultados/;
	rm bs-solctra-seq;

clean:
	rm -f bs-solctra-seq;
	rm resultados/bs-solctra-seq;
