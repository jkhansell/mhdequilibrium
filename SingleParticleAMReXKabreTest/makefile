AMREX_LIBRARY_HOME ?= /home/johan/installdir
CXX=mpic++

SOURCE=Prob.cpp WritePlotFile.cpp Particles.cpp ParticleMesh.cpp main.cpp

LIBDIR := $(AMREX_LIBRARY_HOME)/lib
INCDIR := $(AMREX_LIBRARY_HOME)/include

all: main

main: $(SOURCE)
	$(CXX) -std=c++14 -fopenmp -I$(INCDIR) -o main $(SOURCE) -L$(LIBDIR) -lamrex -lpthread
#mpic++ -DAMReX_MPI=YES -c -I/home/djimenez/amrex/installdir/include -g *.cpp *.h
#mpic++ -o main -L/home/djimenez/amrex/installdir/lib -lamrex main.o

clean:
	rm main 

#COMPILE_CPP_FLAGS ?= $(shell awk '/Cflags:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)
#COMPILE_LIB_FLAGS ?= $(shell awk '/Libs:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)


