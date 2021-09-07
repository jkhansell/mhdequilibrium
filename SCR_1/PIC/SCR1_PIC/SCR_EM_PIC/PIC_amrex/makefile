AMREX_LIBRARY_HOME ?= $(AMReX_ROOT2)

LIBDIR := $(AMREX_LIBRARY_HOME)/lib
INCDIR := $(AMREX_LIBRARY_HOME)/include

COMPILE_CPP_FLAGS ?= $(shell awk '/Cflags:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)
        COMPILE_LIB_FLAGS ?= $(shell awk '/Libs:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)
CFLAGS := -I$(INCDIR) $(COMPILE_CPP_FLAGS)
LFLAGS := -L$(LIBDIR) $(COMPILE_LIB_FLAGS)

all:
	mpic++ -o a.out *.cpp $(CFLAGS) $(LFLAGS)
