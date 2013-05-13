# makefile for linear_elastic_csm

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o
.PHONY: default clean all lecsm

# compiler
CXX= gcc 

# c preprocessor options
CPPFLAGS= -cpp

# compiler options that may vary (user can change)
CXXFLAGS= 
BOOST_ROOT= /usr/local/boost_1_47_0

# linker options
# NOTE: -L is for linking, -Wl,-rpath is for loading
LDFLAGS= -lstdc++ -L$(BOOST_ROOT)/stage/lib -lboost_program_options

# options that DO NOT vary
ALL_CXXFLAGS= -I. $(CXXFLAGS) -I $(BOOST_ROOT)

# source, object, and executable file names
HEADERS= $(wildcard *.hpp)
SOURCES= $(wildcard *.cpp)
OBJS= $(SOURCES:.cpp=.o)
BINARIES= lecsm.bin

# implicit rule
%.o : %.cpp $(HEADERS) Makefile
	@echo "Compiling \""$@"\" from \""$<"\""
	@$(CXX) $(CPPFLAGS) $(ALL_CXXFLAGS) -o $@ -c $<

default: all

all: $(BINARIES)

lecsm.bin: $(OBJS) linear_elastic_csm.o Makefile
	@echo "Compiling \""$@"\" from \""$(OBJS)"\""
	@$(CXX) -o $@ linear_elastic_csm.o $(OBJS) $(LDFLAGS) 

clean:
	@echo "deleting temporary, object, and binary files"
	@rm -f *~
	@rm -f $(BINARIES)
	@rm -f $(OBJS) *.o