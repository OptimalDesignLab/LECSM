# makefile for Kona library

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o
.PHONY: default call tags test_problem clean

# compiler
CXX= g++

# c preprocessor options
CPPFLAGS= -cpp -DDEBUG -DBOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
#CPPFLAGS= -cpp -DNDEBUG -DBOOST_UBLAS_ENABLE_PROXY_SHORTCUTS -DBOOST_UBLAS_TYPE_CHECK=0

# compiler options that may vary (user can change)
CXXFLAGS= -g 
BOOST_ROOT= $(BOOST_HOME)
KONA_ROOT= $(KONA_HOME)

# linker options
LDFLAGS= -lstdc++ -llapack -L. -lboost_program_options -lm -L$(KONA_ROOT) -lkona

# options that DO NOT vary
ALL_CXXFLAGS= -I. $(CXXFLAGS) -I$(BOOST_ROOT) -I$(KONA_ROOT)/src

# directories
SRC=../

# source and object file names
HEADERS= $(wildcard $(SRC)*.hpp *.hpp)
HEADERS_ALL= $(HEADERS)
SOURCES= $(wildcard $(SRC)*.cpp *cpp)
SOURCES_ALL= $(SOURCES)
OBJS= $(SOURCES:.cpp=.o)
OBJS_ALL= $(SOURCES_ALL:.cpp=.o)
BINARIES= test_problem.bin

# implicit rule
%.o : %.cpp $(HEADERS_ALL) Makefile
	@echo "Compiling \""$@"\" from \""$<"\""
	@$(CXX) $(CPPFLAGS) $(ALL_CXXFLAGS) -o $@ -c $<

default: all

all: $(BINARIES)

tags: $(HEADERS) $(SOURCES)
	@echo "Creating TAGS file for emacs"
	@find -maxdepth 2 -iname '*.hpp' -print0 -o \
	-iname '*.cpp' -print0 | xargs -0 etags

test_problem.bin: $(OBJS) Makefile
	@echo "Compiling \""$@"\" from \""$(OBJS)"\""
	@$(CXX) -o $@ $(OBJS) $(LDFLAGS)

clean:
	@echo "deleting temporary, object, and binary files"
	@rm -f *~
	@rm -f $(BINARIES)
	@rm -f $(OBJS) *.o
