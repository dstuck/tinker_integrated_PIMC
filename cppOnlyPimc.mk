CXX=icpc
#CXX=g++
#CXXFLAGS=-O0 -g $(INCLUDES)
CXXFLAGS=-O0 -g
CSOURCES := $(wildcard *.cpp)
COBJECTS := $(patsubst %.cpp,%.o,$(CSOURCES))
CDEPFILES := $(patsubst %.cpp,%.d,$(CSOURCES))

# LDFLAGS=-L/opt/intel/Compiler/11.1/080/Frameworks/mkl/lib/em64t/ -lblas -L/usr/lib/libshell/ -lshell
.PHONY: clean
all: pimcTinker 

# smanzer: Need these soon aomp2.h sort_tuple_list.o  
pimcTinker: $(COBJECTS) $(CDEPFILES)
#	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@	
	$(CXX) $(COBJECTS) -o $@	
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $^
clean: 
	-rm pimcTinker
	-rm *.o

$(CDEPFILES): %.d: %.cpp
	set -e; rm -f $@; \
	$(CXX) -M $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
