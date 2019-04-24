#CXX=icpc
#CXX=g++
#CXXFLAGS=-O0 -g $(INCLUDES)
ARMADILLO_BASE := $(if $(ARMADILLO_BASE),$(ARMADILLO_BASE),'/usr/local/Cellar/armadillo/9.300.2')
ARMAFLAGS=-L$(ARMADILLO_BASE) -I$(ARMADILLO_BASE)/include -larmadillo -framework Accelerate
CXXFLAGS=-O0 -g -DDEBUG $(ARMAFLAGS)
F77FLAGS=-cxxlib -nofor-main -Xlinker -no_compact_unwind -framework Accelerate
#CSOURCES := $(wildcard *.cpp)
CSOURCES := $(wildcard *.cpp) $(wildcard optking/*.cpp)
COBJECTS := $(patsubst %.cpp,%.o,$(CSOURCES))
CDEPFILES := $(patsubst %.cpp,%.d,$(CSOURCES))
#LIBTINKER := -L tinker/libtinker.a
LIBTINKER :=  $(wildcard tinker/*.o)
LIBPOTLIB :=  $(wildcard potlib/*.o)

# LDFLAGS=-L/opt/intel/Compiler/11.1/080/Frameworks/mkl/lib/em64t/ -lblas -L/usr/lib/libshell/ -lshell
.PHONY: clean
all: pimcTinker 

# smanzer: Need these soon aomp2.h sort_tuple_list.o  
pimcTinker: $(COBJECTS) $(CDEPFILES)
#	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@	
#	$(CXX) $(COBJECTS) $(LIBTINKER) -o $@	
	$(F77) $(F77FLAGS) $(LIBTINKER) $(LIBPOTLIB) $(COBJECTS) -o $@	
%.o: %.cpp
#	$(CXX) $(CXXFLAGS) -c $^
	$(CXX) $(CXXFLAGS) -c -o $@ $<
clean: 
	-rm pimcTinker
	-rm optking/*.o
	-rm optking/*.d
	-rm *.o
	-rm *.d

$(CDEPFILES): %.d: %.cpp
	set -e; rm -f $@; \
	$(CXX) -M $(ARMAFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
