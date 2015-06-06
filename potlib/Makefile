#F77 = /opt/intel/Compiler/11.1/080/bin/intel64/ifort
LIBS =
F77FLAGS = -c -axSSSE3 -vec-report0
OPTFLAGS = -O3 -no-ipo -no-prec-div
#OPTFLAGS = -O0 -no-ipo -no-prec-div
LIBFLAGS = -crusv
LINKFLAGS = -static-intel -mmacosx-version-min=10.4

#################################################################
##  Should not be Necessary to Change Things Below this Point  ##
#################################################################

OBJS = h2o2.o

.f.o:
	${F77} ${F77FLAGS} ${OPTFLAGS} $*.f 

#pmpb.o:
#	${F77} ${F77FLAGS} ${OPTFLAGS} pmpb.f

#server.o:
#	${F77} ${F77FLAGS} ${OPTFLAGS} server.f

#%.x: %.o libtinker.a
#	${F77} ${LINKFLAGS} -o $@ $^ ${LIBS}; strip $@

#all:	${EXEFILES}
all:	${OBJS}

clean:
	rm -f *.o *.a *.x
