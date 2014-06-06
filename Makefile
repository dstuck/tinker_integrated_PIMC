#CXX=icpc
#CXX=g++
#CXXFLAGS=-O0 -g $(INCLUDES)
CXXFLAGS=-O0 -g -DDEBUG
F77FLAGS=-cxxlib -nofor-main
CSOURCES := $(wildcard *.cpp)
COBJECTS := $(patsubst %.cpp,%.o,$(CSOURCES))
CDEPFILES := $(patsubst %.cpp,%.d,$(CSOURCES))
#LIBTINKER := -L tinker/libtinker.a
LIBTINKER := tinker/active.o tinker/analysis.o tinker/angles.o tinker/attach.o tinker/basefile.o tinker/beeman.o tinker/bicubic.o tinker/bitors.o tinker/bonds.o tinker/born.o tinker/bounds.o tinker/bussi.o tinker/calendar.o tinker/center.o tinker/chkpole.o tinker/chkring.o tinker/chkxyz.o tinker/cholesky.o tinker/clock.o tinker/cluster.o tinker/column.o tinker/command.o tinker/connect.o tinker/connolly.o tinker/control.o tinker/cspline.o tinker/cutoffs.o tinker/deflate.o tinker/delete.o tinker/diagq.o tinker/diffeq.o tinker/eangang.o tinker/eangang1.o tinker/eangang2.o tinker/eangang3.o tinker/eangle.o tinker/eangle1.o tinker/eangle2.o tinker/eangle3.o tinker/ebond.o tinker/ebond1.o tinker/ebond2.o tinker/ebond3.o tinker/ebuck.o tinker/ebuck1.o tinker/ebuck2.o tinker/ebuck3.o tinker/echarge.o tinker/echarge1.o tinker/echarge2.o tinker/echarge3.o tinker/echgdpl.o tinker/echgdpl1.o tinker/echgdpl2.o tinker/echgdpl3.o tinker/edipole.o tinker/edipole1.o tinker/edipole2.o tinker/edipole3.o tinker/egauss.o tinker/egauss1.o tinker/egauss2.o tinker/egauss3.o tinker/egeom.o tinker/egeom1.o tinker/egeom2.o tinker/egeom3.o tinker/ehal.o tinker/ehal1.o tinker/ehal2.o tinker/ehal3.o tinker/eimprop.o tinker/eimprop1.o tinker/eimprop2.o tinker/eimprop3.o tinker/eimptor.o tinker/eimptor1.o tinker/eimptor2.o tinker/eimptor3.o tinker/elj.o tinker/elj1.o tinker/elj2.o tinker/elj3.o tinker/embed.o tinker/emetal.o tinker/emetal1.o tinker/emetal2.o tinker/emetal3.o tinker/emm3hb.o tinker/emm3hb1.o tinker/emm3hb2.o tinker/emm3hb3.o tinker/empole.o tinker/empole1.o tinker/empole2.o tinker/empole3.o tinker/energy.o tinker/dstuckVibrate.o tinker/dstuckEnergy.o tinker/eopbend.o tinker/eopbend1.o tinker/eopbend2.o tinker/eopbend3.o tinker/eopdist.o tinker/eopdist1.o tinker/eopdist2.o tinker/eopdist3.o tinker/epitors.o tinker/epitors1.o tinker/epitors2.o tinker/epitors3.o tinker/erf.o tinker/erxnfld.o tinker/erxnfld1.o tinker/erxnfld2.o tinker/erxnfld3.o tinker/esolv.o tinker/esolv1.o tinker/esolv2.o tinker/esolv3.o tinker/estrbnd.o tinker/estrbnd1.o tinker/estrbnd2.o tinker/estrbnd3.o tinker/estrtor.o tinker/estrtor1.o tinker/estrtor2.o tinker/estrtor3.o tinker/etors.o tinker/etors1.o tinker/etors2.o tinker/etors3.o tinker/etortor.o tinker/etortor1.o tinker/etortor2.o tinker/etortor3.o tinker/eurey.o tinker/eurey1.o tinker/eurey2.o tinker/eurey3.o tinker/evcorr.o tinker/extra.o tinker/extra1.o tinker/extra2.o tinker/extra3.o tinker/fatal.o tinker/fft3d.o tinker/fftpack.o tinker/field.o tinker/final.o tinker/flatten.o tinker/freeunit.o tinker/geometry.o tinker/getint.o tinker/getkey.o tinker/getmol.o tinker/getmol2.o tinker/getnumb.o tinker/getpdb.o tinker/getprm.o tinker/getref.o tinker/getstring.o tinker/gettext.o tinker/getword.o tinker/getxyz.o tinker/ghmcstep.o tinker/gradient.o tinker/gradrgd.o tinker/gradrot.o tinker/groups.o tinker/grpline.o tinker/gyrate.o tinker/hessian.o tinker/hessrgd.o tinker/hessrot.o tinker/hybrid.o tinker/image.o tinker/impose.o tinker/induce.o tinker/inertia.o tinker/initatom.o tinker/initial.o tinker/initprm.o tinker/initres.o tinker/initrot.o tinker/insert.o tinker/invbeta.o tinker/invert.o tinker/jacobi.o tinker/kangang.o tinker/kangle.o tinker/katom.o tinker/kbond.o tinker/kcharge.o tinker/kdipole.o tinker/kewald.o tinker/kgeom.o tinker/kimprop.o tinker/kimptor.o tinker/kinetic.o tinker/kmetal.o tinker/kmpole.o tinker/kopbend.o tinker/kopdist.o tinker/korbit.o tinker/kpitors.o tinker/kpolar.o tinker/ksolv.o tinker/kstrbnd.o tinker/kstrtor.o tinker/ktors.o tinker/ktortor.o tinker/kurey.o tinker/kvdw.o tinker/lattice.o tinker/lbfgs.o tinker/lights.o tinker/makeint.o tinker/makeref.o tinker/makexyz.o tinker/maxwell.o tinker/mdinit.o tinker/mdrest.o tinker/mdsave.o tinker/mdstat.o tinker/mechanic.o tinker/merge.o tinker/molecule.o tinker/moments.o tinker/mutate.o tinker/nblist.o tinker/nextarg.o tinker/nexttext.o tinker/nose.o tinker/nspline.o tinker/number.o tinker/numeral.o tinker/numgrad.o tinker/ocvm.o tinker/openend.o tinker/optsave.o tinker/orbital.o tinker/orient.o tinker/orthog.o tinker/overlap.o tinker/picalc.o tinker/pmestuff.o tinker/pmpb.o tinker/polymer.o tinker/precise.o tinker/pressure.o tinker/prmkey.o tinker/promo.o tinker/prtdyn.o tinker/prterr.o tinker/prtint.o tinker/prtmol2.o tinker/prtpdb.o tinker/prtprm.o tinker/prtseq.o tinker/prtxyz.o tinker/quatfit.o tinker/random.o tinker/rattle.o tinker/readdyn.o tinker/readgau.o tinker/readint.o tinker/readmol.o tinker/readmol2.o tinker/readpdb.o tinker/readprm.o tinker/readseq.o tinker/readxyz.o tinker/replica.o tinker/respa.o tinker/rgdstep.o tinker/rings.o tinker/rmsfit.o tinker/rotlist.o tinker/rotpole.o tinker/sdstep.o tinker/search.o tinker/server.o tinker/shakeup.o tinker/sigmoid.o tinker/sktstuff.o tinker/sort.o tinker/square.o tinker/suffix.o tinker/surface.o tinker/surfatom.o tinker/switch.o tinker/temper.o tinker/tncg.o tinker/torphase.o tinker/torque.o tinker/torsions.o tinker/trimtext.o tinker/unitcell.o tinker/verlet.o tinker/version.o tinker/volume.o tinker/xyzatm.o tinker/zatom.o

# LDFLAGS=-L/opt/intel/Compiler/11.1/080/Frameworks/mkl/lib/em64t/ -lblas -L/usr/lib/libshell/ -lshell
.PHONY: clean
all: pimcTinker 

# smanzer: Need these soon aomp2.h sort_tuple_list.o  
pimcTinker: $(COBJECTS) $(CDEPFILES)
#	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@	
#	$(CXX) $(COBJECTS) $(LIBTINKER) -o $@	
	$(F77) $(F77FLAGS) $(LIBTINKER) $(COBJECTS) -o $@	
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $^
clean: 
	-rm pimcTinker
	-rm *.o
	-rm *.d

$(CDEPFILES): %.d: %.cpp
	set -e; rm -f $@; \
	$(CXX) -M $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
