#-*- Makefile -*-
## Choose compiler: gfortran
COMPILER=gfortran
FC=$(COMPILER)
##Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=lhapdf
#Choose Analysis: dummy, process specific
## default analysis may require FASTJET package, that has to be installed separately (see below)
ANALYSIS=dummy
#ANALYSIS=default
## For static linking uncomment the following
#STATIC= -static
#

OBJ=obj-$(COMPILER)
OBJDIR:=$(OBJ)

ifeq ("$(COMPILER)","gfortran")	
F77=gfortran -fno-automatic -J$(OBJ) -I$(OBJ)
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
#FFLAGS= -Wall -Wimplicit-interface -fbounds-check
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O2
#OPT=-O0
## For debugging uncomment the following
#DEBUG= -ggdb 
ifdef DEBUG
OPT=-O0
FPE=-ffpe-trap=invalid,zero,overflow,underflow
endif
endif

ifeq ("$(COMPILER)","ifort")
F77 = ifort -save  -extend_source  -module $(OBJ)
#CXX = g++
#LIBS = -limf
#FFLAGS =  -checkm
## For floating point exception trapping  uncomment the following 
#FPE = -fpe0
OPT = -O3 #-fast
## For debugging uncomment the following
#DEBUG= -debug -g
ifdef DEBUG
OPT=-O0 
FPE = -fpe0
endif
endif

PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))
VPATH= ./:../:./DHELAS:./Madlib/:./bvirtuals:$(OBJDIR)

INCLUDE0=$(PWD)
INCLUDE1=$(PWD)/include
INCLUDE2=$(shell dirname $(PWD))/include
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE1) -I$(INCLUDE2)

LIBFILES=$(shell  for dir in ./Madlib ./MODEL ./DHELAS ; do cd $$dir ; echo *.[fF] ' ' | sed 's/[fF] /o /g' ; cd .. ; done  )


ifeq ("$(PDF)","lhapdf")
# insert your local path to LHAPDF here:
LHAPDF_CONFIG=/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lhapdf/6.1.6/bin/lhapdf-config
PDFPACK=lhapdf6if.o lhapdf6ifcc.o
FJCXXFLAGS+= $(shell $(LHAPDF_CONFIG) --cxxflags)
LIBSLHAPDF=-Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF -lstdc++
ifeq  ("$(STATIC)","-static") 
## If LHAPDF has been compiled with gfortran and you want to link it statically, you have to include
## libgfortran as well. The same holds for libstdc++. 
## One possible solution is to use fastjet, since $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
## does perform this inclusion. The path has to be set by the user. 
# LIBGFORTRANPATH= #/usr/lib/gcc/x86_64-redhat-linux/4.1.2
# LIBSTDCPP=/lib64
LIBSLHAPDF+=-L$(LIBGFORTRANPATH) -lgfortranbegin -lgfortran -L$(LIBSTDCPP) -lstdc++
endif
LIBS+=$(LIBSLHAPDF)
else
PDFPACK=mlmpdfif.o hvqpdfpho.o
endif


#ifeq ("$(ANALYSIS)","YOURPROCESS")
ifeq ("$(ANALYSIS)","default")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET+= $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o observables.o pwhg_analysis.o
PWHGANAL+=  fastjetfortran.o
else
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-dummy.o
endif


%.o: %.f $(INCLUDE) | $(OBJDIR)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.F $(INCLUDE) | $(OBJDIR)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.f90 $(INCLUDE) | $(OBJDIR)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.c | $(OBJDIR)
	$(CC) $(DEBUG) -c -o $(OBJ)/$@ $^ 

%.o: %.cc | $(OBJDIR)
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^ $(FJCXXFLAGS)


# genuine bbH code:
LIBS+=-lz
USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o \
     bbh_virt_interface.o real.o  $(PWHGANAL)


PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o \
	setstrongcoupl.o pwhg_analysis_driver.o checkmomzero.o   \
	integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	reshufflemoms.o                                                 \
	sigcollsoft.o sigvirtual.o  ubprojections.o                     \
	pwhgreweight.o setlocalscales.o mint_upb.o opencount.o          \
	validflav.o $(PDFPACK) $(USER) $(FPEOBJ) lhefread.o pwhg_io_interface.o rwl_weightlists.o rwl_setup_param_weights.o

LIBDIRMG=$(OBJ)
LINKMGLIBS=-L$(LIBDIRMG) -lmadgraph -lmodel -ldhelas3 -lvirtb

MADLIBS=libdhelas3.a libmadgraph.a libmodel.a libvirtb.a 

# target to generate LHEF output
pwhg_main: $(PWHG) $(MADLIBS)
	$(FF) $(patsubst %,$(OBJ)/%,$(PWHG)) $(LINKMGLIBS) $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@


LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o		\
     opencount.o powheginput.o printleshouches.o $(PWHGANAL)	\
     lhefread.o pwhg_io_interface.o rwl_weightlists.o newunit.o pwhg_analysis_driver.o $(FPEOBJ)

# target to analyze LHEF output
lhef_analysis:$(LHEF)
	$(FF) $(patsubst %,$(OBJ)/%,$(LHEF)) $(LIBS) $(LIBSFASTJET) $(STATIC)  -o $@ 

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA.o pythia.o boostrot.o powheginput.o \
	$(PWHGANAL) lhefread.o pwhg_io_interface.o rwl_weightlists.o newunit.o printleshouches.o		\
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o	\
	$(FPEOBJ)

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA)) $(LIBS) $(LIBSFASTJET)  $(STATIC) -o $@

# target to cleanup
.PHONY: clean libdhelas3.a libmadgraph.a libmodel.a libvirtb.a 

XFFLAGS=$(DEBUG) $(OPT)
libdhelas3.a:
	cd DHELAS ; make FC="$(F77)" F77="$(F77)" XFFLAGS="$(XFFLAGS)" OBJ="$(OBJ)"

libmadgraph.a:
	cd Madlib ; make FC="$(F77)" F77="$(F77)" XFFLAGS"=$(XFFLAGS)" OBJ="$(OBJ)"


libvirtb.a:
	cd bvirtuals ; make FC="$(F77)" F77="$(F77)" XFFLAGS"=$(XFFLAGS)" OBJ="$(OBJ)"


ifeq ("$(COMPILER)","gfortran")
XFFLAGS +=-ffixed-line-length-132
else
XFFLAGS +=-extend-source
endif

libmodel.a:
	cd MODEL ; make FC="$(F77)" F77="$(F77)" XFFLAGS="$(XFFLAGS)" OBJ="$(OBJ)"


clean:
	rm -f $(patsubst %,$(OBJ)/%,$(USER) $(PWHG) $(LHEF) $(PYTHIA)) pwhg_main lhef_analysis	\
	main-PYTHIA-lhef


veryclean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod $(OBJ)/*.a pwhg_main lhef_analysis 	\
	main-PYTHIA-lhef *.a DHELAS/*.o Madlib/*.o MODEL/*.o bvirtuals/*.o 


# target to generate object directory if it does not exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

##########################################################################

