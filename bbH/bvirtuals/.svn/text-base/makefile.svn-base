#FFLAGS = -Mextend
#FC=pgf77
##FFLAGS = -O0 -fno-emulate-complex -ffast-math
##FC=g77
FC=gfortran
FFLAGS        = $(XFFLAGS)
LIBDIR        = ../$(OBJ)
LIBRARY	      = $(LIBDIR)/libvirtb.a


CODE=bbh_virtual.o lib.o intlib.o libqq.o libglu.o matri.o boxlib.o boxlib_qq.o penta_qq.o  p1new.o p2new.o p3new.o p4new.o p5new.o p6new.o alphas2.o


.f.o: ; $(F77) $(FFLAGS) -c $*.f

all: $(LIBDIR)libvirtb.a

$(LIBDIR)libvirtb.a: $(CODE)
	ar cru libvirtb.a $(CODE)
	ranlib libvirtb.a
	mv libvirtb.a  $(LIBDIR)
