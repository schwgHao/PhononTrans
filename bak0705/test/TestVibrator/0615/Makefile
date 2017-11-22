
.SUFFIXES: .mod .o .F90 .f90 .F .f .cpp .c

#
CXX = g++48 -g
FC = /opt/software/mpi/openmpi-1.6.5/bin/mpif90 -g
CCFLAGS =-O3  -std=c++11
FFLAGS=-O3 -std=f2003
LDFLAGS=/usr/lib/libstdc++.so.6.0.18
#
#
OBJS= TestVib.o vibrator_c.o CalcPhLead.o CalcPhEMT.o EMsolver.o ioFC.o\
	   writensc.o SurfacePhononGF.o CfuncDef.o FuncUtils.o

all: $(OBJS)
	$(FC) -o a.out $(LDFLAGS) $(OBJS)

CfuncDef.o:CfuncDef.f90
TestVib.o: vibrator_c.o CfuncDef.o
CalcPhLead.o: ioFC.o SurfacePhononGF.o
vibrator_c.o: CalcPhLead.o CalcPhEMT.o
writensc.o vibrator_c.o: FuncUtils.o
CalcPhEMT.o: ioFC.o EMsolver.o

.cpp.o:
	$(CXX) -c $(CCFLAGS)  $<

.f90.o:
	$(FC) -c $(FFLAGS)  $<

clean: 
	rm *.o *.mod
