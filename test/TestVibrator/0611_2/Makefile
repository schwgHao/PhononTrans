
.SUFFIXES: .mod .o .F90 .f90 .F .f .cpp .c

#
CXX = g++48
CCFLAGS = -std=c++11
FFLAGS=-O3
LDFLAGS=$(FFLAGS) /usr/lib/libstdc++.so.6.0.18
#
#
OBJS= TestVib.o vibrator_c.o ioFC.o FuncUtils.o SurfacePhononGF.o CalcPhLead.o

all: $(OBJS)
	$(CXX) -o a.out $(LDFLAGS) $(OBJS)

TestVib.o : vibrator_c.o
vibrator_c.o: CalcPhLead.o FuncUtils.o
CalcPhLead.o: ioFC.o SurfacePhononGF.o


.cpp.o:
	$(CXX) -c $(CCFLAGS)  $<

clean: 
	rm *.o *.mod
