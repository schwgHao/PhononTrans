
.SUFFIXES: .mod .o .F90 .f90 .F .f .cpp .c

#
CXX = g++48
CCFLAGS = -std=c++11
FFLAGS=-O3
LDFLAGS=$(FFLAGS) /usr/lib/libstdc++.so.6.0.18
#
#
OBJS= TestVib.o vibrator_c.o ioFC.o FuncUtils.o

all: $(OBJS)
	$(CXX) -o a.out $(LDFLAGS) $(OBJS)

TestVib.o : vibrator_c.o
vibrator_c.o: ioFC.o
vibrator_c.o: FuncUtils.o




.cpp.o:
	$(CXX) -c $(CCFLAGS)  $<

clean: 
	rm *.o *.mod
