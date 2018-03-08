
.SUFFIXES: .mod .o .F90 .f90 .F .f .cpp .c

#
SRCpath=../../Src
include $(SRCpath)/arch.make
#
#
OBJS=vibrator_c.o CalcPhLead.o CalcPhEMT.o EMsolver.o ioFC.o\
     writensc.o SurfacePhononGF.o CfuncDef.o FuncUtils.o

module: $(OBJS)
	cp $(OBJS) libbindc.mod $(SRCpath)
	@echo "phonon Transport Interfaces"

CfuncDef.o: CfuncDef.f90
CalcPhLead.o: ioFC.o SurfacePhononGF.o
vibrator_c.o: CalcPhLead.o CalcPhEMT.o
writensc.o vibrator_c.o: FuncUtils.o
CalcPhEMT.o: ioFC.o EMsolver.o

clean: 
	rm *.o *.mod
