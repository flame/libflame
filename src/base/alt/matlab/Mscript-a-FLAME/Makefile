include /u/rvdg/PFLAME/Make.machine

include $(FLAME_ROOT)/Make.include

MATLABPATH = /lusr/share/software/matlab6p5
MATLABINC = $(MATLABPATH)/extern/include
MATLABBIN = $(MATLABPATH)/bin

MEXFILES = \
	FLA_Axpy.mexglx \
	FLA_Copy.mexglx \
	FLA_Copy_t.mexglx \
	FLA_Scal.mexglx \
	FLA_Inv_scal.mexglx \
	FLA_Swap.mexglx \
	FLA_Iamax.mexglx \
	FLA_Nrm2.mexglx \
	FLA_Dot.mexglx \
	FLA_Dot_s.mexglx \
\
	FLA_Gemv.mexglx \
	FLA_Symv.mexglx \
	FLA_Trmv.mexglx \
	FLA_Trsv.mexglx \
	FLA_Ger.mexglx \
	FLA_Syr.mexglx \
	FLA_Syr2.mexglx \
\
	FLA_Gemm.mexglx \
	FLA_Symm.mexglx \
	FLA_Syrk.mexglx \
	FLA_Syr2k.mexglx \
	FLA_Trmm.mexglx \
	FLA_Trsm.mexglx

UTILS = FLA_Matlab2C.o

all: $(MEXFILES)
	rm *.c *.o

FLA_Matlab2C.o: FLA_Matlab2C.c
	$(CC) $(CFLAGS) -I$(MATLABINC) -c FLA_Matlab2C.c

%.mexglx: %.c $(UTILS)
	$(MATLABBIN)/mex -f $(MATLABBIN)/gccopts.sh \
                -I$(FLAME_ROOT)/INCLUDE \
		$< $(UTILS) $(LIB)

%.c:
	cp old/FLA_Matlab2C.c .  
	cat Create_BLAS | ./interface.pl

clean:
	rm -f *.o *.mexglx *~ *.c
