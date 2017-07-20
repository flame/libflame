FLAMEPATH = /u/vkv/FLAME1.0
FLAMEINC = $(FLAMEPATH)/INCLUDE
CFLAME = $(FLAMEPATH)/libFLAME1.0.a

include $(FLAMEPATH)/Make.include

ATLASPATH = /p/share/ATLAS/Linux_PIIISSE1
ATLASLIBPATH = $(ATLASPATH)/lib
ATLAS = $(ATLASLIBPATH)/libf77blas.a $(ATLASLIBPATH)/libcblas.a $(ATLASLIBPATH)/libatlas.a

MATLABPATH = /lusr/share/software/matlab6p5
MATLABINC = $(MATLABPATH)/extern/include
MATLABBIN = $(MATLABPATH)/bin

MEXFILES = FLA_Axpy.mexglx \
	FLA_Axpy_t.mexglx \
	FLA_Trmvsx.mexglx

UTILS = FLA_Matlab2C.o

all: $(MEXFILES)

FLA_Matlab2C.o: FLA_Matlab2C.c
	$(CC) $(CFLAGS) -I$(MATLABINC) -c FLA_Matlab2C.c

%.mexglx: %.c $(UTILS)
	$(MATLABBIN)/mex -f $(MATLABBIN)/gccopts.sh -I$(FLAMEINC) \
		$< $(UTILS) $(CFLAME) $(ATLAS)

clean:
	rm *.o *.mexglx
