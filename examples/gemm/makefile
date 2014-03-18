BLAS_LIB   := $(HOME)/flame/lib/libgoto.a
#BLAS_LIB   := $(HOME)/flame/lib/libblas.a
FLAME_LIB  := $(HOME)/flame/lib/libflame.a
FLAME_INC  := $(HOME)/flame/include

# indicate where the object files are to be created
CC         := gcc
FC         := gfortran
LINKER     := $(CC)
CFLAGS     := -O3 -Wall -I$(FLAME_INC)
FFLAGS     := $(CFLAGS)

LDFLAGS    := -lm -lpthread -L/lusr/opt/gcc-4.2.2/lib/gcc/i686-pc-linux-gnu/4.2.2 -L/lusr/opt/gcc-4.2.2/lib/gcc/i686-pc-linux-gnu/4.2.2/../../.. -lgfortranbegin -lgfortran -lm

# indicate where the FLAME include files reside

TEST_OBJS  := driver.o \
	Gemm_unb_var1.o Gemm_blk_var1.o \
	Gemm_unb_var2.o Gemm_blk_var2.o \
	Gemm_unb_var3.o Gemm_blk_var3.o \
	Gepp_blk_var1.o Gebp.o

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

driver.x: $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME_LIB) $(BLAS_LIB) -o driver.x

test:   driver.x
	echo "3 6.8 128 50 500 50" | ./driver.x 
clean:
	rm -f *.o *~ core *.x

