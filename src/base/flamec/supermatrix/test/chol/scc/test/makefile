#
# test directory makefile
#

RCCEROOT=$(HOME)/rcce/rcce
include $(RCCEROOT)/common/symbols

FNAME        := supermatrix

SRC_PATH     := ..
OBJ_PATH     := .

LIB_PATH     := $(HOME)/flame/lib
INC_PATH     := $(HOME)/flame/include

FLAME        := $(LIB_PATH)/libflame.a
MKL_PATH     := /shared/intel/mkl/8.1.1/lib/32
BLAS         := $(MKL_PATH)/libmkl_lapack.a $(MKL_PATH)/libmkl_ia32.a 

CC           := icpc
LINKER       := $(CC)
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -I$(RCCEROOT)/include -O3
LDFLAGS      := -lm -openmp -static

TEST_BIN     := $(FNAME).x
TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

compile: $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(FLAME) $(LDFLAGS) $(BLAS) $(ARCHIVE) -o $(TEST_BIN)

run:
	$(RCCEROOT)/rccerun -clock 0.533 -nue 4 -f $(RCCEROOT)/hosts/rc.hosts ./$(TEST_BIN) 1 10 100 1 1

clean:
	rm -f *.o *~ core *.x
