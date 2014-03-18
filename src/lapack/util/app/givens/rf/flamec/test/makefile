#
# test directory makefile
#

FNAME        := Apply_G

SRC_PATH     := ..
OBJ_PATH     := .

INST_PATH    := $(HOME)/flame
LIB_PATH     := $(INST_PATH)/lib
INC_PATH     := $(INST_PATH)/include
FLAME_LIB    := $(LIB_PATH)/libflame.a
#BLAS_LIB     := $(LIB_PATH)/libblas.a
#BLAS_LIB     := $(LIB_PATH)/libgoto.a
BLAS_LIB     := $(LIB_PATH)/liblapack.a $(LIB_PATH)/libgoto.a
#BLAS_LIB     := -L/opt/intel/mkl/10.2.2.025/lib/em64t/ -lmkl_sequential -lmkl_core -lmkl_intel_lp64
#BLAS_LIB     := -L/opt/intel/mkl/10.2.2.025/lib/em64t/ -lmkl_sequential -lmkl_core -lmkl_intel_lp64

CC           := gcc
FC           := gfortran
LINKER       := $(CC)
#CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O3 -Wall -Wno-comment -march=nocona -msse3 -std=c99 -g
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O2 -m64 -Wall -Wno-comment -march=core2 -msse3 -std=c99
FFLAGS       := $(CFLAGS)
AFLAGS       := -S $(CFLAGS)
LDFLAGS      := -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lgfortranbegin -lgfortran -lm
LDFLAGS      += -lpthread -fopenmp

TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))
TEST_ASMS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.s, $(wildcard $(SRC_PATH)/*.c))
TEST_OBJS    += test_$(FNAME).o \
                time_$(FNAME).o
TEST_BIN     := test_$(FNAME).x

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_PATH)/%.s: $(SRC_PATH)/%.c
	$(CC) $(AFLAGS) -c $< -o $@

all: test_$(FNAME) $(TEST_ASMS)

test_$(FNAME): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME_LIB) $(LAPACK_LIB) $(BLAS_LIB) -o $(TEST_BIN)

clean:
	rm -f *.o $(OBJ_PATH)/*.o *~ core *.x *.s

