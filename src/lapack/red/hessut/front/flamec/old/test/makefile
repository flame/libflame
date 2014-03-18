#
# test directory makefile
#

FNAME        := Hess

SRC_PATH     := ..
OBJ_PATH     := .

LIB_PATH     := $(HOME)/flame/lib
INC_PATH     := $(HOME)/flame/include
FLAME_BASE   := $(LIB_PATH)/libflame-base_ia64.a
FLAME_LAPACK := $(LIB_PATH)/libflame-lapack_ia64.a
LAPACK2FLAME   := $(LIB_PATH)/liblapack2flame_ia64.a
BLAS_LIB     := $(LIB_PATH)/libgoto_ia64-it2.a
#BLAS_LIB     := $(LIB_PATH)/libgoto_ia64-it2p.a
LAPACK_LIB   := $(LIB_PATH)/liblapack_ia64-v3.0.a
#LAPACK_LIB   := $(LIB_PATH)/liblapack_ia64-v3.1.0.a
#BLAS_LIB     := -L/opt/intel/mkl/9.0/lib/64 -lmkl
#LAPACK_LIB   := -L/opt/intel/mkl/9.0/lib/64 -lmkl_lapack -lmkl_ipf -lguide -lmkl

CC           := icc
FC           := ifort
LINKER       := $(CC)
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -g #-O3
FFLAGS       := $(CFLAGS)
LDFLAGS      := -lm -lpthread -L/opt/intel/fc/ia64/9.1.036/lib -L/usr/lib/gcc-lib/ia64-unknown-linux-gnu/3.3.2/ -L/usr/lib/gcc-lib/ia64-unknown-linux-gnu/3.3.2/../../../ -lifport -lifcore -limf -lm -lipgo -lirc -lunwind -lgcc_s -lirc_s -ldl -openmp

TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))
TEST_OBJS    += test_$(FNAME).o \
                time_$(FNAME).o \
                REF_$(FNAME).o
TEST_BIN     := test_$(FNAME).x

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME_LAPACK) $(FLAME_BASE) $(LAPACK_LIB) $(BLAS_LIB) -o $(TEST_BIN)

clean:
	rm -f *.o $(OBJ_PATH)/*.o *~ core *.x

