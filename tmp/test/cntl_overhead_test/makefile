#
# test directory makefile
#

FNAME      := Chol

SRC_PATH   := ../..
OBJ_PATH   := .

LIB_PATH   := $(HOME)/flame/lib
INC_PATH   := $(HOME)/flame/include
FLAME_BASE := $(LIB_PATH)/libflame-base_x86_64.a
FLAME_BLAS := $(LIB_PATH)/libflame-blas_x86_64.a
FLAME_LAPACK := $(LIB_PATH)/libflame-lapack_x86_64.a
BLAS_LIB   := $(LIB_PATH)/libgoto_x86_64.a
LAPACK_LIB := $(LIB_PATH)/liblapack_x86_64-v3.0.a

CC         := gcc
FC         := gfortran
LINKER     := $(CC)
CFLAGS     := -I$(SRC_PATH) -I$(INC_PATH) -g -O3 -Wall
FFLAGS     := $(CFLAGS)
LDFLAGS    := -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.1.2 -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.1.2/../../../../x86_64-pc-linux-gnu/lib -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.1.2/../../.. -lgfortranbegin -lgfortran -lm -lpthread

TEST_OBJS  := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))
TEST_OBJS  += test_$(FNAME).o \
              time_$(FNAME).o \
              REF_$(FNAME).o \
              FLA_Chol_l_blk_var3.o \
              FLA_Syrk_ln_blk_var2.o \
              FLA_Gemm_nn_blk_var1.o \
              FLA_Trsm_rlt_blk_var3.o
TEST_BIN   := test_$(FNAME).x

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME_LAPACK) $(FLAME_BLAS) $(FLAME_BASE) $(LAPACK_LIB) $(BLAS_LIB) -o $(TEST_BIN)

clean:
	rm -f *.o $(OBJ_PATH)/*.o *~ core *.x

