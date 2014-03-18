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

CC           := gcc
FC           := gfortran
LINKER       := $(CC)
#CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -g -O3 -Wall -Wno-comment
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O3 -Wall -Wno-comment -march=native -msse3 -malign-double -std=c99 -g #-mfpmath=sse
AFLAGS       := -S $(CFLAGS)
FFLAGS       := $(CFLAGS)
LDFLAGS      := -L/lusr/opt/gcc-4.2.2/lib/gcc/i686-pc-linux-gnu/4.2.2 -L/lusr/opt/gcc-4.2.2/lib/gcc/i686-pc-linux-gnu/4.2.2/../../.. -lgfortranbegin -lgfortran -lm
LDFLAGS      += -lpthread

TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))
TEST_ASMS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.s, $(wildcard $(SRC_PATH)/*.c))
TEST_OBJS    += test_$(FNAME).o \
                time_$(FNAME).o \
                FLA_Apply_G_rf_opt_var1.o \
                FLA_Apply_G_rf_asm_var1.o \
                FLA_Apply_G_rf_blk_var2.o \
                FLA_Apply_G_rf_asm_var2.o \
                FLA_Apply_G_rf_blk_var3.o \
                FLA_Apply_G_rf_asm_var3.o \
                FLA_Apply_G_rf_blk_var6.o \
                FLA_Apply_G_rf_asm_var6.o \
                FLA_Apply_G_rf_blk_var9.o \
                FLA_Apply_G_rf_asm_var9.o

TEST_BIN     := test_$(FNAME).x

#$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
#	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_PATH)/%.s: $(SRC_PATH)/%.c
	$(CC) $(AFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_OBJS) $(TEST_ASMS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME_LIB) $(LAPACK_LIB) $(BLAS_LIB) -o $(TEST_BIN)

clean:
	rm -f $(TEST_OBJS) $(TEST_ASMS) *.x

