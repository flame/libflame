#
# test directory makefile
#

FNAME        := Eig_gest

SRC_PATH     := ../..
OBJ_PATH     := .

INST_PATH    := /dockerx/libflame_temp
LIB_PATH     := $(INST_PATH)/lib/x86_64-unknown-linux-gnu/
INC_PATH     := $(INST_PATH)/include/x86_64-unknown-linux-gnu/
FLAME_LIB    := $(LIB_PATH)/libflame.so
#BLAS_LIB     := $(LIB_PATH)/libblas.a
#BLAS_LIB     := $(LIB_PATH)/libgoto.a
BLAS_LIB     := -lopenblas

CC           := gcc
FC           := gfortran
LINKER       := $(CC)
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -I/opt/rocm/include -g -O3 -Wall -Wno-comment
FFLAGS       := $(CFLAGS)
LDFLAGS      := -lm
LDFLAGS      += -lpthread

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
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME_LIB) $(LAPACK_LIB) $(BLAS_LIB) -o $(TEST_BIN)

clean:
	rm -f *.o $(OBJ_PATH)/*.o *~ core *.x

