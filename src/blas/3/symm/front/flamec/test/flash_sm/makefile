#
# test directory makefile
#

FNAME        := Symm

SRC_PATH     := ../..
OBJ_PATH     := .

INST_PATH    := $(HOME)/flame
LIB_PATH     := $(INST_PATH)/lib
INC_PATH     := $(INST_PATH)/include
FLAME_LIB    := $(LIB_PATH)/libflame.a
#BLAS_LIB     := $(LIB_PATH)/libblas.a
BLAS_LIB     := $(LIB_PATH)/libgoto.a
#BLAS_LIB     := -L/opt/intel/mkl/10.0.3.020/lib/em64t -lmkl -lguide -lpthread

CC           := gcc
FC           := gfortran
LINKER       := $(CC)
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -g -O3 -Wall
FFLAGS       := $(CFLAGS)
#LDFLAGS      := -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.2.3 -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.2.3/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.2.3/../../../../x86_64-pc-linux-gnu/lib -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.2.3/../../.. -lgfortranbegin -lgfortran -lm -lpthread
LDFLAGS      := -L/lusr/opt/gcc-4.2.2/lib/gcc/i686-pc-linux-gnu/4.2.2 -L/lusr/opt/gcc-4.2.2/lib/gcc/i686-pc-linux-gnu/4.2.2/../../.. -lgfortranbegin -lgfortran -lm -lpthread

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
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME_LIB) $(BLAS_LIB) -o $(TEST_BIN)

clean:
	rm -f *.o $(OBJ_PATH)/*.o *~ core *.x

