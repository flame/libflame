#
# test directory makefile
#

FNAME        := Hevd_lv

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
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O3 -Wall -Wno-comment -march=native -msse3 -mfpmath=sse -malign-double
FFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O3 -Wall
LDFLAGS      := -L/lusr/opt/gcc-4.2.2/lib/gcc/i686-pc-linux-gnu/4.2.2 -L/lusr/opt/gcc-4.2.2/lib/gcc/i686-pc-linux-gnu/4.2.2/../../.. -lgfortranbegin -lgfortran -lm
LDFLAGS      += -lpthread

TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))
TEST_OBJS    += test_$(FNAME).o \
                time_$(FNAME).o zsteqr.o \
                REF_$(FNAME).o REF_Hevdd_lv.o REF_Hevdr_lv.o zheevd.o zunmqr.o
TEST_BIN     := test_$(FNAME).x

test_$(FNAME): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME_LIB) $(LAPACK_LIB) $(BLAS_LIB) -o $(TEST_BIN)

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o $(OBJ_PATH)/*.o *~ core *.x

