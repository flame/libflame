#
# test directory makefile
#

FNAME        := Gemm_nn

SRC_PATH     := ../..
OBJ_PATH     := ..

LIB_PATH     := $(HOME)/flame/lib
INC_PATH     := $(HOME)/flame/include 
BLAS_LIB     := $(LIB_PATH)/libgoto.a
FLAME_BASE   := $(LIB_PATH)/libflame-base.a
FLAME_BLAS   := $(LIB_PATH)/libflame-blas.a

CC           := icc
FC           := ifort
LINKER       := $(CC)
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O3 -openmp
FFLAGS       := $(CFLAGS)
LDFLAGS      := -lm -openmp

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
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(BLAS_LIB) $(FLAME_BLAS) $(FLAME_BASE) -o $(TEST_BIN)

clean:
	rm -f *.o $(OBJ_PATH)/*.o $(SRC_PATH)/*.o *~ core *.x

