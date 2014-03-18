#
# test directory makefile
#

FNAME      := query_blocksizes

SRC_PATH   := .
OBJ_PATH   := .

LIB_PATH   := $(HOME)/flame/lib
INC_PATH   := $(HOME)/flame/include 
BLAS_LIB   := $(LIB_PATH)/libblas_clover.a
FLAME_BLAS := $(LIB_PATH)/libflame.a

CC         := gcc
FC         := gfortran
LINKER     := $(CC)
CFLAGS     := -I$(SRC_PATH) -I$(INC_PATH)
FFLAGS     := $(CFLAGS)
LDFLAGS    := -lm -lpthread -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.2.3 -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.2.3/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.2.3/../../../../x86_64-pc-linux-gnu/lib -L/usr/lib/gcc/x86_64-pc-linux-gnu/4.2.3/../../.. -lgfortranbegin -lgfortran -lm

TEST_OBJS  += $(OBJ_PATH)/$(FNAME).o
TEST_BIN   := $(FNAME).x

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJ_PATH)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME_BLAS) $(FLAME_BASE) $(BLAS_LIB) -o $(TEST_BIN)

clean:
	rm -f *.o *~ core *.x

