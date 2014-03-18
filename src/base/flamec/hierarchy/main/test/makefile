#
# test directory makefile
#

FNAME        := FLASH_Obj

SRC_PATH     := ..
OBJ_PATH     := .

LIB_PATH     := $(HOME)/flame/lib
INC_PATH     := $(HOME)/flame/include
FLAME_BASE   := $(LIB_PATH)/libflame-base.a
FLAME_LAPACK := $(LIB_PATH)/libflame-lapack.a
BLAS_LIB     := $(LIB_PATH)/libgoto.a
LAPACK_LIB   := $(LIB_PATH)/liblapack.a

CC           := gcc
FC           := gfortran
LINKER       := $(CC)
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -g -O1 -Wall
FFLAGS       := $(CFLAGS)
LDFLAGS      := -lm -lpthread

TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))
TEST_OBJS    += $(OBJ_PATH)/test_$(FNAME).o
TEST_BIN     := test_$(FNAME).x

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJ_PATH)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME_BLAS) $(FLAME_BASE) $(BLAS_LIB) -o $(TEST_BIN)

clean:
	rm -f *.o *~ core *.x

