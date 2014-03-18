#
# test directory makefile
#

FNAME        := Symm_ru

SRC_PATH     := ..
OBJ_PATH     := .

LIB_PATH     := $(HOME)/flame/lib
INC_PATH     := $(HOME)/flame/include
FLAME_BASE   := $(LIB_PATH)/libflame-base_ia64.a
FLAME_LAPACK := $(LIB_PATH)/libflame-lapack_ia64.a
BLAS_LIB     := $(LIB_PATH)/libgoto_ia64-it2.a
#BLAS_LIB     := $(LIB_PATH)/libgoto_ia64-it2p.a
LAPACK_LIB   := $(LIB_PATH)/liblapack_ia64-v3.0.a

CC           := icc
FC           := ifort
LINKER       := $(CC)
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O3
FFLAGS       := $(CFLAGS)
LDFLAGS      := -lm -lpthread

TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))
TEST_OBJS    += $(OBJ_PATH)/test_$(FNAME).o \
              $(OBJ_PATH)/time_$(FNAME).o \
              $(OBJ_PATH)/REF_$(FNAME).o
TEST_BIN     := test_$(FNAME).x

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJ_PATH)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME_BLAS) $(FLAME_BASE) $(BLAS_LIB) -o $(TEST_BIN)

clean:
	rm -f *.o *~ core *.x

