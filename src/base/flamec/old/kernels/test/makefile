#
# test directory makefile
#

FNAME        := Gemm_pp_nn

SRC_PATH     := ..
OBJ_PATH     := .
BLAS_PATH    := /home/rvdg/FLAME_GOTO/
LIB_PATH     := $(HOME)/flame/lib
INC_PATH     := $(HOME)/flame/include 
BLAS_LIB     := $(BLAS_PATH)/libgoto_pentium4_512-r0.96.a
FLAME_BASE   := $(LIB_PATH)/libflame-base.a
FLAME_BLAS   := $(LIB_PATH)/libflame-blas.a

CC           := gcc
FC           := g77
LINKER       := $(CC)
CFLAGS       :=  -I$(SRC_PATH) -I$(INC_PATH) \
              -g -O3 -fomit-frame-pointer -funroll-all-loops -finline-functions \
              -Wmissing-prototypes
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

