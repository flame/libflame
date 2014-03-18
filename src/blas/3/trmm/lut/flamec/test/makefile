#
# test directory makefile
#

FNAME        := Trmm_lut

SRC_PATH     := ..
OBJ_PATH     := .

LIB_PATH     := $(HOME)/flame/lib
INC_PATH     := $(HOME)/flame/include
FLAME_BASE   := $(LIB_PATH)/libflame.a
BLAS_LIB     := $(LIB_PATH)/libgoto_ia64.a

CC           := icc
FC           := ifort
LINKER       := $(CC)
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O3
FFLAGS       := $(CFLAGS)
LDFLAGS      := -L/opt/intel/fc/ia64/10.0.026/lib -L/usr/lib/gcc/ia64-unknown-linux-gnu/4.1.2/ -L/usr/lib/gcc/ia64-unknown-linux-gnu/4.1.2/../../../ -lifport -lifcore -limf -lm -lipgo -lirc -lirc_s -ldl -lm -lpthread

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

