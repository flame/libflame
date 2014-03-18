#
# test directory makefile
#

FNAME        := plasma

SRC_PATH     := ..
OBJ_PATH     := .

HOME         := /home/rvdg

LIB_PATH     := $(HOME)/plasma/2.0.0/lib
INC_PATH     := $(HOME)/plasma/2.0.0/include

PLASMA       := $(LIB_PATH)/libplasma.a $(LIB_PATH)/libcoreblas.a $(LIB_PATH)/libcorelapack.a
FLAME        := $(LIB_PATH)/libflame.a
GOTO         := $(LIB_PATH)/libgoto.a
MKL          := -L/opt/intel/mkl/8.1/lib/64/ -lmkl -lguide

CC           := icc
FC           := ifort
LINKER       := $(CC)
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O3
LDFLAGS      := -lm -lpthread -openmp -L/opt/intel/fc/9.0/lib -L/usr/lib -lifport -lifcore -limf -lirc -lunwind -lirc_s -ldl

TEST_BIN     := $(FNAME).x
TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

compile: $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(PLASMA) $(FLAME) $(MKL) $(GOTO) -o $(TEST_BIN)

run:
	dplace -x2 -e -c0-15 ./$(TEST_BIN) < input

clean:
	rm -f *.o *~ core *.x

remove:
	rm ./results/*.m
