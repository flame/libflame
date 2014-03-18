#
# test directory makefile
#

FNAME        := supermatrix

SRC_PATH     := ..
OBJ_PATH     := .

LIB_PATH     := $(HOME)/flame/lib
INC_PATH     := $(HOME)/flame/include

FLAME        := $(LIB_PATH)/libflame.a
BLAS         := $(LIB_PATH)/libgoto.a

CC           := gcc
LINKER       := $(CC)
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O3
LDFLAGS      := -lm -lpthread

TEST_BIN     := $(FNAME).x
TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

compile: $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(FLAME) $(BLAS) -o $(TEST_BIN)

run:
	./$(TEST_BIN) < input

clean:
	rm -f *.o *~ core *.x

remove:
	rm ./results/*.m
