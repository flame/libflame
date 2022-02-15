# 
# Makefile for libflame test suite driver.
#


CC             := gcc
LINKER         := $(CC)
RM_F           := rm -f
CFLAGS         := -g -O2 -Wall -Wno-comment
LDFLAGS        := -lm -lpthread -fopenmp


#
# --- BLAS library -------------------------------------------------------------
#

# A BLAS library must be given in order to run the libflame test suite.
LIBBLAS        := 


#
# --- LIBFLAME library ---------------------------------------------------------
#

LIB_PATH       := ../../lib/x86_64-unknown-linux-gnu/
LIBFLAME       := $(LIB_PATH)/libflame.a
#LIBFLAME       := $(LIB_PATH)/libflame.so


#
# --- ILP64 --------------------------------------------------------------------
#

ILP64      := 
ILP64FLAG := ENABLE_ILP64=$(ILP64)


#
# --- General build definitions ------------------------------------------------
#

TEST_SRC_PATH  := src
TEST_OBJ_PATH  := obj


CFLAGS         += -I$(TEST_SRC_PATH)

ifeq ($(ILP64),1)
	CPPFLAGS	+= -D$(ILP64FLAG)
endif

FNAME          := libflame

TEST_OBJS      := $(patsubst $(TEST_SRC_PATH)/%.c, \
                             $(TEST_OBJ_PATH)/%.o, \
                             $(wildcard $(TEST_SRC_PATH)/*.c))
TEST_BIN       := test_$(FNAME).x

$(TEST_OBJ_PATH)/%.o: $(TEST_SRC_PATH)/%.c
	mkdir -p $(TEST_OBJ_PATH)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LIBFLAME) $(LIBBLAS) $(LDFLAGS) -o $(TEST_BIN)

clean:
	$(RM_F) $(TEST_OBJS) $(TEST_BIN)