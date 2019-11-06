

#
# Makefile
#
# 
#



#
# --- Include libflame config makefile fragment --------------------------------
#

# Determine the path to the libflame config makefile fragment. We'll use
# several variables defined there.
BUILD_DIR       := ../build
CONFIG_DIR      := ../config
LIB_DIR         := ../lib
INCLUDE_DIR     := ../include
SRC_PATH         := ../src
HOST            := $(shell sh $(BUILD_DIR)/ac-utils/config.guess)
CONFIG_MK_FILE  := $(CONFIG_DIR)/$(HOST)/config.mk
LIB_PATH        := $(LIB_DIR)/$(HOST)/
INC_PATH        := $(INCLUDE_DIR)/$(HOST)/

# Include the definitions in the config makefile fragment.
-include $(CONFIG_MK_FILE)



#
# --- Optional overrides -------------------------------------------------------
#

# Uncomment and modify these definitions if you wish to override the values
# present in the master config makefile fragment.
CC             := g++
LINKER         := $(CC)
CFLAGS         += -g -O2 -Wall -Wno-comment -fopenmp
LDFLAGS        := -lm -fopenmp -lgfortran 
# INSTALL_PREFIX := $(HOME)/flame



#
# --- BLAS and LAPACK implementations ------------------------------------------
#

# BLAS implementation path. A BLAS library must be given in order to run
# the libflame test suite. Modify these definitions if needed.
#LIBBLAS_PATH   := /home/pradeep/amd/aocl/2.0/libs
LIBBLAS_PATH   := /home/managala/2019/blis/lib/zen
LIBBLAS        := $(LIBBLAS_PATH)/libblis.a

# LAPACK implementation path. These values only matter if libflame was
# configured with the external-lapack-interfaces option enabled. Modify
# these definitions if needed.
LIBLAPACK      :=



#
# --- General build definitions ------------------------------------------------
#

TEST_SRC_PATH  := src
TEST_SRC := libflame_potrf
TEST_OBJ_PATH  := obj
TEST_BIN_PATH := bin

FLA_LIB_PATH   := $(LIB_PATH)
FLA_INC_PATH   := $(INC_PATH)
LIBFLAME       := $(FLA_LIB_PATH)/libflame.a
#LIBFLAME       := $(FLA_LIB_PATH)/libflame.so
LAPACKE        := $(SRC_PATH)/lapacke/liblapacke.a

CFLAGS         += -I$(FLA_INC_PATH) -I$(SRC_PATH)/src_cpp -I$(SRC_PATH)/lapacke/LAPACKE/include -I$(TEST_SRC_PATH) 
FNAME          := libflame

TEST_OBJS      := $(patsubst $(TEST_SRC_PATH)/$(TEST_SRC).cc, \
                             $(TEST_OBJ_PATH)/$(TEST_SRC).o, \
                             $(wildcard $(TEST_SRC_PATH)/$(TEST_SRC).cc))
TEST_BIN       := $(TEST_SRC)

$(TEST_OBJ_PATH)/%.o: $(TEST_SRC_PATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

test_$(TEST_SRC): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LAPACKE) $(LIBFLAME) $(LIBBLAS) $(LDFLAGS) -o $(TEST_BIN)

clean:
	$(RM_F) $(TEST_OBJS) $(TEST_BIN)

