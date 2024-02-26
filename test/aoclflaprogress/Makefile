

#
# Makefile
#
# Makefile for libflame_aocl_progress driver.
#



#
# --- Include libflame config makefile fragment --------------------------------
#

# Determine the path to the libflame config makefile fragment. We'll use
# several variables defined there.
BUILD_DIR       := ../../build
CONFIG_DIR      := ../../config
LIB_DIR         := ../../lib
INCLUDE_DIR     := ../../include
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
# CC             := gcc
# LINKER         := $(CC)
# CFLAGS         := -g -O2 -Wall -Wno-comment
LDFLAGS        := -lm -lpthread
# INSTALL_PREFIX := $(HOME)/flame



#
# --- BLAS and LAPACK implementations ------------------------------------------
#

# BLAS implementation path. A BLAS library must be given in order to run
# the libflame test suite. Modify these definitions if needed.
LIBBLAS_PATH   := $(INSTALL_LIBDIR)
#LIBBLAS        := $(LIBBLAS_PATH)/libblas.a
#LIBBLAS        := $(LIBBLAS_PATH)/libopenblas.a
#LIBBLAS        := $(HOME)/blis/lib/libblis.a
LIBBLAS        := 

# LAPACK implementation path. These values only matter if libflame was
# configured with the external-lapack-interfaces option enabled. Modify
# these definitions if needed.
LIBLAPACK_PATH := $(INSTALL_LIBDIR)
LIBLAPACK      := 

#
# ----- AOCL-Utils library path
#
LIBAOCLUTILS_LIBRARY_PATH := 


#
# --- General build definitions ------------------------------------------------
#

TEST_SRC_PATH  := src
TEST_OBJ_PATH  := obj

#FLA_LIB_PATH   := $(INSTALL_PREFIX)/lib
#FLA_INC_PATH   := $(INSTALL_PREFIX)/include
FLA_LIB_PATH   := $(LIB_PATH)
FLA_INC_PATH   := $(INC_PATH)
LIBFLAME       := $(FLA_LIB_PATH)/libflame.a
#LIBFLAME       := $(FLA_LIB_PATH)/libflame.so

CFLAGS         += -I$(FLA_INC_PATH) -I$(TEST_SRC_PATH)

FNAME          := libFLAME_aocl

TEST_OBJS      := $(patsubst $(TEST_SRC_PATH)/%.c, \
                             $(TEST_OBJ_PATH)/%.o, \
                             $(wildcard $(TEST_SRC_PATH)/*.c))
TEST_BIN       := test_$(FNAME).x

$(TEST_OBJ_PATH)/%.o: $(TEST_SRC_PATH)/%.c
	mkdir -p $(TEST_OBJ_PATH)
	$(CC) $(CFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) $(LIBAOCLUTILS_LIBRARY_PATH) $(LDFLAGS) -o $(TEST_BIN)

clean:
	$(RM_F) $(TEST_OBJS) $(TEST_BIN)

