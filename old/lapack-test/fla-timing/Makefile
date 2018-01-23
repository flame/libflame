
#
# Makefile
#
# Field G. Van Zee
# 
# Makefile for libflame test suite driver.
#

#
# --- Include libflame config makefile fragment --------------------------------
#

# Determine the path to the libflame config makefile fragment. We'll use
# several variables defined there.
BUILD_DIRPATH       := ../../build
CONFIG_DIRPATH      := ../../config
HOST                := $(shell sh $(BUILD_DIRPATH)/ac-utils/config.guess)
CONFIG_MK_FRAGMENT  := $(CONFIG_DIRPATH)/$(HOST)/config.mk

# Include the definitions in the config makefile fragment.
-include $(CONFIG_MK_FRAGMENT)

#
# --- Optional overrides -------------------------------------------------------
#

# Uncomment and modify these definitions if you wish to override the values
# present in the master config makefile fragment.
# CC             := gcc
# LINKER         := $(CC)
# CFLAGS         := -g -O2 -Wall -Wno-comment
# LDFLAGS        := 
INSTALL_PREFIX := /Users/kyungjoo/Lib/FLAME



#
# --- BLAS and LAPACK implementations ------------------------------------------
#

# BLAS implementation path. A BLAS library must be given in order to run
# the libflame test suite. Modify these definitions if needed.
LIBBLAS        := -L/opt/local/lib -lcblas -lf77blas -latlas -lf2c 

# LAPACK implementation path. These values only matter if libflame was
# configured with the external-lapack-interfaces option enabled. Modify
# these definitions if needed.
LIBLAPACK      := -L/opt/local/lib -llapack

#
# --- General build definitions ------------------------------------------------
#

TEST_SRC_PATH  := src
TEST_OBJ_PATH  := obj

FLA_LIB_PATH   := ../../lib/$(HOST)
FLA_INC_PATH   := ../../include
LIBFLAME       := $(FLA_LIB_PATH)/libflame.a

CFLAGS         += -I$(FLA_INC_PATH) -I$(TEST_SRC_PATH)

FNAME          := lapack2flame

TEST_OBJS      := $(patsubst $(TEST_SRC_PATH)/%.c, \
                             $(TEST_OBJ_PATH)/%.o, \
                             $(wildcard $(TEST_SRC_PATH)/*.c))

TEST_BIN       := test_$(FNAME).x

$(TEST_OBJ_PATH)/%.o: $(TEST_SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o $(TEST_BIN)

clean:
	$(RM_F) $(TEST_OBJS) $(TEST_BIN) *~

