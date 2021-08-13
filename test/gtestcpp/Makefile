###############################################################################
# Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################

#
# --- Makefile PHONY target definitions ----------------------------------------
#

.PHONY: all \
        libflame \
        clean cleanx

#
# --- Include libflame config makefile fragment --------------------------------
#

# Determine the path to the libflame config makefile fragment. We'll use
# several variables defined there.
BUILD_DIR       := ../../build
CONFIG_DIR      := ../../config
LIB_DIR         := ../../lib
INCLUDE_DIR     := ../../src/aocl_dtl -I../../windows/include -I../../config/x86_64-unknown-linux-gnu
SRC_PATH        := ../../src
HOST            := $(shell sh $(BUILD_DIR)/ac-utils/config.guess)
CONFIG_MK_FILE  := $(CONFIG_DIR)/$(HOST)/config.mk
LIB_PATH        := $(LIB_DIR)/$(HOST)
INC_PATH        := $(INCLUDE_DIR)


# Include the definitions in the config makefile fragment.
-include $(CONFIG_MK_FILE)

#
# --- Optional overrides -------------------------------------------------------
#

# Uncomment and modify these definitions if you wish to override the values
# present in the master config makefile fragment.
CXX            := g++
LINKER         := $(CXX)

CXXFLAGS       += -g -Wno-unused-result -Wno-comment -fopenmp
LDFLAGS        := -lm -fopenmp -lgfortran 

#
# --- BLAS and LAPACK implementations ------------------------------------------
#

# BLAS implementation path. A BLAS library must be given in order to run
# the libflame test suite. Modify these definitions if needed.
LIBBLAS_PATH   := ../../blis/lib/
LIBBLAS        := $(LIBBLAS_PATH)/libblis-mt.a

# LAPACK implementation path. These values only matter if libflame was
# configured with the external-lapack-interfaces option enabled. Modify
# these definitions if needed.
LIBLAPACK      :=

all: libflame
libflame: libflame_hbev.x \
libflame_hbev_2stage.x \
libflame_hbevd.x \
libflame_hbevd_2stage.x \
libflame_hbevx.x \
libflame_hbevx_2stage.x \
libflame_hbgst.x \
libflame_hbgv.x \
libflame_hbgvd.x \
libflame_hbgvx.x \
libflame_hbtrd.x \
libflame_heequb.x \
libflame_heev.x \
libflame_heev_2stage.x \
libflame_heevd.x \
libflame_heevd_2stage.x \
libflame_heevr.x \
libflame_heevr_2stage.x \
libflame_heevx.x \
libflame_heevx_2stage.x

# --- General build definitions ------------------------------------------------
#

TEST_SRC_PATH := src
TEST_OBJ_PATH := obj
TEST_BIN_PATH := bin

#  GTEST Related build variables
GTESTDIR       := ../../googletest/
GTESTLIBDIR    := $(GTESTDIR)/build/lib
GTESTINC       := $(GTESTDIR)/include  
GTESTLIB       := $(GTESTLIBDIR)/libgtest.a

FLA_LIB_PATH   := $(LIB_PATH)
FLA_INC_PATH   := $(INC_PATH)
LIBFLAMELIB    := ../../lib/x86_64-unknown-linux-gnu/libflame.a

LIBRARIES      = $(LIBFLAMELIB) $(LIBBLAS) $(GTESTLIB) 
CXXFLAGS      += -I$(FLA_INC_PATH) -I$(TEST_SRC_PATH) -I$(SRC_PATH)/src_cpp  -I$(GTESTINC)

TEST_OBJS      := $(patsubst $(TEST_SRC_PATH)/%.c, \
                             $(TEST_OBJ_PATH)/%.o, \
                             $(wildcard $(TEST_SRC_PATH)/*.c))

$(TEST_OBJ_PATH)/%.o: $(TEST_SRC_PATH)/%.cc
	@mkdir -p $(TEST_OBJ_PATH)
	@$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_OBJ_PATH)/libflame_%.o: libflame_%.cc
	@$(CXX) $(CXXFLAGS) -c $< -o $@

libflame_%.x: $(TEST_OBJ_PATH)/libflame_%.o $(TEST_OBJ_PATH)/main.o
	$(LINKER) $^ $(TEST_OBJS) -ldl $(LIBRARIES) $(LDFLAGS) -o $@
	./$@

# -- Clean rules --

clean: cleanx

cleanx:
	- $(RM_RF) $(TEST_OBJ_PATH)
	- $(RM_F) *.x
