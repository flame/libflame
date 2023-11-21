###############################################################################
# Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################
# 
# Makefile for libflame test suite driver.
#

CC             := gcc
LINKER         := $(CC)
RM_RF          := rm -rf
OMP_LIB        := -fopenmp
CPPFLAGS       :=
CFLAGS         := -std=c11 -D_GNU_SOURCE -O2 -Wall -Wno-comment ${OMP_LIB}
LDFLAGS        := -lstdc++ -lm -lpthread $(OMP_LIB)

#
# --- LIBFLAME \ LAPACK library ------------------------------------------------
#
LIB_PATH          := ../../lib/x86_64-unknown-linux-gnu/
LIBFLAME          := $(LIB_PATH)/libflame.a
LAPACK_INC_PATH   := ../../include/x86_64-unknown-linux-gnu/

#
# --- BLAS library -------------------------------------------------------------
#
# A BLAS library must be given in order to run the libflame test suite.
LIBBLAS           :=
BLAS_HEADER_PATH  :=

#
# ----- AOCL-Utils library path
#
LIBAOCLUTILS_LIBRARY_PATH := 

#
# --- Uncomment the four lines below and set MKLROOT & LDFLAGS -----------------
#     appropriately to link to MKL libs
#
#LIBFLAME       :=
#LIBBLAS        :=
#MKLROOT        :=
#LDFLAGS        := -Wl,--start-group ${MKLROOT}/intel64/libmkl_intel_lp64.a \
                   ${MKLROOT}/intel64/libmkl_sequential.a \
                   ${MKLROOT}/intel64/libmkl_core.a -Wl,--end-group \
                   -lpthread -lm -fopenmp -ldl -lstdc++

#
# --- ILP64 --- Set to 1 to enable ILP64 ---------------------------------------
#
ifeq ($(ILP64),1)
	CPPFLAGS   += -DFLA_ENABLE_ILP64=1
endif

#
# --- Flag to enable uppercase/lowercase API calling conventions --------------
#     Set API_CALLING_CONVENTION = upper  for uppercase without underscore
#                                = upper_ for uppercase with underscore
#                                = lower  for lowercase without underscore
#                                = default for default case (lower_)
#     NOTE: LAPACK default is lowercase with underscore
#
API_CALLING_CONVENTION := default

#
# -- API_CALLING_CONVENTION set to "upper" enables uppercase w/o underscore ----
#
ifeq ($(API_CALLING_CONVENTION),upper)
    CPPFLAGS   += -DUPPER=1
endif

#
# -- API_CALLING_CONVENTION set to "upper_" enables uppercase w/ underscore ----
#
ifeq ($(API_CALLING_CONVENTION),upper_)
    CPPFLAGS   += -DUPPER_=1
endif

#
# -- API_CALLING_CONVENTION set to "lower" enables lowercase w/o underscore ----
#
ifeq ($(API_CALLING_CONVENTION),lower)
    CPPFLAGS   += -DLOWER=1
endif

# --- enable dynamic memory  unaligned  addres logic --------------------------
ifeq ($(MEM_UNALN),1)
	CPPFLAGS   += -DFLA_MEM_UNALIGNED=1
endif
#
# --- General build definitions ------------------------------------------------
#
TEST_SRC_PATH  := src
TEST_OBJ_PATH  := src/obj
TEST_VALIDATE_SRC_PATH := validate_src
TEST_VALIDATE_OBJ_PATH := validate_src/obj

CFLAGS         += -I$(BLAS_HEADER_PATH) -I$(LAPACK_INC_PATH) -I$(TEST_SRC_PATH) -I$(TEST_VALIDATE_SRC_PATH)

FNAME          := lapack

TEST_VALIDATE_OBJS := $(patsubst $(TEST_VALIDATE_SRC_PATH)/%.c, \
                             $(TEST_VALIDATE_OBJ_PATH)/%.o, \
                             $(wildcard $(TEST_VALIDATE_SRC_PATH)/*.c))

TEST_OBJS      := $(patsubst $(TEST_SRC_PATH)/%.c, \
                             $(TEST_OBJ_PATH)/%.o, \
                             $(wildcard $(TEST_SRC_PATH)/*.c))

TEST_BIN       := test_$(FNAME).x

$(TEST_VALIDATE_OBJ_PATH)/%.o: $(TEST_VALIDATE_SRC_PATH)/%.c
	mkdir -p $(TEST_VALIDATE_OBJ_PATH)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TEST_OBJ_PATH)/%.o: $(TEST_SRC_PATH)/%.c
	mkdir -p $(TEST_OBJ_PATH)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_VALIDATE_OBJS) $(TEST_OBJS)
	$(LINKER) $(TEST_VALIDATE_OBJS) $(TEST_OBJS) $(LIBFLAME) $(LIBBLAS) $(LIBAOCLUTILS_LIBRARY_PATH) \
	$(LDFLAGS) -o $(TEST_BIN)

clean:
	$(RM_RF) $(TEST_VALIDATE_OBJ_PATH) $(TEST_OBJ_PATH) $(TEST_BIN)
