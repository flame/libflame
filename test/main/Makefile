###############################################################################
# Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################
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

#
# --- ILP64 --------------------------------------------------------------------
#
ILP64      := 
ILP64FLAG := FLA_ENABLE_ILP64=$(ILP64)

#
# --- General build definitions ------------------------------------------------
#
TEST_SRC_PATH  := src
TEST_OBJ_PATH  := obj
TEST_LIBRARY_SRC_PATH := ../library
TEST_LIBRARY_OBJ_PATH := ../library/obj

CFLAGS         += -I$(TEST_SRC_PATH) -I$(TEST_LIBRARY_SRC_PATH)

ifeq ($(ILP64),1)
	CPPFLAGS	+= -D$(ILP64FLAG)
endif

FNAME          := lapack

TEST_LIBRARY_OBJS := $(patsubst $(TEST_LIBRARY_SRC_PATH)/%.c, \
                             $(TEST_LIBRARY_OBJ_PATH)/%.o, \
                             $(wildcard $(TEST_LIBRARY_SRC_PATH)/*.c))

TEST_OBJS      := $(patsubst $(TEST_SRC_PATH)/%.c, \
                             $(TEST_OBJ_PATH)/%.o, \
                             $(wildcard $(TEST_SRC_PATH)/*.c))

TEST_BIN       := test_$(FNAME).x

$(TEST_LIBRARY_OBJ_PATH)/%.o: $(TEST_LIBRARY_SRC_PATH)/%.c
	mkdir -p $(TEST_LIBRARY_OBJ_PATH)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TEST_OBJ_PATH)/%.o: $(TEST_SRC_PATH)/%.c
	mkdir -p $(TEST_OBJ_PATH)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_LIBRARY_OBJS) $(TEST_OBJS)
	$(LINKER) $(TEST_LIBRARY_OBJS) $(TEST_OBJS) $(LIBFLAME) $(LIBBLAS) $(LDFLAGS) -o $(TEST_BIN)

clean:
	$(RM_F) $(TEST_LIBRARY_OBJS) $(TEST_OBJS) $(TEST_BIN)
