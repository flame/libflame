
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
LIBLAPACK      :=  
# -L/opt/local/lib -llapack

#
# --- General build definitions ------------------------------------------------
#

FLA_LIB_PATH   := ../../lib/$(HOST)
FLA_INC_PATH   := ../../include_local
LIBFLAME       := $(FLA_LIB_PATH)/libflame.a

CFLAGS         += -I$(FLA_INC_PATH) -I.

TEST  = tridiag_realify_diagonals
TESTS = qr rq lq ql qr_piv \
	pivot normalize bidiag_realify_diagonals tridiag_realify_diagonals \
	extract_component givens constant \
	bidiag bidiag_ext bidiag_external tridiag tridiag_external ls_biqr \
	svd bsvd  bsvd_ext svd_ext

%.o : %.c
	@echo "Compiling $<"
	$(CC) $(CFLAGS) -c $< -o $@

one : $(TEST).o
	@echo "Linking   $<"
	$(LINKER) -o $(TEST).x $(TEST).o  $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS)

all :
	@for e in $(TESTS) ; do \
		make one TEST=$$e ;\
	done

clean :
	@rm -f *.o *~
	@for t in $(TESTS) ; do rm -f $$t.x ; done
	@for d in $(DIRS); do (cd $$d; make clean ); done

