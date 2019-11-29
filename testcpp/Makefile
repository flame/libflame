

#
# Makefile
#
# 
#


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
BUILD_DIR       := ../build
CONFIG_DIR      := ../config
LIB_DIR         := ../lib
INCLUDE_DIR     := ../src/base/flamec/include
SRC_PATH        := ../src
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
CC             := g++
LINKER         := $(CC)


CXXFLAGS       += -O2 -Wall -Wno-comment -fopenmp
LDFLAGS        := -lm -fopenmp -lgfortran 



#
# --- BLAS and LAPACK implementations ------------------------------------------
#

# BLAS implementation path. A BLAS library must be given in order to run
# the libflame test suite. Modify these definitions if needed.
LIBBLAS_PATH   := 
LIBBLAS        := $(LIBBLAS_PATH)/libblis.a

# LAPACK implementation path. These values only matter if libflame was
# configured with the external-lapack-interfaces option enabled. Modify
# these definitions if needed.
LIBLAPACK      :=

all: libflame
libflame: libflame_potrf.x \
	libflame_potf2.x \
	libflame_getrf.x \
	libflame_getf2.x \
	libflame_geqrf.x \
	libflame_geqr2.x \
	libflame_geqpf.x \
	libflame_geqp3.x \
	libflame_gelqf.x \
	libflame_gelq2.x \
	libflame_gelsd.x \
	libflame_gelss.x \
	libflame_lauum.x \
	libflame_lauu2.x \
	libflame_potri.x \
	libflame_trtri.x \
	libflame_trti2.x \
	libflame_trsyl.x \
	libflame_gehrd.x \
	libflame_gehd2.x \
	libflame_sytrd.x \
	libflame_hetrd.x \
	libflame_sytd2.x \
	libflame_hetd2.x \
	libflame_gebrd.x \
	libflame_gebd2.x \
	libflame_sygst.x \
	libflame_hegst.x \
	libflame_sygs2.x \
	libflame_hegs2.x \
	libflame_larft.x \
	libflame_larfg.x \
	libflame_larfgp.x \
	libflame_orgqr.x \
	libflame_ungqr.x \
	libflame_ormqr.x \
	libflame_unmqr.x \
	libflame_orm2r.x \
	libflame_unm2r.x \
	libflame_orglq.x \
	libflame_unglq.x \
	libflame_ormlq.x \
	libflame_unmlq.x \
	libflame_orml2.x \
	libflame_unml2.x \
	libflame_orgtr.x \
	libflame_ungtr.x \
	libflame_ormtr.x \
	libflame_unmtr.x \
	libflame_orgbr.x \
	libflame_ungbr.x \
	libflame_ormbr.x \
	libflame_unmbr.x \
	libflame_steqr.x \
	libflame_stedc.x \
	libflame_stemr.x \
	libflame_syev.x \
	libflame_heev.x \
	libflame_syevd.x \
	libflame_heevd.x \
	libflame_syevr.x \
	libflame_heevr.x \
	libflame_bdsqr.x \
	libflame_bdsdc.x \
	libflame_gesvd.x \
	libflame_gesdd.x \
	libflame_laswp.x \
	libflame_laset.x

#
# --- General build definitions ------------------------------------------------
#

TEST_SRC_PATH  := src
TEST_OBJ_PATH  := obj
TEST_BIN_PATH := bin

FLA_LIB_PATH   := $(LIB_PATH)
FLA_INC_PATH   := $(INC_PATH)
LIBFLAME       := $(FLA_LIB_PATH)/libflame.a
LAPACKE        := $(SRC_PATH)/lapacke/liblapacke.a

CXXFLAGS         += -I$(FLA_INC_PATH) -I$(SRC_PATH)/testcpp/src -I$(SRC_PATH)/src_cpp -I$(SRC_PATH)/lapacke/LAPACKE/include -I$(TEST_SRC_PATH) -I$(SRC_PATH)/base/flamec/blis/include -I$(SRC_PATH)/base/flamec/control  -I$(SRC_PATH)/base/flamec/include
FNAME          := libflame

TEST_OBJS      := $(patsubst $(TEST_SRC_PATH)/%.c, \
                             $(TEST_OBJ_PATH)/%.o, \
                             $(wildcard $(TEST_SRC_PATH)/*.c))

$(TEST_OBJ_PATH)/%.o: $(TEST_SRC_PATH)/%.cc
	@mkdir -p $(TEST_OBJ_PATH)
	@$(CC) $(CXXFLAGS) -c $< -o $@

$(TEST_OBJ_PATH)/libflame_%.o: libflame_%.cc
	@$(CC) $(CXXFLAGS) -c $< -o $@

libflame_%.x: $(TEST_OBJ_PATH)/libflame_%.o 
	@$(LINKER) $^ $(TEST_OBJS) $(LAPACKE) $(LIBFLAME) $(LIBBLAS) $(LDFLAGS) -o $@
	./$@

# -- Clean rules --

clean: cleanx

cleanx:
	- $(RM_RF) $(TEST_OBJ_PATH)
	- $(RM_F) *.x
