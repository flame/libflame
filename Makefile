
#
# Makefile
#
# Field G. Van Zee
#
# Top-level makefile for libflame linear algebra library.
#
#

#
# --- Makefile PHONY target definitions ----------------------------------------
#

.PHONY: all \
        libs libflame \
        check-env check-env-config check-env-fragments \
        flat-headers \
        test check checklibflame\
	checkcpp cleancpptest \
        install-headers install-libs install-lib-symlinks \
        clean cleanmk cleanh cleanlib cleanleaves distclean \
        install \
        uninstall-libs uninstall-lib-symlinks uninstall-headers


# Accept an abbreivated request for verbosity (e.g. 'make V=1 ...')
ifeq ($(V),1)
ENABLE_VERBOSE := yes
endif

ifeq ($(V),0)
ENABLE_VERBOSE := no
endif



#
# --- Include architecture-specific variable definitions ----------------------
#

# Avoid catting the config.sys_type file unless it exists. This makes the
# output of things like 'make distclean' (when the directory is already
# clean) less confusing.
ifneq ($(wildcard config.sys_type),)
# The host string will uniquely identify the current host (as much
# as is reasonable) for purposes of separating the configure products and
# object files of one architecture from another.
HOST            := $(shell cat config.sys_type)
else
HOST            := unknown-generic
endif

# --- Important directories for source code and build products ---

SRC_DIR         := src
CONFIG_DIR      := config
OBJ_DIR         := obj
LIB_DIR         := lib
INC_DIR         := include
LAPACKE_DIR     := lapacke
AOCLDTL_DIR     := aocl_dtl
LAPACKE_SRC_DIR := LAPACKE/src
LAPACKE_UTIL_DIR:= LAPACKE/utils
TEST_DIR        := test
CPP_TEST_DIR    := testcpp

# Use the system type to name the config, object, and library directories.
# These directories are special in that they will contain products specific
# to this particular architecture.
BASE_CONFIG_PATH := ./$(CONFIG_DIR)/$(HOST)
BASE_OBJ_PATH    := ./$(OBJ_DIR)/$(HOST)/$(SRC_DIR)
BASE_LIB_PATH    := ./$(LIB_DIR)/$(HOST)
BASE_INC_PATH    := ./$(INC_DIR)/$(HOST)

# Pathname to the makefile fragment containing lots of various definitions,
# most of which were substituted via configure.
CONFIG_MK_FILE   := $(BASE_CONFIG_PATH)/config.mk

# Include the definitions in the config makefile fragment.
-include $(CONFIG_MK_FILE)

# Detect whether we actually got the config makefile fragment. If we didn't,
# then it is likely that the user has not yet generated it (via configure).
ifeq ($(strip $(CONFIG_MK_INCLUDED)),yes)
CONFIG_MK_PRESENT := yes
IS_CONFIGURED     := yes
else
CONFIG_MK_PRESENT := no
IS_CONFIGURED     := no
endif



#
# --- Library name and local paths ---------------------------------------------
#

# --- Shared library extension ---

# The shared (dynamic) library file suffix is different for Linux and OS X.
ifeq ($(OS_NAME),Darwin)
SHLIB_EXT          := dylib
else
SHLIB_EXT          := so
endif

# --- Library names ---

LIBFLAME             := libflame

LIBFLAME_A           := $(LIBFLAME).a
LIBFLAME_SO          := $(LIBFLAME).$(SHLIB_EXT)

LAPACKE_A	     := liblapacke.a
AOCLDTL_A            := libaocldtl.a
AOCLDTL_SO           := libaocldtl.so

# --- Library filepaths ---

# Append the base library path to the library names.
LIBFLAME_A_PATH      := $(BASE_LIB_PATH)/$(LIBFLAME_A)
LIBFLAME_SO_PATH     := $(BASE_LIB_PATH)/$(LIBFLAME_SO)

ifeq ($(OS_NAME),Darwin)
# OS X shared library extensions.
LIBFLAME_SO_MAJ_EXT  := $(SO_MAJOR).$(SHLIB_EXT)
LIBFLAME_SO_MMB_EXT  := $(SO_MMB).$(SHLIB_EXT)
else
# Linux shared library extensions.
LIBFLAME_SO_MAJ_EXT  := $(SHLIB_EXT).$(SO_MAJOR)
LIBFLAME_SO_MMB_EXT  := $(SHLIB_EXT).$(SO_MMB)
endif

LIBFLAME_SONAME      := $(LIBFLAME).$(LIBFLAME_SO_MAJ_EXT)
LIBFLAME_SO_MAJ_PATH := $(BASE_LIB_PATH)/$(LIBFLAME_SONAME)

LAPACKE_A_PATH       := $(SRC_DIR)/$(LAPACKE_DIR)/$(LAPACKE_A)
AOCLDTL_A_PATH       := $(SRC_DIR)/$(AOCLDTL_DIR)/$(AOCLDTL_A)
AOCLDTL_SO_PATH      := $(SRC_DIR)/$(AOCLDTL_DIR)/$(AOCLDTL_SO)
AOCLDTL_obj_PATH     := $(SRC_DIR)/$(AOCLDTL_DIR)/*.o
AOCLDTL_gch_PATH     := $(SRC_DIR)/$(AOCLDTL_DIR)/*.gch
LAPACKE_S_OBJS_PATH  := $(SRC_DIR)/$(LAPACKE_DIR)/$(LAPACKE_SRC_DIR)/*.o
LAPACKE_U_OBJS_PATH  := $(SRC_DIR)/$(LAPACKE_DIR)/$(LAPACKE_UTIL_DIR)/*.o
# Construct the output path when building a shared library.
LIBFLAME_SO_OUTPUT_NAME := $(LIBFLAME_SO_PATH)



#
# --- Shared library flags ----------------------------------------------------
#

# Specify the shared library's 'soname' field.
# NOTE: The flag for creating shared objects is different for Linux and OS X.
ifeq ($(OS_NAME),Darwin)
# OS X shared library link flags.
SOFLAGS    := -dynamiclib
SOFLAGS    += -Wl,-install_name,$(LIBFLAME_SONAME)
else
SOFLAGS    := -shared
# Linux shared library link flags.
SOFLAGS    += -Wl,-soname,$(LIBFLAME_SONAME)
endif



#
# --- Main target variable definitions ----------------------------------------
#

MK_BASE_FLAMEC_SRC                    :=
MK_BASE_FLAMEC_OBJS                   :=

MK_BLAS_FLAMEC_SRC                    :=
MK_BLAS_FLAMEC_OBJS                   :=

MK_LAPACK_FLAMEC_SRC                  :=
MK_LAPACK_FLAMEC_OBJS                 :=


MK_MAP_LAPACK2FLAMEC_SRC              :=
MK_MAP_LAPACK2FLAMEC_OBJS             :=

MK_MAP_LAPACK2FLAMEC_F2C_SRC          :=
MK_MAP_LAPACK2FLAMEC_F2C_OBJS         :=

MK_MAP_LAPACK2FLAMEC_F2C_FLAMEC_SRC   :=
MK_MAP_LAPACK2FLAMEC_F2C_FLAMEC_OBJS  :=

MK_MAP_LAPACK2FLAMEC_F2C_INSTALL_SRC  :=
MK_MAP_LAPACK2FLAMEC_F2C_INSTALL_OBJS :=

MK_FLABLAS_F2C_SRC                    :=
MK_FLABLAS_F2C_OBJS                   :=

# --- Define install target names for static libraries ---

LIBFLAME_A_INST       := $(INSTALL_LIBDIR)/$(LIBFLAME_A)

# --- Define install target names for dynamic libraries ---

LIBFLAME_SO_INST      := $(INSTALL_LIBDIR)/$(LIBFLAME_SO)
LIBFLAME_SO_MAJ_INST  := $(INSTALL_LIBDIR)/$(LIBFLAME_SONAME)
LIBFLAME_SO_MMB_INST  := $(INSTALL_LIBDIR)/$(LIBFLAME).$(LIBFLAME_SO_MMB_EXT)
# --- Determine which libraries to build ---

MK_LIBS                   :=
MK_LIBS_INST              :=
MK_LIBS_SYML              :=

ifeq ($(FLA_ENABLE_STATIC_BUILD),yes)
MK_LIBS                   += $(LIBFLAME_A_PATH)
MK_LIBS_INST              += $(LIBFLAME_A_INST)
MK_LIBS_SYML              +=
endif

MK_LIBS                   += $(LAPACKE_A_PATH) \
			     $(AOCLDTL_A_PATH)

ifeq ($(FLA_ENABLE_DYNAMIC_BUILD),yes)
MK_LIBS                   += $(LIBFLAME_SO_PATH) \
                             $(LIBFLAME_SO_MAJ_PATH)
MK_LIBS_INST              += $(LIBFLAME_SO_MMB_INST)
MK_LIBS_SYML              += $(LIBFLAME_SO_INST) \
			     $(LIBFLAME_SO_MAJ_INST)
endif

# Strip leading, internal, and trailing whitespace.
MK_LIBS_INST              := $(strip $(MK_LIBS_INST))
MK_LIBS_SYML              := $(strip $(MK_LIBS_SYML))

# Set the path to the subdirectory of the include installation directory.
MK_INCL_DIR_INST          := $(INSTALL_INCDIR)



#
# --- Test definitions --------------------------------------------------------
#





#
# --- Include sub-tree makefile fragments --------------------------------------
#

# Makefile fragment name.
FRAGMENT_MK     := .fragment.mk

# Initialize our list of directory paths to makefile fragments with the empty list.
FRAGMENT_DIR_PATHS :=

# The only fragment subdirectory that we build from the top-level directory is the
# source directory.
SRC_PATH        := $(DIST_PATH)/$(SRC_DIR)
SRC_FRAG_PATH   := ./$(OBJ_DIR)/$(HOST)/$(SRC_DIR)

# These variables are used by the include statements as they recursively include
# one another. We initialize them to the distribution path and the base object
# path.
PARENT_SRC_PATH := $(DIST_PATH)
PARENT_PATH     := ./$(OBJ_DIR)/$(HOST)

# Recursively include all the makefile fragments.
#-include $(addsuffix /$(FRAGMENT_MK), $(FRAGMENT_SUB_DIRS))
-include $(addsuffix /$(FRAGMENT_MK), $(SRC_FRAG_PATH))

# Create a list of the makefile fragments.
MAKEFILE_FRAGMENTS := $(addsuffix /$(FRAGMENT_MK), $(FRAGMENT_DIR_PATHS))


# Detect whether we actually got any makefile fragments. If we didn't, then it
# is likely that the user has not yet generated them (via configure).
ifeq ($(strip $(MAKEFILE_FRAGMENTS)),)
MAKEFILE_FRAGMENTS_PRESENT := no
else
MAKEFILE_FRAGMENTS_PRESENT := yes
endif



#
# --- Compiler include path definitions ----------------------------------------
#

# A script (originating from BLIS) that creates a monolithic header file from
# many header files that are recursively #included from one another.
BUILD_PATH := $(DIST_PATH)/build
FLATTEN_H  := $(PYTHON) $(BUILD_PATH)/flatten-headers.py

# The path to the main header files.
FLAF2C_H_SRC_PATH := $(DIST_PATH)/src/base/flamec/include/FLA_f2c.h
FLAME_H_SRC_PATH  := $(DIST_PATH)/src/base/flamec/include/FLAME.h
BLIS1_H_SRC_PATH  := $(DIST_PATH)/src/base/flamec/blis/include/blis1.h

# Construct the path to what will be the intermediate flattened/monolithic
# header files.
FLAME_H       := FLAME.h
FLAME_H_FLAT  := $(BASE_INC_PATH)/$(FLAME_H)

BLIS1_H       := blis1.h
BLIS1_H_FLAT  := $(BASE_INC_PATH)/$(BLIS1_H)

FLAF2C_H      := FLA_f2c.h
FLAF2C_H_FLAT := $(BASE_INC_PATH)/$(FLAF2C_H)

#Define path of CPP Template header files
CPP_TEMPLATE_H_PATH := ./$(SRC_DIR)/src_cpp

# Expand the fragment paths that contain .h files to attain the set of header
# files present in all fragment paths.
MK_HEADER_FILES := $(foreach frag_path, $(FRAGMENT_DIR_PATHS), $(wildcard $(frag_path)/*.h))
MK_HEADER_FILES += $(BASE_CONFIG_PATH)/FLA_config.h

# Strip the leading, internal, and trailing whitespace from our list of header
# files. This makes the "make install-headers" much more readable.
MK_HEADER_FILES := $(strip $(MK_HEADER_FILES))

# Expand the fragment paths that contain .h files, and take the first expansion.
# Then, strip the header filename to leave the path to each header location.
# Notice this process even weeds out duplicates! Add the config directory manually
# since it contains FLA_config.h.
MK_HEADER_DIR_PATHS := $(dir $(foreach frag_path, $(FRAGMENT_DIR_PATHS), \
                                       $(firstword $(wildcard $(frag_path)/*.h))))
MK_HEADER_DIR_PATHS += $(BASE_CONFIG_PATH)

# Define a list of headers to flatten. We have to flatten blis1.h and FLA_f2c.h
# because a few files #include only those files, but they aren't needed after
# libflame is compiled.
HEADERS_TO_FLATTEN := $(FLAME_H_FLAT) $(BLIS1_H_FLAT) $(FLAF2C_H_FLAT)

# Define a list of headers to install and their installation path. The default
# is to only install FLAME.h.
HEADERS_TO_INSTALL := $(FLAME_H_FLAT)
HEADERS_TO_INSTALL += $(CPP_TEMPLATE_H_PATH)/*.hh
#LAPACKE headers for cpp interface
LAPACKE_HEADERS    := $(SRC_DIR)/$(LAPACKE_DIR)/LAPACKE/include/lapacke.h
LAPACKE_HEADERS    += $(SRC_DIR)/$(LAPACKE_DIR)/LAPACKE/include/lapacke_mangling.h
LAPACKE_HEADERS    += $(SRC_DIR)/$(LAPACKE_DIR)/LAPACKE/include/lapack.h

HEADERS_TO_INSTALL += $(LAPACKE_HEADERS)
HEADERS_INST       := $(addprefix $(MK_INCL_DIR_INST)/, $(notdir $(HEADERS_TO_INSTALL)))

# Add -I to each header path so we can specify our include search paths to the
# C and Fortran compilers. NOTE: This is primarily for access to the monolithic
# (flattened) FLAME.h file, athough a few other files are placed there too (see
# above).
#INCLUDE_PATHS   := $(strip $(patsubst %, -I%, $(MK_HEADER_DIR_PATHS)))
INCLUDE_PATHS   := $(strip $(patsubst %, -I%, $(BASE_INC_PATH)))

# When lapack2flame is enabled, we need to add a -I option for the directory
# in which the lapack2flame headers reside.
ifeq ($(FLA_ENABLE_LAPACK2FLAME),yes)
L2F_FRAG_DIR_PATHS   := $(filter $(DIST_PATH)/src/map/lapack2flamec%,$(FRAGMENT_DIR_PATHS))
L2F_HEADER_DIR_PATHS := $(dir $(foreach frag_path, $(L2F_FRAG_DIR_PATHS), \
                                       $(firstword $(wildcard $(frag_path)/*.h))))
INCLUDE_PATHS   += $(strip $(patsubst %, -I%, $(L2F_HEADER_DIR_PATHS)))
endif

# Add the include flags determined above to various compiler flags variables.
CFLAGS          := $(CFLAGS) $(INCLUDE_PATHS)
CFLAGS_NOOPT    := $(CFLAGS_NOOPT) $(INCLUDE_PATHS)
CPPFLAGS        := $(CPPFLAGS) $(INCLUDE_PATHS)
FFLAGS          := $(FFLAGS) $(INCLUDE_PATHS)



#
# --- Library object definitions -----------------------------------------------
#

# Convert source file paths to object file paths by replaying the base source
# directory with the base object directory, and also replacing the source file
# suffix (ie: '.c' or '.f') with '.o'.

MK_FLABLAS_F2C_OBJS                   := $(patsubst $(SRC_PATH)/%.c, $(BASE_OBJ_PATH)/%.o, \
                                                    $(filter %.c, $(MK_FLABLAS_F2C_SRC)))

MK_BASE_FLAMEC_OBJS                   := $(patsubst $(SRC_PATH)/%.c, $(BASE_OBJ_PATH)/%.o, \
                                                    $(filter %.c, $(MK_BASE_FLAMEC_SRC)))

MK_BLAS_FLAMEC_OBJS                   := $(patsubst $(SRC_PATH)/%.c, $(BASE_OBJ_PATH)/%.o, \
                                                    $(filter %.c, $(MK_BLAS_FLAMEC_SRC)))

MK_LAPACK_FLAMEC_OBJS                 := $(patsubst $(SRC_PATH)/%.c, $(BASE_OBJ_PATH)/%.o, \
                                                    $(filter %.c, $(MK_LAPACK_FLAMEC_SRC)))

MK_MAP_LAPACK2FLAMEC_OBJS             := $(patsubst $(SRC_PATH)/%.c, $(BASE_OBJ_PATH)/%.o, \
                                                    $(filter %.c, $(MK_MAP_LAPACK2FLAMEC_SRC)))

MK_MAP_LAPACK2FLAMEC_F2C_OBJS         := $(patsubst $(SRC_PATH)/%.c, $(BASE_OBJ_PATH)/%.o, \
                                                    $(filter %.c, $(MK_MAP_LAPACK2FLAMEC_F2C_SRC)))

MK_MAP_LAPACK2FLAMEC_F2C_FLAMEC_OBJS  := $(patsubst $(SRC_PATH)/%.c, $(BASE_OBJ_PATH)/%.o, \
                                                    $(filter %.c, $(MK_MAP_LAPACK2FLAMEC_F2C_FLAMEC_SRC)))

MK_MAP_LAPACK2FLAMEC_F2C_INSTALL_OBJS := $(patsubst $(SRC_PATH)/%.c, $(BASE_OBJ_PATH)/%.o, \
                                                    $(filter %.c, $(MK_MAP_LAPACK2FLAMEC_F2C_INSTALL_SRC)))

# Combine the base, blas, and lapack libraries.
MK_ALL_FLAMEC_OBJS        := $(MK_BASE_FLAMEC_OBJS) \
                             $(MK_BLAS_FLAMEC_OBJS) \
                             $(MK_LAPACK_FLAMEC_OBJS)

# Prepend the flablas source code files, if requested

# LAPACK
ifeq ($(FLA_ENABLE_LAPACK2FLAME),yes)
MK_ALL_FLAMEC_OBJS        := $(MK_MAP_LAPACK2FLAMEC_OBJS) \
                             $(MK_MAP_LAPACK2FLAMEC_F2C_OBJS) \
                             $(MK_MAP_LAPACK2FLAMEC_F2C_FLAMEC_OBJS) \
                             $(MK_MAP_LAPACK2FLAMEC_F2C_INSTALL_OBJS) \
                             $(MK_ALL_FLAMEC_OBJS)
endif

# BLAS
ifeq ($(FLA_ENABLE_BUILTIN_BLAS),yes)
ifeq ($(FLA_ENABLE_LAPACK2FLAME),no)
MK_FLABLAS_F2C_OBJS       := $(MK_MAP_LAPACK2FLAMEC_F2C_INSTALL_OBJS) \
                             $(MK_FLABLAS_F2C_OBJS)
endif
MK_ALL_FLAMEC_OBJS        := $(MK_FLABLAS_F2C_OBJS) \
                             $(MK_ALL_FLAMEC_OBJS)
endif

### Kyungjoo 2015.10.21
#AR_CHUNK_SIZE=4096
AR_CHUNK_SIZE=1024

#
# --- Targets/rules ------------------------------------------------------------
#

# --- Primary targets ---

all: libs

libs: libflame

check: checklibflame

install: libs install-libs install-lib-symlinks install-headers

uninstall: uninstall-libs uninstall-lib-symlinks uninstall-headers

clean: cleanh cleanlib


# --- Environment check rules ---

check-env: check-env-config check-env-fragments

check-env-config:
ifeq ($(CONFIG_MK_PRESENT),no)
	$(error Cannot proceed: config.mk not detected! Run configure first)
endif

check-env-fragments: check-env-config
ifeq ($(MAKEFILE_FRAGMENTS_PRESENT),no)
	$(error Cannot proceed: makefile fragments not detected! Run configure first)
endif

# --- Cosolidated header creation ---

flat-headers: check-env $(FLAME_H_FLAT) $(BLIS1_H_FLAT) $(FLAF2C_H_FLAT)

# Consolidated FLAME.h header creation
$(FLAME_H_FLAT): $(MK_HEADER_FILES)
ifeq ($(ENABLE_VERBOSE),yes)
	$(FLATTEN_H) -c -v1 $(FLAME_H_SRC_PATH) $@ $(BASE_INC_PATH) "$(MK_HEADER_DIR_PATHS)"
else
	@echo -n "Generating monolithic $(@)"
	@$(FLATTEN_H) -c -v1 $(FLAME_H_SRC_PATH) $@ $(BASE_INC_PATH) "$(MK_HEADER_DIR_PATHS)"
	@echo "Generated monolithic $@"
endif

# Consolidated blis1.h header creation
$(BLIS1_H_FLAT): $(MK_HEADER_FILES) $(FLAME_H_FLAT)
ifeq ($(ENABLE_VERBOSE),yes)
	$(FLATTEN_H) -c -v1 $(BLIS1_H_SRC_PATH) $@ $(BASE_INC_PATH) "$(MK_HEADER_DIR_PATHS)"
else
	@echo -n "Generating monolithic $(@)"
	@$(FLATTEN_H) -c -v1 $(BLIS1_H_SRC_PATH) $@ $(BASE_INC_PATH) "$(MK_HEADER_DIR_PATHS)"
	@echo "Generated monolithic $@"
endif

# Consolidated FLA_f2c.h header creation
# NOTE: This file doesn't actually need any inlining, but we do need to
# "install" it to $(BASE_INC_PATH), so we opt to go through the same
# motions as for FLAME.h and FLA_f2c.h.
$(FLAF2C_H_FLAT): $(MK_HEADER_FILES) $(FLAME_H_FLAT) $(BLIS1_H_FLAT)
ifeq ($(ENABLE_VERBOSE),yes)
	$(FLATTEN_H) -c -v1 $(FLAF2C_H_SRC_PATH) $@ $(BASE_INC_PATH) "$(MK_HEADER_DIR_PATHS)"
else
	@echo -n "Generating monolithic $(@)"
	@$(FLATTEN_H) -c -v1 $(FLAF2C_H_SRC_PATH) $@ $(BASE_INC_PATH) "$(MK_HEADER_DIR_PATHS)"
	@echo "Generated monolithic $@"
endif


# --- Special source code / object code rules ---

FLA_SLAMCH=base/flamec/util/lapack/mch/fla_slamch
$(BASE_OBJ_PATH)/$(FLA_SLAMCH).o: $(SRC_PATH)/$(FLA_SLAMCH).c $(CONFIG_MK_FILE) $(HEADERS_TO_FLATTEN)
ifeq ($(ENABLE_VERBOSE),yes)
	$(CC) $(CFLAGS_NOOPT) -c $< -o $@
else
	@echo "Compiling $<"
	@$(CC) $(CFLAGS_NOOPT) -c $< -o $@
endif
ifeq ($(OS_NAME),Darwin)
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
	@echo $@ >> $(AR_OBJ_LIST_FILE)
endif
endif

FLA_DLAMCH=base/flamec/util/lapack/mch/fla_dlamch
$(BASE_OBJ_PATH)/$(FLA_DLAMCH).o: $(SRC_PATH)/$(FLA_DLAMCH).c $(CONFIG_MK_FILE) $(HEADERS_TO_FLATTEN)
ifeq ($(ENABLE_VERBOSE),yes)
	$(CC) $(CFLAGS_NOOPT) -c $< -o $@
else
	@echo "Compiling $<"
	@$(CC) $(CFLAGS_NOOPT) -c $< -o $@
endif
ifeq ($(OS_NAME),Darwin)
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
	@echo $@ >> $(AR_OBJ_LIST_FILE)
endif
endif

# --- General source code / object code rules ---

# Default compilation rules
$(BASE_OBJ_PATH)/%.o: $(SRC_PATH)/%.c $(CONFIG_MK_FILE) $(HEADERS_TO_FLATTEN)
ifeq ($(ENABLE_VERBOSE),yes)
	$(CC) $(CFLAGS) -c $< -o $@
else
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -c $< -o $@
endif
ifeq ($(OS_NAME),Darwin)
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
	@echo $@ >> $(AR_OBJ_LIST_FILE)
endif
endif

# --- All-purpose library rule (static and shared) ---

libflame: check-env $(MK_LIBS)


# --- Static library archiver rules ---
$(LAPACKE_A_PATH):
ifeq ($(ENABLE_VERBOSE),yes)
	$(MAKE) -e -C $(SRC_DIR)/$(LAPACKE_DIR)/LAPACKE
else
	@echo -n "Generating LAPACKE library"
	$(MAKE) -e -C $(SRC_DIR)/$(LAPACKE_DIR)/LAPACKE
	@echo "Generated LAPACKE library"
endif
	@echo $(LAPACKE_S_OBJS_PATH) >> $(AR_OBJ_LIST_FILE)
	@echo $(LAPACKE_U_OBJS_PATH) >> $(AR_OBJ_LIST_FILE)
$(AOCLDTL_A_PATH):
ifeq ($(ENABLE_VERBOSE),yes)
	$(MAKE) -e -C $(SRC_DIR)/$(AOCLDTL_DIR)
else
	@echo -n "Generating AOCLDTL library"
	$(MAKE) -e -C $(SRC_DIR)/$(AOCLDTL_DIR)
	@echo "Generated AOCLDTL library"
endif

$(LIBFLAME_A_PATH): $(MK_ALL_FLAMEC_OBJS)
ifeq ($(ENABLE_VERBOSE),yes)
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
ifeq ($(OS_NAME),Darwin)
### Kyungjoo 2015.10.21
	$(CAT) $(AR_OBJ_LIST_FILE) | xargs -n$(AR_CHUNK_SIZE) $(AR) $(ARFLAGS) $@
### Previous hack (works on linux, not on osx; osx's ar does not support @file)
#	echo $(ARFLAGS) $@ > $(AR_ARG_LIST_FILE)
#	$(CAT) $(AR_OBJ_LIST_FILE) >> $(AR_ARG_LIST_FILE)
#	$(AR) @$(AR_ARG_LIST_FILE)
else
	$(file > $@.in,$^)
	$(AR) $(ARFLAGS) $@ @$@.in
	$(RM_F) $@.in
endif 
else
#	NOTE: Can't use $^ automatic variable as long as $(AR_OBJ_LIST_FILE) is in
#	the list of prerequisites.
	$(AR) $(ARFLAGS) $@ $^
endif
	$(RANLIB) $@
#	$(MKDIR) include_local
#	cp -f $(MK_HEADER_FILES) include_local
else
	@echo "Archiving $@"
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
ifeq ($(OS_NAME),Darwin)
### Kyungjoo 2015.10.21
	@$(CAT) $(AR_OBJ_LIST_FILE) | xargs -n$(AR_CHUNK_SIZE) $(AR) $(ARFLAGS) $@
### Previous hack (works on linux, not on osx; osx's ar does not support @file)
#	@echo $(ARFLAGS) $@ > $(AR_ARG_LIST_FILE)
#	@$(CAT) $(AR_OBJ_LIST_FILE) >> $(AR_ARG_LIST_FILE)
#	@$(AR) @$(AR_ARG_LIST_FILE)
else
	@$(file > $@.in,$^)
	@$(AR) $(ARFLAGS) $@ @$@.in
	@$(RM_F) $@.in
endif
else
#	NOTE: Can't use $^ automatic variable as long as $(AR_OBJ_LIST_FILE) is in
#	the list of prerequisites.
	@$(AR) $(ARFLAGS) $@ $^
endif
	@$(RANLIB) $@
#	@$(MKDIR) include_local
#	@cp -f $(MK_HEADER_FILES) include_local
endif


# --- Shared library linker rules ---

$(LIBFLAME_SO_PATH): $(MK_ALL_FLAMEC_OBJS)
ifeq ($(ENABLE_VERBOSE),yes)
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
ifeq ($(OS_NAME),Darwin)
	$(CAT) $(AR_OBJ_LIST_FILE) | xargs -n$(AR_CHUNK_SIZE) $(AR) $(ARFLAGS) $(LIBFLAME_A)
	$(LINKER) $(SOFLAGS) -o $@ -Wl,-force_load,$(LIBFLAME_A) $(LDFLAGS)
else
	$(file > $@.in,$^)
	$(LINKER) $(SOFLAGS) -o $(LIBFLAME_SO_OUTPUT_NAME) @$@.in $(LDFLAGS)
	$(RM_F) $@.in
endif
else
#	NOTE: Can't use $^ automatic variable as long as $(AR_OBJ_LIST_FILE) is in
#	the list of prerequisites.
	$(LINKER) $(SOFLAGS) -o $(LIBFLAME_SO_OUTPUT_NAME) $^ $(LDFLAGS)
endif
else
	@echo "Dynamically linking $@"
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
ifeq ($(OS_NAME),Darwin)
	@$(CAT) $(AR_OBJ_LIST_FILE) | xargs -n$(AR_CHUNK_SIZE) $(AR) $(ARFLAGS) $(LIBFLAME_A)
	@$(LINKER) $(SOFLAGS) -o $@ -Wl,-force_load,$(LIBFLAME_A) $(LDFLAGS)
else
	@$(file > $@.in,$^)
	@$(LINKER) $(SOFLAGS) -o $(LIBFLAME_SO_OUTPUT_NAME) @$@.in $(LDFLAGS)
	@$(RM_F) $@.in
endif
else
#	NOTE: Can't use $^ automatic variable as long as $(AR_OBJ_LIST_FILE) is in
#	the list of prerequisites.
	@$(LINKER) $(SOFLAGS) -o $(LIBFLAME_SO_OUTPUT_NAME) $^ $(LDFLAGS)
endif
endif

# Local symlink for shared library.
# NOTE: We use a '.loc' suffix to avoid filename collisions in case this
# rule is executed concurrently with the install-lib-symlinks rule, which
# also creates symlinks in the current directory (before installing them).
$(LIBFLAME_SO_MAJ_PATH): $(LIBFLAME_SO_PATH)
ifeq ($(ENABLE_VERBOSE),yes)
	$(SYMLINK) $(<F) $(@F).loc
	$(MV) $(@F).loc $(BASE_LIB_PATH)/$(@F)
else # ifeq ($(ENABLE_VERBOSE),no)
	@echo "Creating symlink $@"
	@$(SYMLINK) $(<F) $(@F).loc
	@$(MV) $(@F).loc $(BASE_LIB_PATH)/$(@F)
endif

# Original implementation of the rule above.
# FGVZ: This rule has been observed to not work on at least one system, where
# it appears the ".in" file is not fully written out, or written out at all,
# prior to the shared library link command being executed.
#$(MK_ALL_FLAMEC_DLL): $(MK_ALL_FLAMEC_OBJS)
#ifeq ($(ENABLE_VERBOSE),yes)
#ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
#	$(file > $@.in,$^)
#	$(LINKER) -shared -Wl,-soname,libflame.so $(LDFLAGS) -o $@ @$@.in
#	$(RM) $@.in
#else
#	$(LINKER) -shared -Wl,-soname,libflame.so $(LDFLAGS) -o $@ $^
#endif
#else
#	@echo "Dynamically linking $@"
#ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
#	@$(file > $@.in,$^)
#	@$(LINKER) -shared -Wl,-soname,libflame.so $(LDFLAGS) -o $@ @$@.in
#	@$(RM) $@.in
#else
#	@$(LINKER) -shared -Wl,-soname,libflame.so $(LDFLAGS) -o $@ $^
#endif
#endif




# --- Test suite rules ---

# The test binary executable filename.
#
#   $(TEST_WRAPPER) ./$(TEST_BIN) ... > $(TEST_OUT_FILE)
#
# "TEST_WRAPPER can be used to set additional environment variables
#  such as LD_LIBRARY_PATH, numactl and others"
#
TEST_BIN      := test_$(LIBFLAME).x
TEST_WRAPPER  ?=
TEST_CHECK    := check-flametest.sh

# The location of the script that checks the BLIS testsuite output.
TEST_CHECK_PATH    := $(DIST_PATH)/$(TEST_DIR)/$(TEST_CHECK)

test-bin: check-env $(TEST_BIN)

TEST_OUT_FILE := output.test

# Check the results of the LIBLFLAME tests.
$(TEST_BIN):
ifeq ($(ENABLE_VERBOSE),yes)
	$(MAKE) -C $(TEST_DIR)
else
	@$(MAKE) -C $(TEST_DIR)
endif

test-run: test-bin
ifeq ($(ENABLE_VERBOSE),yes)
	cd $(TEST_DIR) && $(TEST_WRAPPER) ./$(TEST_BIN) > $(TEST_OUT_FILE)
else
	@echo "Running $(TEST_BIN) with output redirected to '$(TEST_OUT_FILE)'"
	@echo "Please wait for tests to complete..."
	@cd $(TEST_DIR) && $(TEST_WRAPPER) ./$(TEST_BIN) > $(TEST_OUT_FILE)
endif

checklibflame:test-run
ifeq ($(ENABLE_VERBOSE),yes)
	-sh $(TEST_CHECK_PATH) $(TEST_DIR)/$(TEST_OUT_FILE)
else
	-@sh $(TEST_CHECK_PATH) $(TEST_DIR)/$(TEST_OUT_FILE)
endif


# --- Install rules ---


# --- Install header rules ---

install-headers: check-env $(HEADERS_INST)

$(HEADERS_INST): $(HEADERS_TO_INSTALL) $(CONFIG_MK_FILE)
ifeq ($(ENABLE_VERBOSE),yes)
	$(MKDIR) $(MK_INCL_DIR_INST)
	$(INSTALL) -m 0644 $(HEADERS_TO_INSTALL) $(MK_INCL_DIR_INST)
else
	@$(MKDIR) $(MK_INCL_DIR_INST)
	@echo "Installing $(notdir $(HEADERS_TO_INSTALL)) into $(MK_INCL_DIR_INST)/"
	@$(INSTALL) -m 0644 $(HEADERS_TO_INSTALL) $(MK_INCL_DIR_INST)
endif


# Run CPP Tests
checkcpp:
	$(MAKE) -e -C $(CPP_TEST_DIR)

# --- Install library rules ---

install-libs: check-env $(MK_LIBS_INST)

# Install static library.
$(INSTALL_LIBDIR)/%.a: $(BASE_LIB_PATH)/%.a $(CONFIG_MK_FILE)
ifeq ($(ENABLE_VERBOSE),yes)
	$(MKDIR) $(@D)
	$(INSTALL) -m 0644 $< $@
	$(INSTALL) -m 0644 $(LAPACKE_A_PATH) $(INSTALL_LIBDIR)/$(LAPACKE_A)
	$(INSTALL) -m 0644 $(AOCLDTL_A_PATH) $(INSTALL_LIBDIR)/$(AOCLDTL_A)
else
	@echo "Installing $(@F) into $(INSTALL_LIBDIR)/"
	@$(MKDIR) $(@D)
	@$(INSTALL) -m 0644 $< $@
	@$(INSTALL) -m 0644 $(LAPACKE_A_PATH) $(INSTALL_LIBDIR)/$(LAPACKE_A)
	@$(INSTALL) -m 0644 $(AOCLDTL_A_PATH) $(INSTALL_LIBDIR)/$(AOCLDTL_A)
endif

# Install shared library.
$(INSTALL_LIBDIR)/%.$(LIBFLAME_SO_MMB_EXT): $(BASE_LIB_PATH)/%.$(SHLIB_EXT) $(CONFIG_MK_FILE)
ifeq ($(ENABLE_VERBOSE),yes)
	$(MKDIR) $(@D)
	$(INSTALL) -m 0755 $< $@
	$(INSTALL) -m 0644 $(AOCLDTL_SO_PATH) $(INSTALL_LIBDIR)/$(AOCLDTL_SO)
else
	@echo "Installing $(@F) into $(INSTALL_LIBDIR)/"
	@$(MKDIR) $(@D)
	@$(INSTALL) -m 0755 $< $@
	@$(INSTALL) -m 0644 $(AOCLDTL_SO_PATH) $(INSTALL_LIBDIR)/$(AOCLDTL_SO)
endif


# --- Install-symlinks rules ---

install-lib-symlinks: check-env $(MK_LIBS_SYML)

$(INSTALL_LIBDIR)/%.$(SHLIB_EXT): $(INSTALL_LIBDIR)/%.$(LIBFLAME_SO_MMB_EXT)
ifeq ($(ENABLE_VERBOSE),yes)
	$(SYMLINK) $(<F) $(@F)
	$(MV) $(@F) $(INSTALL_LIBDIR)/
else
	@echo "Installing symlink $(@F) into $(INSTALL_LIBDIR)/"
	@$(SYMLINK) $(<F) $(@F)
	@$(MV) $(@F) $(INSTALL_LIBDIR)/
endif

# Install shared library symlink containing only .so major version.
$(INSTALL_LIBDIR)/%.$(LIBFLAME_SO_MAJ_EXT): $(INSTALL_LIBDIR)/%.$(LIBFLAME_SO_MMB_EXT)
ifeq ($(ENABLE_VERBOSE),yes)
	$(SYMLINK) $(<F) $(@F)
	$(MV) $(@F) $(INSTALL_LIBDIR)/
else
	@echo "Installing symlink $(@F) into $(INSTALL_LIBDIR)/"
	@$(SYMLINK) $(<F) $(@F)
	@$(MV) $(@F) $(INSTALL_LIBDIR)/
endif


# --- Clean rules ---

cleanmk:
ifeq ($(IS_CONFIGURED),yes)
ifeq ($(ENABLE_VERBOSE),yes)
	- $(FIND) $(SRC_FRAG_PATH) -name "$(FRAGMENT_MK)" | $(XARGS) $(RM_F)
else
	@echo "Removing makefile fragments from $(SRC_FRAG_PATH)"
	@$(FIND) $(SRC_FRAG_PATH) -name "$(FRAGMENT_MK)" | $(XARGS) $(RM_F)
endif
endif

cleanh:
ifeq ($(IS_CONFIGURED),yes)
ifeq ($(ENABLE_VERBOSE),yes)
	- $(RM_F) $(HEADERS_TO_FLATTEN)
else
	@echo "Removing flattened header files from $(BASE_INC_PATH)"
	@$(RM_F) $(HEADERS_TO_FLATTEN)
endif
endif

cleanlib:
ifeq ($(IS_CONFIGURED),yes)
ifeq ($(ENABLE_VERBOSE),yes)
	- $(FIND) $(BASE_OBJ_PATH) -name "*.o" | $(XARGS) $(RM_F)
	- $(FIND) $(LAPACKE_S_OBJS_PATH) -name "*.o" | $(XARGS) $(RM_F)
	- $(FIND) $(LAPACKE_U_OBJS_PATH) -name "*.o" | $(XARGS) $(RM_F)
	- $(RM_F) $(LAPACKE_A_PATH)
	- $(RM_F) $(AOCLDTL_A_PATH)
	- $(RM_F) $(AOCLDTL_SO_PATH)
	- $(RM_F) $(AOCLDTL_obj_PATH)
	- $(RM_F) $(AOCLDTL_gch_PATH)
	- $(RM_F) $(BASE_LIB_PATH)/*
else
	@echo "Removing object files from $(BASE_OBJ_PATH)"
	@$(FIND) $(BASE_OBJ_PATH) -name "*.o" | $(XARGS) $(RM_F)
	@$(FIND) $(LAPACKE_S_OBJS_PATH) -name "*.o" | $(XARGS) $(RM_F)
	@$(FIND) $(LAPACKE_U_OBJS_PATH) -name "*.o" | $(XARGS) $(RM_F)
	@echo "Removing libraries from $(BASE_LIB_PATH)"
	@$(RM_F) $(LAPACKE_A_PATH)
	@$(RM_F) $(AOCLDTL_A_PATH)
	@$(RM_F) $(AOCLDTL_SO_PATH)
	@$(RM_F) $(AOCLDTL_obj_PATH)
	@$(RM_F) $(AOCLDTL_gch_PATH)
	@$(RM_F) $(BASE_LIB_PATH)/*
endif
endif

distclean: cleanmk cleanh cleanlib
ifeq ($(IS_CONFIGURED),yes)
ifeq ($(ENABLE_VERBOSE),yes)
ifeq ($(OS_NAME),Darwin)
	- $(RM_F) $(AR_OBJ_LIST_FILE)
endif 
	- $(RM_RF) $(CONFIG_DIR)
	- $(RM_RF) $(OBJ_DIR)
	- $(RM_RF) $(LIB_DIR)
	- $(RM_RF) $(INC_DIR)
	- $(RM_RF) config.log
	- $(RM_RF) aclocal.m4
	- $(RM_RF) autom4te.cache
	- $(RM_RF) config.status
	- $(RM_RF) config.sys_type
	- $(RM_RF) config.dist_path
else
ifeq ($(OS_NAME),Darwin)
	@echo "Removing $(AR_OBJ_LIST_FILE)"
	@$(RM_F) $(AR_OBJ_LIST_FILE)
endif 
	@echo "Removing $(CONFIG_DIR)"
	@$(RM_RF) $(CONFIG_DIR)
	@echo "Removing $(OBJ_DIR)"
	@$(RM_RF) $(OBJ_DIR)
	@echo "Removing $(LIB_DIR)"
	@$(RM_RF) $(LIB_DIR)
	@echo "Removing $(INC_DIR)"
	@$(RM_RF) $(INC_DIR)
	@echo "Removing intermediate configure files"
	@$(RM_RF) config.log
	@$(RM_RF) aclocal.m4
	@$(RM_RF) autom4te.cache
	@$(RM_RF) config.status
	@$(RM_RF) config.sys_type
	@$(RM_RF) config.dist_path
endif
endif

cleanleaves:
ifeq ($(IS_CONFIGURED),yes)
ifeq ($(ENABLE_VERBOSE),yes)
	- $(FIND) $(BASE_OBJ_PATH) -name "*.[osx]" | $(XARGS) $(RM_F)
else
	@echo "Removing leaf-level build objects from source tree"
	@$(FIND) $(BASE_OBJ_PATH) -name "*.[osx]" | $(XARGS) $(RM_F)
endif
endif

cleancpptest:
ifeq ($(ENABLE_VERBOSE),yes)
	- $(MAKE) -C $(CPP_TEST_DIR) clean
else
	@echo "Clean up CPP tests"
	@@$(MAKE) -C $(CPP_TEST_DIR) clean
endif

# --- Uninstall rules ---

# NOTE: We can't write these uninstall rules directly in terms of targets
# $(MK_LIBS_INST) and $(MK_INCL_DIR_INST)
# because those targets are already defined in terms of rules that *build*
# those products.

uninstall-libs: check-env
ifeq ($(ENABLE_VERBOSE),yes)
	- $(RM_F) $(MK_LIBS_INST)
else
	@echo "Uninstalling libraries $(notdir $(MK_LIBS_INST)) from $(dir $(firstword $(MK_LIBS_INST)))."
	@- $(RM_F) $(MK_LIBS_INST)
endif

uninstall-lib-symlinks: check-env
ifeq ($(ENABLE_VERBOSE),yes)
	- $(RM_F) $(MK_LIBS_SYML)
else
	@echo "Uninstalling symlinks $(notdir $(MK_LIBS_SYML)) from $(dir $(firstword $(MK_LIBS_SYML)))."
	@- $(RM_F) $(MK_LIBS_SYML)
endif

uninstall-headers: check-env
ifeq ($(ENABLE_VERBOSE),yes)
	- $(RM_F) $(HEADERS_INST)
else
	@echo "Uninstalling headers '$(notdir $(HEADERS_TO_INSTALL))' from $(MK_INCL_DIR_INST)."
#	@- $(RM_F) $(addprefix $(MK_INCL_DIR_INST)/, $(notdir $(HEADERS_TO_INSTALL)))
	@- $(RM_F) $(HEADERS_INST)
endif

