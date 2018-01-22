
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

.PHONY: all libs docs test install clean \
        check check-config check-fragments \
        install-libs install-headers install-docs \
        install-lib-symlinks install-header-symlinks \
        install-without-symlinks \
        cleanmost distclean cleanmk cleanleaves \
        send-thanks



#
# --- Makefile fragment initialization -----------------------------------------
#

# Makefile fragment name.
FRAGMENT_MK     := .fragment.mk

# Locations of important files.
BUILD_DIR       := ./build
CONFIG_DIR      := ./config
SRC_DIR         := ./src
OBJ_DIR         := ./obj
LIB_DIR         := ./lib
INCLUDE_LOCAL   := ./include_local

# The host string will uniquely identify the current host (as much
# as is reasonable) for purposes of separating the configure products and
# object files of one architecture from another.
#HOST            := $(shell sh $(BUILD_DIR)/ac-utils/config.guess)
HOST            := $(shell cat config.sys_type)

# Use the system type to name the config, object, and library directories.
# These directories are special in that they will contain products specific
# to this particular architecture.
BASE_CONFIG_DIR := $(CONFIG_DIR)/$(HOST)
BASE_OBJ_DIR    := $(OBJ_DIR)/$(HOST)
BASE_LIB_DIR    := $(LIB_DIR)/$(HOST)



#
# --- Include architecture-specific variable definitions ----------------------
#

# Pathnames to makefile fragment containing lots of various definitions.
CONFIG_MK_FRAGMENT   := $(BASE_CONFIG_DIR)/config.mk

# Include the definitions in the config makefile fragment.
-include $(CONFIG_MK_FRAGMENT)

# Detect whether we actually got the config makefile fragment. If we didn't,
# then it is likely that the user has not yet generated it (via configure).
ifeq ($(strip $(CONFIG_MK_INCLUDED)),yes)
CONFIG_MK_PRESENT := yes
else
CONFIG_MK_PRESENT := no
endif



#
# --- Main target variable definitions ----------------------------------------
#

# Construct the architecture-version string, which will be used to name the
# libraries upon installation.
VERSION                        := $(shell cat version)
ARCH_VERS                      := $(ARCH)-$(VERSION)

# --- Library names ---
ALL_FLAMEC_LIB_NAME            := libflame.a
ALL_FLAMEC_DLL_NAME            := libflame.so

# --- FLAME/C variable names ---
MK_ALL_FLAMEC_LIB                     := $(BASE_LIB_DIR)/$(ALL_FLAMEC_LIB_NAME)
MK_ALL_FLAMEC_DLL                     := $(BASE_LIB_DIR)/$(ALL_FLAMEC_DLL_NAME)

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
MK_FLAMEC_LIBS                    := $(MK_ALL_FLAMEC_LIB)
MK_FLAMEC_LIBS_INST               := $(patsubst $(BASE_LIB_DIR)/%.a, \
                                                $(INSTALL_PREFIX)/lib/%.a, \
                                                $(MK_FLAMEC_LIBS))
MK_FLAMEC_LIBS_INST_W_ARCH_VERS   := $(patsubst $(BASE_LIB_DIR)/%.a, \
                                                $(INSTALL_PREFIX)/lib/%-$(ARCH_VERS).a, \
                                                $(MK_FLAMEC_LIBS))
#MK_FLAMEC_LIBS_INST_W_ARCH_ONLY   := $(patsubst $(BASE_LIB_DIR)/%.a, \
#                                                $(INSTALL_PREFIX)/lib/%-$(ARCH).a, \
#                                                $(MK_FLAMEC_LIBS))

# --- Define install target names for dynamic libraries ---
MK_FLAMEC_DLLS                    := $(MK_ALL_FLAMEC_DLL)
MK_FLAMEC_DLLS_INST               := $(patsubst $(BASE_LIB_DIR)/%.so, \
                                                $(INSTALL_PREFIX)/lib/%.so, \
                                                $(MK_FLAMEC_DLLS))
MK_FLAMEC_DLLS_INST_W_ARCH_VERS   := $(patsubst $(BASE_LIB_DIR)/%.so, \
                                                $(INSTALL_PREFIX)/lib/%-$(ARCH_VERS).so, \
                                                $(MK_FLAMEC_DLLS))
#MK_FLAMEC_DLLS_INST_W_ARCH_ONLY   := $(patsubst $(BASE_LIB_DIR)/%.so, \
#                                                $(INSTALL_PREFIX)/lib/%-$(ARCH).so, \
#                                                $(MK_FLAMEC_DLLS))

# --- Determine which libraries to build ---
MK_LIBS                           :=
MK_LIBS_INST                      :=
MK_LIBS_INST_W_ARCH_VERS          :=
#MK_LIBS_INST_W_ARCH_ONLY          :=

ifeq ($(FLA_ENABLE_STATIC_BUILD),yes)
MK_LIBS                           += $(MK_FLAMEC_LIBS)
MK_LIBS_INST                      += $(MK_FLAMEC_LIBS_INST)
MK_LIBS_INST_W_ARCH_VERS          += $(MK_FLAMEC_LIBS_INST_W_ARCH_VERS)
#MK_LIBS_INST_W_ARCH_ONLY          += $(MK_FLAMEC_LIBS_INST_W_ARCH_ONLY)
endif

ifeq ($(FLA_ENABLE_DYNAMIC_BUILD),yes)
MK_LIBS                           += $(MK_FLAMEC_DLLS)
MK_LIBS_INST                      += $(MK_FLAMEC_DLLS_INST)
MK_LIBS_INST_W_ARCH_VERS          += $(MK_FLAMEC_DLLS_INST_W_ARCH_VERS)
#MK_LIBS_INST_W_ARCH_ONLY          += $(MK_FLAMEC_DLLS_INST_W_ARCH_ONLY)
endif

# --- Set the include directory names ---
MK_INCL_DIR_INST                  := $(INSTALL_PREFIX)/include
MK_INCL_DIR_INST_W_ARCH_VERS      := $(INSTALL_PREFIX)/include-$(ARCH_VERS)
#MK_INCL_DIR_INST_W_ARCH_ONLY      := $(INSTALL_PREFIX)/include-$(ARCH)



#
# --- Test definitions --------------------------------------------------------
#





#
# --- Include sub-tree makefile fragments --------------------------------------
#

# Initialize our list of directory paths to makefile fragments with the empty list.
FRAGMENT_DIR_PATHS :=

# The only fragment subdirectory that we build from the top-level directory is the
# source directory.
FRAGMENT_SUB_DIRS  := $(SRC_DIR)

# This variable is used by the include statements as they recursively include one
# another. We initialize it to the current directory.
PARENT_PATH        := .

# Recursively include all the makefile fragments.
-include $(addsuffix /$(FRAGMENT_MK), $(FRAGMENT_SUB_DIRS))

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

# Expand the fragment paths that contain .h files to attain the set of header
# files present in all fragment paths.
MK_HEADER_FILES := $(foreach frag_path, $(FRAGMENT_DIR_PATHS), $(wildcard $(frag_path)/*.h))
MK_HEADER_FILES += $(BASE_CONFIG_DIR)/FLA_config.h

# Strip the leading, internal, and trailing whitespace from our list of header
# files. This makes the "make install-headers" much more readable.
MK_HEADER_FILES := $(strip $(MK_HEADER_FILES))

# Expand the fragment paths that contain .h files, and take the first expansion.
# Then, strip the header filename to leave the path to each header location.
# Notice this process even weeds out duplicates! Add the config directory manually
# since it contains FLA_config.h.
MK_HEADER_DIR_PATHS := $(dir $(foreach frag_path, $(FRAGMENT_DIR_PATHS), \
                                       $(firstword $(wildcard $(frag_path)/*.h))))
MK_HEADER_DIR_PATHS += $(BASE_CONFIG_DIR)

# Add -I to each header path so we can specify our include search paths to the
# C and Fortran compilers.
INCLUDE_PATHS   := $(strip $(patsubst %, -I%, $(MK_HEADER_DIR_PATHS)))
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
MK_FLABLAS_F2C_OBJS                   := $(patsubst $(SRC_DIR)/%.c, $(BASE_OBJ_DIR)/%.o, \
                                                    $(filter %.c, $(MK_FLABLAS_F2C_SRC)))

MK_BASE_FLAMEC_OBJS                   := $(patsubst $(SRC_DIR)/%.c, $(BASE_OBJ_DIR)/%.o, \
                                                    $(filter %.c, $(MK_BASE_FLAMEC_SRC)))

MK_BLAS_FLAMEC_OBJS                   := $(patsubst $(SRC_DIR)/%.c, $(BASE_OBJ_DIR)/%.o, \
                                                    $(filter %.c, $(MK_BLAS_FLAMEC_SRC)))

MK_LAPACK_FLAMEC_OBJS                 := $(patsubst $(SRC_DIR)/%.c, $(BASE_OBJ_DIR)/%.o, \
                                                    $(filter %.c, $(MK_LAPACK_FLAMEC_SRC)))

MK_MAP_LAPACK2FLAMEC_OBJS             := $(patsubst $(SRC_DIR)/%.c, $(BASE_OBJ_DIR)/%.o, \
                                                    $(filter %.c, $(MK_MAP_LAPACK2FLAMEC_SRC)))

MK_MAP_LAPACK2FLAMEC_F2C_OBJS         := $(patsubst $(SRC_DIR)/%.c, $(BASE_OBJ_DIR)/%.o, \
                                                    $(filter %.c, $(MK_MAP_LAPACK2FLAMEC_F2C_SRC)))

MK_MAP_LAPACK2FLAMEC_F2C_FLAMEC_OBJS  := $(patsubst $(SRC_DIR)/%.c, $(BASE_OBJ_DIR)/%.o, \
                                                    $(filter %.c, $(MK_MAP_LAPACK2FLAMEC_F2C_FLAMEC_SRC)))

MK_MAP_LAPACK2FLAMEC_F2C_INSTALL_OBJS := $(patsubst $(SRC_DIR)/%.c, $(BASE_OBJ_DIR)/%.o, \
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
AR_CHUNK_SIZE=4096

#
# --- Targets/rules ------------------------------------------------------------
#

# --- Primary targets ---
all: libs

libs: check $(MK_LIBS)

test: check

install: libs install-libs install-headers \
         install-lib-symlinks install-header-symlinks

install-without-symlinks: libs install-libs install-headers

clean: cleanmost



# --- Environment check rules ---
check: check-config check-fragments

check-config:
ifeq ($(CONFIG_MK_PRESENT),no)
	$(error Cannot proceed: config.mk not detected! Run configure first)
endif

check-fragments: check-config
ifeq ($(MAKEFILE_FRAGMENTS_PRESENT),no)
	$(error Cannot proceed: makefile fragments not detected! Run configure first)
endif

# --- Special source code / object code rules ---

FLA_SLAMCH=base/flamec/util/lapack/mch/fla_slamch
$(BASE_OBJ_DIR)/$(FLA_SLAMCH).o: $(SRC_DIR)/$(FLA_SLAMCH).c $(CONFIG_MK_FRAGMENT)
ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
	$(CC) $(CFLAGS_NOOPT) -c $< -o $@
else
	@echo "Compiling $<"
	@$(CC) $(CFLAGS_NOOPT) -c $< -o $@
endif
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
	@echo $@ >> $(AR_OBJ_LIST_FILE)
endif

FLA_DLAMCH=base/flamec/util/lapack/mch/fla_dlamch
$(BASE_OBJ_DIR)/$(FLA_DLAMCH).o: $(SRC_DIR)/$(FLA_DLAMCH).c $(CONFIG_MK_FRAGMENT)
ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
	$(CC) $(CFLAGS_NOOPT) -c $< -o $@
else
	@echo "Compiling $<"
	@$(CC) $(CFLAGS_NOOPT) -c $< -o $@
endif
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
	@echo $@ >> $(AR_OBJ_LIST_FILE)
endif

# --- General source code / object code rules ---

# Default compilation rules
$(BASE_OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(CONFIG_MK_FRAGMENT)
ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
	$(CC) $(CFLAGS) -c $< -o $@
else
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -c $< -o $@
endif
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
	@echo $@ >> $(AR_OBJ_LIST_FILE)
endif



# --- Static library archiver rules for libflame ---
$(MK_ALL_FLAMEC_LIB): $(MK_ALL_FLAMEC_OBJS)

ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
### Kyungjoo 2015.10.21
	$(CAT) $(AR_OBJ_LIST_FILE) | xargs -n$(AR_CHUNK_SIZE) $(AR) $(ARFLAGS) $@
### Previous hack (works on linux, not on osx)
#	echo $(ARFLAGS) $@ > $(AR_ARG_LIST_FILE)
#	$(CAT) $(AR_OBJ_LIST_FILE) >> $(AR_ARG_LIST_FILE)
#	$(AR) @$(AR_ARG_LIST_FILE)
	$(RANLIB) $@
else
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@
endif
	mkdir -p include_local
	cp -f $(MK_HEADER_FILES) include_local
else
	@echo "Archiving $@"
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
### Kyungjoo 2015.10.21
	@$(CAT) $(AR_OBJ_LIST_FILE) | xargs -n$(AR_CHUNK_SIZE) $(AR) $(ARFLAGS) $@
### Previous hack (works on linux, not on osx)
#	@echo $(ARFLAGS) $@ > $(AR_ARG_LIST_FILE)
#	@$(CAT) $(AR_OBJ_LIST_FILE) >> $(AR_ARG_LIST_FILE)
#	@$(AR) @$(AR_ARG_LIST_FILE)
	@$(RANLIB) $@
else
	@$(AR) $(ARFLAGS) $@ $?
	@$(RANLIB) $@
endif
	@mkdir -p include_local
	@cp -f $(MK_HEADER_FILES) include_local
endif



# --- Dynamic library linker rules for libflame ---
$(MK_ALL_FLAMEC_DLL): $(MK_ALL_FLAMEC_OBJS)
ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
	$(LINKER) -shared $(LDFLAGS) -o $@ $?
else
	@echo "Dynamically linking $@"
ifeq ($(FLA_ENABLE_MAX_ARG_LIST_HACK),yes)
	@$(file > $@.in,$^)
	@$(LINKER) -shared $(LDFLAGS) -o $@ @$@.in
	@$(RM) $@.in
else
	@$(LINKER) -shared $(LDFLAGS) -o $@ $?
endif
endif



# --- Test suite rules ---




# --- Install rules ---
install-libs: check $(MK_LIBS_INST_W_ARCH_VERS)

install-headers: check $(MK_INCL_DIR_INST_W_ARCH_VERS)

$(MK_INCL_DIR_INST_W_ARCH_VERS): $(MK_HEADER_FILES)
ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
	$(INSTALL) -m 0755 -d $(@)
	$(INSTALL) -m 0644 $(MK_HEADER_FILES) $(@)
else
	@$(INSTALL) -m 0755 -d $(@)
	@echo "Installing C header files into $(@)"
	@$(INSTALL) -m 0644 $(MK_HEADER_FILES) $(@)
endif

$(INSTALL_PREFIX)/lib/%-$(ARCH_VERS).a: $(BASE_LIB_DIR)/%.a
ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
	$(INSTALL) -m 0755 -d $(@D)
	$(INSTALL) -m 0644 $< $@
else
	@echo "Installing $(@F) into $(INSTALL_PREFIX)/lib/"
	@$(INSTALL) -m 0755 -d $(@D)
	@$(INSTALL) -m 0644 $< $@
endif

$(INSTALL_PREFIX)/lib/%-$(ARCH_VERS).so: $(BASE_LIB_DIR)/%.so
ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
	$(INSTALL) -m 0755 -d $(@D)
	$(INSTALL) -m 0644 $< $@
else
	@echo "Installing $(@F) into $(INSTALL_PREFIX)/lib/"
	@$(INSTALL) -m 0755 -d $(@D)
	@$(INSTALL) -m 0644 $< $@
endif



# --- Install-symlinks rules ---
install-lib-symlinks: check-config $(MK_LIBS_INST)

install-header-symlinks: check-config $(MK_INCL_DIR_INST)

$(MK_INCL_DIR_INST): $(MK_INCL_DIR_INST_W_ARCH_VERS)
ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
	$(SYMLINK) $(<F) $(@F)
	$(MV) $(@F) $(INSTALL_PREFIX)
else
	@echo "Installing symlink $(@F) into $(INSTALL_PREFIX)/"
	@$(SYMLINK) $(<F) $(@F)
	@$(MV) $(@F) $(INSTALL_PREFIX)
endif

$(INSTALL_PREFIX)/lib/%.a: $(INSTALL_PREFIX)/lib/%-$(ARCH_VERS).a
ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
	$(SYMLINK) $(<F) $(@F)
	$(MV) $(@F) $(INSTALL_PREFIX)/lib/
else
	@echo "Installing symlink $(@F) into $(INSTALL_PREFIX)/lib/"
	@$(SYMLINK) $(<F) $(@F)
	@$(MV) $(@F) $(INSTALL_PREFIX)/lib/
endif

#$(INSTALL_PREFIX)/lib/%-$(ARCH).a: $(INSTALL_PREFIX)/lib/%-$(ARCH_VERS).a
#ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
#	$(SYMLINK) $(<F) $(@F)
#	$(MV) $(@F) $(INSTALL_PREFIX)/lib/
#else
#	@echo "Installing symlink $(@F) into $(INSTALL_PREFIX)/lib/"
#	@$(SYMLINK) $(<F) $(@F)
#	@$(MV) $(@F) $(INSTALL_PREFIX)/lib/
#endif

$(INSTALL_PREFIX)/lib/%.so: $(INSTALL_PREFIX)/lib/%-$(ARCH_VERS).so
ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
	$(SYMLINK) $(<F) $(@F)
	$(MV) $(@F) $(INSTALL_PREFIX)/lib/
else
	@echo "Installing symlink $(@F) into $(INSTALL_PREFIX)/lib/"
	@$(SYMLINK) $(<F) $(@F)
	@$(MV) $(@F) $(INSTALL_PREFIX)/lib/
endif

#$(INSTALL_PREFIX)/lib/%-$(ARCH).so: $(INSTALL_PREFIX)/lib/%-$(ARCH_VERS).so
#ifeq ($(FLA_ENABLE_VERBOSE_MAKE_OUTPUT),yes)
#	$(SYMLINK) $(<F) $(@F)
#	$(MV) $(@F) $(INSTALL_PREFIX)/lib/
#else
#	@echo "Installing symlink $(@F) into $(INSTALL_PREFIX)/lib/"
#	@$(SYMLINK) $(<F) $(@F)
#	@$(MV) $(@F) $(INSTALL_PREFIX)/lib/
#endif



# --- Clean rules ---
cleanmost: check-config
	- $(FIND) $(BASE_OBJ_DIR) -name "*.o" | $(XARGS) $(RM_F) 
	- $(FIND) $(BASE_LIB_DIR) -name "*.a" | $(XARGS) $(RM_F) 
	- $(FIND) $(BASE_LIB_DIR) -name "*.so" | $(XARGS) $(RM_F) 
	- $(RM_F) $(AR_OBJ_LIST_FILE)
#	- $(RM_F) $(AR_ARG_LIST_FILE)
	- $(RM_F) $(INCLUDE_LOCAL)/*.h

distclean: check-config cleanmost cleanmk
	- $(RM_RF) $(CONFIG_DIR)
	- $(RM_RF) $(OBJ_DIR)
	- $(RM_RF) $(LIB_DIR)
	- $(RM_RF) config.log
	- $(RM_RF) aclocal.m4
	- $(RM_RF) autom4te.cache
	- $(RM_RF) config.status
	- $(RM_RF) config.sys_type

cleanmk: check-config
	- $(FIND) $(SRC_DIR) -name "$(FRAGMENT_MK)" | $(XARGS) $(RM_F) 

cleanleaves: check-config
	- $(FIND) $(SRC_DIR) -name "*.[osx]" | $(XARGS) $(RM_F)



# --- Send thanks to FLAME group ---
send-thanks: check-config
	@echo $(THANKS_MSG_BODY) | $(MAIL) -s $(THANKS_MSG_SUBJECT) $(THANKS_MSG_EMAIL)

