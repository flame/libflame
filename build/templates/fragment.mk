
#
# fragment.mk 
#
# This is an automatically-generated makefile fragment and will likely get
# overwritten or deleted if the user is not careful. Modify at your own risk.
#

# These two mmakefile variables need to be set in order for the recursive
# include process to work!
CURRENT_DIR_NAME := _mkfile_fragment_curr_dir_name_
CURRENT_SUB_DIRS := _mkfile_fragment_sub_dir_names_

# Source files local to this fragment
LOCAL_SRC_FILES  := _mkfile_fragment_local_src_files_

# Add the fragment's local source files to the _global_variable_ variable.
_mkfile_fragment_src_var_name_ += $(addprefix $(PARENT_SRC_PATH)/$(CURRENT_DIR_NAME)/, $(LOCAL_SRC_FILES))




# -----------------------------------------------------------------------------
# NOTE: The code below is generic and should remain in all fragment.mk files!
# -----------------------------------------------------------------------------

# Add the current fragment to the global list of fragments so the top-level
# Makefile knows which directories are participating in the build.
FRAGMENT_DIR_PATHS  += $(PARENT_SRC_PATH)/$(CURRENT_DIR_NAME)

# Recursively descend into other subfragments' local makefiles and include them.
ifneq ($(strip $(CURRENT_SUB_DIRS)),)
key1                := $(key1).x
key2                := $(key2).y
stack_$(key1)       := $(PARENT_PATH)
stack_$(key2)       := $(PARENT_SRC_PATH)
PARENT_PATH         := $(PARENT_PATH)/$(CURRENT_DIR_NAME)
PARENT_SRC_PATH     := $(PARENT_SRC_PATH)/$(CURRENT_DIR_NAME)
FRAGMENT_SUB_DIRS   := $(addprefix $(PARENT_PATH)/, $(CURRENT_SUB_DIRS))
-include  $(addsuffix /$(FRAGMENT_MK), $(FRAGMENT_SUB_DIRS))
PARENT_PATH         := $(stack_$(key1))
PARENT_SRC_PATH     := $(stack_$(key2))
key1                := $(basename $(key1))
key2                := $(basename $(key2))
endif
