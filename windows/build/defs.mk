

#
# --- General build system options --------------------------------------------
#

# Uncomment this for verbose output from nmake.
# VERBOSE = 1

# Assign this varible to be the full path to the directory to which you would
# like the libflame build products to be installed upon running "nmake install".
# The nmake install target will create the install directory and all requisite
# subdirectories if they do not already exist (in which case the user must have
# permission to create these directories).
INSTALL_PREFIX = c:\field\lib


#
# --- Important build system filenames ----------------------------------------
#

# DLL link arguments. The contents of this file should be customized when
# building a dynamically-linked library. The lines of the file should contain
# linker options, library names, and library paths. Note that the library
# paths must be declared in the following form:
#
#   /link /LIBPATH:<path1>
#   /link /LIBPATH:<path2>
#   /link /LIBPATH:<path3>
#
# where <path1>, <path2>, and <path3> are library paths to add to the list
# of paths to search when the linker attempts to locate other libraries
# listed in the file.
LINKARGS_FILENAME = linkargs.txt
LINKARGS_FILEPATH = $(PWD)\$(LINKARGS_FILENAME)

# Various log file names that capture standard output when VERBOSE is undefined.
CC_LOG_FILE   = nmake-cc.log
FC_LOG_FILE   = nmake-fc.log
COPY_LOG_FILE = nmake-copy.log


#
# --- General name and directory definitions -----------------------------------
#

# The relative and absolute locations of the top-level Windows build directory.
# This is the directory in which nmake is run (not the directory named "build").
TOP_BUILD_DIR_REL = .
TOP_BUILD_DIR_ABS = $(PWD)

# The revision string.
REV_STR           = r$(REVISION)

# The names of the libraries.
LIBFLAME_NAME_ONLY        = libflame
LIBFLAME                  = $(LIBFLAME_NAME_ONLY)-$(ARCH_STR)-$(REV_STR)

# Directories that reside within the top-level Windows directory.
CNF_DIRNAME       = config
INC_DIRNAME       = include
SRC_DIRNAME       = src
OBJ_DIRNAME       = obj
LIB_DIRNAME       = lib
DLL_DIRNAME       = dll

# Leaves of interest for Windows.
FLA_DIRNAME       = flamec
L2F_DIRNAME       = lapack2flamec

# Relative directory paths to each of the above subdirectories.
INC_DIRPATH       = $(TOP_BUILD_DIR_REL)\$(INC_DIRNAME)
SRC_DIRPATH       = $(TOP_BUILD_DIR_REL)\$(SRC_DIRNAME)
OBJ_DIRPATH       = $(TOP_BUILD_DIR_REL)\$(OBJ_DIRNAME)
LIB_DIRPATH       = $(TOP_BUILD_DIR_REL)\$(LIB_DIRNAME)
DLL_DIRPATH       = $(TOP_BUILD_DIR_REL)\$(DLL_DIRNAME)

# We only have header files for flamec leaves.
INC_FLA_DIRPATH   = $(INC_DIRPATH)\$(FLA_DIRNAME)

# We have source code for flamec and lapack2flamec leaves.
SRC_FLA_DIRPATH   = $(SRC_DIRPATH)\$(FLA_DIRNAME)
SRC_L2F_DIRPATH   = $(SRC_DIRPATH)\$(L2F_DIRNAME)

# And we have object file paths corresponding to those source leaves defined
# above.
OBJ_FLA_DIRPATH   = $(OBJ_DIRPATH)\$(FLA_DIRNAME)\$(ARCH_STR)\$(BUILD_STR)
OBJ_L2F_DIRPATH   = $(OBJ_DIRPATH)\$(L2F_DIRNAME)\$(ARCH_STR)\$(BUILD_STR)

# Separate directories into which we'll move object files when we create the
# static libraries.
LIB_LIBFLAME_DIRPATH = $(LIB_DIRPATH)\$(ARCH_STR)\$(BUILD_STR)

# Separate directories into which we'll move object files when we create the
# dynamic libraries.
DLL_LIBFLAME_DIRPATH = $(DLL_DIRPATH)\$(ARCH_STR)\$(BUILD_STR)

# The install subdirectories.
INSTALL_PREFIX_LIB = $(INSTALL_PREFIX)\libflame\lib
INSTALL_PREFIX_DLL = $(INSTALL_PREFIX)\libflame\dll
INSTALL_PREFIX_INC = $(INSTALL_PREFIX)\libflame\include-$(ARCH_STR)-$(REV_STR)

# Definitions for important header files used in the install-headers rule.
BUILD_DIRNAME      = build
FLA_CONFIG_H       = FLA_config.h
FLAME_H            = FLAME.h


#
# --- General shell definitions ------------------------------------------------
#

CD     = cd
DIR    = dir
COPY   = copy
DEL    = del /F /Q
MKDIR  = mkdir
RMDIR  = rd /S /Q
ECHO   = echo


#
# --- Helper scripts -----------------------------------------------------------
#

NMAKE_HELP = .\build\nmake-help.cmd



#
# --- Compiler-related definitions ---------------------------------------------
#

# --- C compiler definitions ---

!if "$(CCOMPILER_STR)"=="icl"

!if "$(BUILD_STR)"=="debug"
CDEBUG = /Zi
COPTIM = /Od
!elseif "$(BUILD_STR)"=="release"
CDEBUG =
COPTIM = /Ox
!endif

CC            = icl.exe
CMISCFLAGS    = /nologo
CLANGFLAGS    =
CPPROCFLAGS   = /I.\build /I$(INC_FLA_DIRPATH)
CWARNFLAGS    = /w
CDBGFLAGS     = $(CDEBUG)
COPTFLAGS     = $(COPTIM)
CRTIMEFLAGS   = /MT
CMTHREADFLAGS = /Qopenmp
CFLAGS        = $(CMISCFLAGS) $(CLANGFLAGS) $(CPPROCFLAGS) $(CWARNFLAGS) \
                $(CDBGFLAGS) $(COPTFLAGS) $(CRTIMEFLAGS) $(CMTHREADFLAGS)

!elseif "$(CCOMPILER_STR)"=="cl"

!if "$(BUILD_STR)"=="debug"
CDEBUG = /Zi
COPTIM = /Od
!elseif "$(BUILD_STR)"=="release"
CDEBUG =
COPTIM = /Ox
!endif

CC            = cl.exe
CMISCFLAGS    = /nologo
CLANGFLAGS    =
CPPROCFLAGS   = /I.\build /I$(INC_FLA_DIRPATH)
CWARNFLAGS    = /w
CDBGFLAGS     = $(CDEBUG)
COPTFLAGS     = $(COPTIM)
CRTIMEFLAGS   = /MT
CMTHREADFLAGS = /openmp
CFLAGS        = $(CMISCFLAGS) $(CLANGFLAGS) $(CPPROCFLAGS) $(CWARNFLAGS) \
                $(CDBGFLAGS) $(COPTFLAGS) $(CRTIMEFLAGS) $(CMTHREADFLAGS)

!endif



#
# --- Library-related definitions ----------------------------------------------
#

# --- Static library definitions ---

LIBFLAME_LIB          = $(LIBFLAME).lib

LIB                   = lib
LIB_OPTIONS           = /nologo
LIB_FLA_OUTPUT_ARG    = /out:$(LIBFLAME_LIB)
LIB_FLA_INPUT_ARGS    = *.obj

# --- Dynamic library definitions ---

LIBFLAME_DLL          = $(LIBFLAME).dll

GENDLL                = $(TOP_BUILD_DIR_ABS)\gendll.cmd
OBJ_LIST_FILE         = libflame-objects.txt

SYM_DEF_FILEPATH      = $(TOP_BUILD_DIR_ABS)\$(BUILD_DIRNAME)\libflame-symbols.def

