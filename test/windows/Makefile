

#
# --- Include environment- and build-specific definitions ----------------------
#

CONFIG_MK = ..\..\windows\config\config.mk
DEFS_MK   = ..\..\windows\build\defs.mk

# Include build definitions
!if exist ( $(CONFIG_MK) )
!include $(CONFIG_MK)
!else
!error nmake: $(CONFIG_MK) does not exist! Run configure.cmd first.
!endif

# Include build definitions
!if exist ( $(DEFS_MK) )
!include $(DEFS_MK)
!else
!error nmake: $(DEFS_MK) does not exist! Your libflame distribution may be incomplete.
!endif



#
# --- General build definitions ------------------------------------------------
#

# Path to test suite source code.
SRC_PATH         = ..\src

# Path to directory where test suite object files will be stored.
OBJ_PATH         = obj

# Path to libflame and local headers.
LIBFLAME_HEADERS = $(INSTALL_PREFIX_INC)
LOCAL_HEADERS    = $(SRC_PATH)

# The installation paths for the static and dynamic libflame library builds.
LIB_LIBPATH      = $(INSTALL_PREFIX_LIB)
DLL_LIBPATH      = $(INSTALL_PREFIX_DLL)

# We have to reconstruct the C compiler flags with the new include paths.
CPPROCFLAGS      = /I$(LIBFLAME_HEADERS) /I$(LOCAL_HEADERS)
CFLAGS           = $(CMISCFLAGS) $(CLANGFLAGS) $(CPPROCFLAGS) $(CWARNFLAGS) \
                   $(CDBGFLAGS) $(COPTFLAGS) $(CRTIMEFLAGS) $(CMTHREADFLAGS)

# We use the C compiler as our linker.
LINKER           = $(CC)

# The linker flags. This might need to be modified depending on your build
# environment.
LDFLAGS          = /nologo /link \
                   /LIBPATH:"C:\Program Files\Microsoft SDKs\Windows\v6.0A\Lib" \
                   /LIBPATH:"C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\lib" \
                   /nodefaultlib:libcmt /nodefaultlib:libc /nodefaultlib:libmmt \
                   msvcrt.lib \
                   /LIBPATH:"C:\Program Files (x86)\Intel\Compiler\11.1\048\lib\ia32" \
                   /LIBPATH:"C:\Program Files (x86)\Intel\Compiler\11.1\048\mkl\ia32\lib" \
                   mkl_intel_c.lib \
                   mkl_sequential.lib \
                   mkl_core.lib

TEST_OBJS        = $(OBJ_PATH)\test_libflame.obj \
                   $(OBJ_PATH)\test_chol.obj \
                   $(OBJ_PATH)\test_lu_nopiv.obj \
                   $(OBJ_PATH)\test_lu_piv.obj \
                   $(OBJ_PATH)\test_lu_incpiv.obj \
                   $(OBJ_PATH)\test_qrut.obj \
                   $(OBJ_PATH)\test_qrutinc.obj \
                   $(OBJ_PATH)\test_lqut.obj \
                   $(OBJ_PATH)\test_uddateut.obj \
                   $(OBJ_PATH)\test_uddateutinc.obj \
                   $(OBJ_PATH)\test_apqudut.obj \
                   $(OBJ_PATH)\test_apqudutinc.obj \
                   $(OBJ_PATH)\test_hessut.obj \
                   $(OBJ_PATH)\test_tridiagut.obj \
                   $(OBJ_PATH)\test_bidiagut.obj \
                   $(OBJ_PATH)\test_trinv.obj \
                   $(OBJ_PATH)\test_spdinv.obj \
                   $(OBJ_PATH)\test_sylv.obj

# Binary path. We put it one directory up because that is where the input
# paremeter files reside.
TEST_BIN         = ..\test_libflame.exe

# The inference rule to compile source to object files.
{$(SRC_PATH)}.c{$(OBJ_PATH)}.obj:
	$(CC) $(CFLAGS) /c $< /Fo$@

# The rule to build the test suite driver by linking against a static
# build of libflame. This is the default.
with-lib: $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) /Fe$(TEST_BIN) $(LDFLAGS) $(LIB_LIBPATH)\$(LIBFLAME_LIB)

# The rule to build the test suite driver by linking against a dynamic
# build of libflame.
with-dll: $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) /Fe$(TEST_BIN) $(LDFLAGS) $(DLL_LIBPATH)\$(LIBFLAME_LIB)

clean:
	del /F /Q $(OBJ_PATH)\*.obj
	del /F /Q $(TEST_BIN)

