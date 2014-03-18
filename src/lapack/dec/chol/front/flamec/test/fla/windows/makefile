#
# Windows nmake makefile
#

FNAME        = Chol

SRC_PATH     = ..
OBJ_PATH     = .

INST_PATH    = c:\field\lib\libflame
INC_PATH     = $(INST_PATH)\include-x86-r2683
FLAME_LIB    = $(INST_PATH)\dll\libflame-x86-r2683.lib

CC           = cl.exe
FC           = ifort.exe
LINKER       = link.exe
CFLAGS       = /nologo /I$(INC_PATH) /O2
FFLAGS       = $(CFLAGS)
LDFLAGS      = /nologo 

TEST_OBJS    = $(OBJ_PATH)\test_$(FNAME).obj \
               $(OBJ_PATH)\time_$(FNAME).obj \
               $(OBJ_PATH)\REF_$(FNAME).obj
TEST_BIN     = test_$(FNAME).exe

{$(SRC_PATH)}.c{$(OBJ_PATH)}.obj:
	$(CC) $(CFLAGS) /c $< /Fo$@

{$(SRC_PATH)}.f{$(OBJ_PATH)}.obj:
	$(FC) $(FFLAGS) /c $< /Fo$@

all: $(TEST_OBJS)
	$(LINKER) $(LDFLAGS) $(TEST_OBJS) $(FLAME_LIB) 

clean:
	del /F /Q $(OBJ_PATH)\*.obj
	del /F /Q $(TEST_BIN)*.manifest
	del /F /Q $(TEST_BIN)

