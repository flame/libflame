#
# test directory makefile
#

FNAME        := Hevd_lv_components

SRC_PATH     := .
OBJ_PATH     := .

INST_PATH    := $(HOME)/flame
LIB_PATH     := $(INST_PATH)/lib
INC_PATH     := $(INST_PATH)/include
FLAME_LIB    := $(LIB_PATH)/libflame.a
LAPACK_LIB   := $(LIB_PATH)/liblapack.a
BLAS_LIB     := $(LIB_PATH)/libopenblas.a

CC           := gcc
FC           := gfortran
LINKER       := $(CC)
CFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O3 -Wall -Wno-comment -march=native -msse3 -malign-double -g
FFLAGS       := -I$(SRC_PATH) -I$(INC_PATH) -O3 -Wall
LDFLAGS      := -L/usr/lib/gcc/x86_64-linux-gnu/5 -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/u/field/intel/compilers_and_libraries_2017.3.191/linux/compiler/lib/intel64_lin -L/u/field/intel/compilers_and_libraries_2017.3.191/linux/mkl/lib/intel64_lin -L/usr/lib/gcc/x86_64-linux-gnu/5/../../.. -lgfortran -lm -lquadmath
LDFLAGS      += -lpthread -fopenmp

#TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))
TEST_OBJS    := test_$(FNAME).o \
                time_$(FNAME).o \
                REF_Hevd_lv.o \
                REF_Hevdd_lv.o \
                REF_Hevdr_lv.o \
                REF_Hevd_lv_components.o \
                REF_Hevdd_lv_components.o \
                REF_Hevdr_lv_components.o \
                FLA_Hevd_lv_var1_components.o \
                FLA_Hevd_lv_var2_components.o
TEST_BIN     := test_$(FNAME).x

test_$(FNAME): $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(FLAME_LIB) $(LAPACK_LIB) $(BLAS_LIB) $(LDFLAGS) -o $(TEST_BIN)

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o $(OBJ_PATH)/*.o *~ core *.x

