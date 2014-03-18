#
# makefile for testing reference QR v. FLAMEC QR
#

LIB_PATH     = /home/field/flame/lib

BLAS         = $(LIB_PATH)/libgoto.so
#BLAS         = $(LIB_PATH)/libblas.a
LAPACK       = $(LIB_PATH)/liblapack.a

FLAME_BASE   = $(LIB_PATH)/libflame-base.a
FLAME_BLAS   = $(LIB_PATH)/libflame-blas.a
FLAME_LAPACK = $(LIB_PATH)/libflapack.a
LAPACK2FLAME = $(LIB_PATH)/liblapack2flame.a


FC           = g77
FFLAGS       = -O3 -inline-functions -fno-f2c

QR_OBJS      = qr_driver.o rand.o


all: ref fla

ref: $(QR_OBJS)
	g77 -o qr_driver_$@.x \
               $(QR_OBJS) \
               $(LAPACK) $(BLAS)

fla: $(QR_OBJS)
	g77 -o qr_driver_$@.x \
               $(QR_OBJS) \
               $(LAPACK2FLAME) \
               $(FLAME_LAPACK) $(FLAME_BLAS) $(FLAME_BASE) \
               $(LAPACK) $(BLAS)

clean:
	rm -f *.o *.x *~

