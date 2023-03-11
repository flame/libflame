
LAPACK_PRES := $(shell if [ -e ../liblapack.so ]; then echo "yes"; else echo "no"; fi)
BLAS_PRES   := $(shell if [ -e ../libblas.so   ]; then echo "yes"; else echo "no"; fi)

all:
ifeq ($(LAPACK_PRES),yes)
	$(info Found 'liblapack.so'.)
else
	$(error Missing liblapack.so. Cannot continue. Please create a symlink named 'liblapack.so' and then re-run 'make'.)
endif
ifeq ($(BLAS_PRES),yes)
	$(info Found 'libblas.so'.)
else
	$(error Missing libblas.so. Cannot continue. Please create a symlink named 'libblas.so' and then re-run 'make'.)
endif
	@echo "Netlib LAPACK test suite prerequisites present."

clean:
	@echo "No objects to clean in SRC."
