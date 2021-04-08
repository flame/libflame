
LAPACK_PRES := $(shell if [ -e ../liblapack.a ]; then echo "yes"; else echo "no"; fi)
BLAS_PRES   := $(shell if [ -e ../libblas.a   ]; then echo "yes"; else echo "no"; fi)

all:
ifeq ($(LAPACK_PRES),yes)
	$(info Found 'liblapack.a'.)
else
	$(error Missing liblapack.a. Cannot continue. Please create a symlink named 'liblapack.a' and then re-run 'make'.)
endif
ifeq ($(BLAS_PRES),yes)
	$(info Found 'libblas.a'.)
else
	$(error Missing libblas.a. Cannot continue. Please create a symlink named 'libblas.a' and then re-run 'make'.)
endif
	@echo "Netlib LAPACK test suite prerequisites present."

clean:
	@echo "No objects to clean in SRC."
