#!/bin/bash

# The install prefix specifies the root directory of the installation.
export FLAME_INSTALL_PREFIX=$HOME/flame

# Override the default compiler search order and default compiler flags.
export CC=gcc
#export CC=icc

# Any options specified in CFLAGS will be preprended to the flags that
# libflame automatically determines.
export CFLAGS="-march=native"

# Run configure.
./configure --prefix=${FLAME_INSTALL_PREFIX} \
            --disable-verbose-make-output \
            --enable-static-build \
            --enable-dynamic-build \
            --enable-autodetect-f77-ldflags \
            --enable-autodetect-f77-name-mangling \
            --disable-max-arg-list-hack \
            --enable-non-critical-code \
            --disable-builtin-blas \
            --enable-lapack2flame \
            --disable-external-lapack-for-subproblems \
            --enable-external-lapack-interfaces \
            --disable-blas3-front-end-cntl-trees \
            --enable-multithreading=pthreads \
            --enable-supermatrix \
            --disable-gpu \
            --disable-vector-intrinsics \
            --enable-memory-alignment=16 \
            --enable-ldim-alignment \
            --enable-optimizations \
            --enable-warnings \
            --enable-debug \
            --disable-profiling \
            --enable-internal-error-checking=full \
            --disable-memory-leak-counter \
            --enable-blis-use-of-fla-malloc \
            --disable-goto-interfaces \
            --disable-cblas-interfaces \
            --enable-portable-timer \
            --disable-windows-build \
            --disable-scc \
            --disable-tidsp

