#!/bin/bash

# This code allows us to invoke configure with the same relative
# path-to-the-top-level-directory component used to invoke this script.
script_name=${0##*/}
dist_path=${0%/run-conf/${script_name}}
configure_path=${dist_path}/configure

# The install prefix specifies the root directory of the installation.
# If you enable this option, generally speaking, libraries will be
# installed to $PREFIX_INST/lib and header files to $PREFIX_INST/include.
export PREFIX_INST="--prefix=$HOME/flame"

# If you want finer control over where libraries and headers are installed,
# change the paths passed to the --libdir and --includedir options to your
# desired installation paths. You can also comment these lines out if you
# are fine using the default installation path, which is based on the
# '/usr/local' prefix.
#export LIBDIR_INST="    --libdir=$HOME/flame/lib"
#export INCDIR_INST="--includedir=$HOME/flame/include"

# Override the default compiler search order and default compiler flags.
#export CC=clang
#export CC=icc

# Any options specified in CFLAGS will be preprended to the flags that
# libflame automatically determines.
export CFLAGS="-march=native"

# Run configure.
${configure_path} \
            ${PREFIX_INST} \
            ${LIBDIR_INST} \
            ${INCDIR_INST} \
            --disable-verbose-make-output \
            --enable-static-build \
            --enable-dynamic-build \
            --enable-autodetect-f77-ldflags \
            --enable-autodetect-f77-name-mangling \
            --disable-max-arg-list-hack \
            --enable-non-critical-code \
            --disable-builtin-blas \
            --disable-lapack2flame \
            --disable-lapack2flash \
            --disable-external-lapack-for-subproblems \
            --enable-external-lapack-interfaces \
            --disable-blas3-front-end-cntl-trees \
            --enable-multithreading=openmp \
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

