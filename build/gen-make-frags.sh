#!/usr/bin/env bash

#
# gen-make-frags.sh
#
# A simple wrapper for the makefile fragment generator script.
#

# Generate makefile fragments recursively, level 1 verbosity.
./build/gen-make-frag.sh -rhv0 build/templates/fragment.mk src

