#!/bin/bash

#
# gen-make-frags.sh
#
# A simple wrapper for the makefile fragment generator script.
#

# Generate makefile fragments recursively, level 1 verbosity.
./build/gen-make-frag.sh -rhv1 build/templates/fragment.mk src

