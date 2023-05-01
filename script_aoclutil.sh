#!/bin/bash

#Set directory and git paths of libaoclutils library
LIBAOCLUTILS_DIR=libaoclutils
LIBAOCLUTILS_OBJ_DIR=libaoclutils_objdir
LIBAOCLUTILS_STATICLIB=libaoclutils.a
LIBAOCLUTILS_GIT_URL=git@github.amd.com:AOCL/aocl-utils.git
LIBAOCLUTILS_GIT_TAG=amd-main
LIBAOCLUTILS_GIT_REPO=aocl-utils

#Check of libaoclutils object directory already exists and has object files
if [[ -d $LIBAOCLUTILS_OBJ_DIR  &&  ! -z `ls $LIBAOCLUTILS_OBJ_DIR/*.o` ]];
then
	echo "libaoclutils object files exist."
else
	#If object files not found, clone the aoclutils repository
	#and build it
	rm -rf $LIBAOCLUTILS_DIR
	mkdir -p $LIBAOCLUTILS_DIR
	mkdir -p $LIBAOCLUTILS_OBJ_DIR

	cd $LIBAOCLUTILS_DIR
	git clone $LIBAOCLUTILS_GIT_URL -b $LIBAOCLUTILS_GIT_TAG
	
	cd $LIBAOCLUTILS_GIT_REPO
	mkdir build_temp
	cd build_temp
	cmake ..
	make -j
	
	cp $LIBAOCLUTILS_STATICLIB ../../../$LIBAOCLUTILS_OBJ_DIR
	cd ../../../$LIBAOCLUTILS_OBJ_DIR

	#Extract the object files
	ar -x $LIBAOCLUTILS_STATICLIB
	cd ..
fi
