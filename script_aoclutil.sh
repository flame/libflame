#!/bin/bash

#Set directory and git paths of libaoclutils library
LIBAOCLUTILS_DIR=libaoclutils
LIBAOCLUTILS_OBJ_DIR=libaoclutils_objdir
LIBAOCLUTILS_STATICLIB=libaoclutils.a
LIBAOCLUTILS_GIT_URL=https://github.com/amd/aocl-utils.git
LIBAOCLUTILS_GIT_TAG="4.1"
LIBAOCLUTILS_GIT_REPO=aocl-utils

for ARG in "$@"
do
   VAR=$(echo $ARG | cut -f1 -d=)
   DATA=$(echo $ARG | cut -f2 -d=)

   case "$VAR" in
         LIBAOCLUTILS_GIT_TAG)          if [[ ! -z "${DATA}" ]]
					then
       						LIBAOCLUTILS_GIT_TAG=${DATA}
					fi
					;;
         LIBAOCLUTILS_GIT_URL)          if [[ ! -z "${DATA}" ]]
					then
       						LIBAOCLUTILS_GIT_URL=${DATA}
					fi
					;;
         *)
   esac
done

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
	echo "Cloning libaoclutils from URL $LIBAOCLUTILS_GIT_URL and tag $LIBAOCLUTILS_GIT_TAG"
	git clone $LIBAOCLUTILS_GIT_URL -b $LIBAOCLUTILS_GIT_TAG

	#Check if clone succeeded
	if [ $? -eq 0 ]; 
	then	
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
	else
		echo "Cloning libaoclutils failed! Please check if clone path is correct."
		echo "Else provide path to prebuilt libaoclutils library using LIBAOCLUTILS_LIBRARY_PATH and LIBAOCLUTILS_INCLUDE_PATH."
	fi
fi
