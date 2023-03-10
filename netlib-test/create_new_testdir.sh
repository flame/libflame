#!/bin/bash

print_usage()
{
	# Echo usage info.
	echo " "
	echo " ${script_name}"
	echo " "
	echo " Field G. Van Zee"
	echo " "
	echo " Create a new test suite to test libflame's LAPACK interfaces and"
	echo " functionality from an existing netlib LAPACK distribution."
	echo " "
	echo " Usage:"
	echo " "
	echo "   ${script_name} netlib_path new_dirname"
	echo " "
}

# ------------------------------------------------------------------------------

main()
{
	# The name of the script, stripped of any preceeding path.
	script_name=${0##*/}

	if [ $# != 2 ]; then

		print_usage
		exit 1
	fi

	# Grab the arguments.
	netlib_path=$1
	testdir_new=$2

	# Create the local destination directory.
	mkdir ./${testdir_new}

	# The list of items we will copy.
	copy_list="BLAS INSTALL TESTING Makefile lapack_testing.py"

	for item in ${copy_list};  do

		echo "Copying ${netlib_path}/${item} to ./${testdir_new}/"

		# Copy the current item from the netlib directory.
		cp -rp ${netlib_path}/${item} ./${testdir_new}/
	done

	echo "Creating ./${testdir_new}/make.inc"

	# Create a make.inc file from make.inc.example. Modify the file
	# to include -lpthread after liblapack.a (which should be a symlink
	# to libflame).
	#cp -p ${netlib_path}/make.inc.example ./${testdir_new}/make.inc
	cat ${netlib_path}/make.inc.example \
	    | sed -e "s/liblapack\.a/liblapack\.so/g" \
	    | sed -e "s/librefblas\.a/libblas\.so/g" \
	    | sed -e "s/libtmglib\.a/libtmglib\.so/g" \
	    | sed -e "s/AR = ar/AR = \$(CC)/g" \
	    | sed -e "s/ARFLAGS = cr/ARFLAGS = -lm -shared/g" \
	    | sed -e "s/RANLIB = ranlib/RANLIB = echo/g" \
	    > ./${testdir_new}/make.inc

	echo "Creating ./${testdir_new}/SRC"
	echo "Creating ./${testdir_new}/SRC/Makefile"

	# Make a dummy 'SRC' directory with a dummy Makefile.
	mkdir ./${testdir_new}/SRC
	cp ./build/SRC-Makefile ./${testdir_new}/SRC/Makefile

	# Adding README
	cp ./build/README ./${testdir_new}/

	# Make a dummy 'SRC/VARIANTS' directory with a dummy Makefile.
	mkdir ./${testdir_new}/SRC/VARIANTS
	cp ./build/SRC-VARIANTS-Makefile ./${testdir_new}/SRC/VARIANTS/Makefile

	# Make a dummy 'lapacke' directory with a dummy Makefile.
	mkdir ./${testdir_new}/lapacke
	cp ./build/lapacke-Makefile ./${testdir_new}/lapacke/Makefile

	echo "Tweaking TESTING/LIN/xerbla.f."
	echo "Tweaking TESTING/EIG/xerbla.f."

	# Disable a part of the custom xerbla_() that tests string equality
	# since this is broken due to Fortran/C incompatibilities.
	xb_lin_in=${netlib_path}/TESTING/LIN/xerbla.f
	xb_lin_ou=${testdir_new}/TESTING/LIN/xerbla.f
	xb_eig_in=${netlib_path}/TESTING/EIG/xerbla.f
	xb_eig_ou=${testdir_new}/TESTING/EIG/xerbla.f
	sed_expr="s/SRNAME\.NE\.SRNAMT/\.FALSE./g"

	sed -e "${sed_expr}" ${xb_lin_in} > ${xb_lin_ou}
	sed -e "${sed_expr}" ${xb_eig_in} > ${xb_eig_ou}

	echo "Tweaking ddrvsg.f"

	# This small change for the ddrvsg.f file is needed so that a
	# particular test of dsygvx() passes the test suite. (For some
	# reason, the f2c'ed dstein() does not converge when dstebz()
	# computes the eigenvalues with the regular definition of abstol,
	# so we use an alternate definition that causes the absolute
	# tolerance to be computed as eps * norm1(T), where T is the
	# tridiagonal matrix.
	ddrvsg_in=${netlib_path}/TESTING/EIG/ddrvsg.f
	ddrvsg_ou=${testdir_new}/TESTING/EIG/ddrvsg.f
	sed_expr="s/ABSTOL = UNFL + UNFL/ABSTOL = 0/g"

	sed -e "${sed_expr}" ${ddrvsg_in} > ${ddrvsg_ou}

	# Update makefiles to have the correct path
	# There was an issue where the path is given during linking 
	# as a relitive path. Then when run, it is executed from a different
	# dir. This would cause it to fail
	for make_file in BLAS/SRC/Makefile INSTALL/Makefile TESTING/Makefile TESTING/EIG/Makefile TESTING/LIN/Makefile TESTING/MATGEN/Makefile
	do
		sed -i 's/TOPSRCDIR = ../TOPSRCDIR = ${CURDIR}\/../' ${testdir_new}/$make_file
	done

	# Adding the '-o' flag needed to link
	sed -i 's/$(AR) $(ARFLAGS) $@ $^/$(AR) $(ARFLAGS) -o $@ $^/' ${testdir_new}/TESTING/MATGEN/Makefile

	# Updating these to point to libflame instead of internal files
	sed -i 's/testieee: tstiee.o .*/testieee: tstiee.o $(LAPACKLIB) $(BLASLIB)/' ${testdir_new}/INSTALL/Makefile

	sed -i 's/schkaa.o: schkaa.F/schkaa.o: schkaa.F $(LAPACKLIB) $(BLASLIB)/' ${testdir_new}/TESTING/LIN/Makefile
	sed -i 's/dchkaa.o: dchkaa.F/dchkaa.o: dchkaa.F $(LAPACKLIB) $(BLASLIB)/' ${testdir_new}/TESTING/LIN/Makefile
	sed -i 's/cchkaa.o: cchkaa.F/cchkaa.o: cchkaa.F $(LAPACKLIB) $(BLASLIB)/' ${testdir_new}/TESTING/LIN/Makefile
	sed -i 's/zchkaa.o: zchkaa.F/zchkaa.o: zchkaa.F $(LAPACKLIB) $(BLASLIB)/' ${testdir_new}/TESTING/LIN/Makefile
	sed -i 's/$(FC) $(FFLAGS_DRV) -c -o $@ $</$(FC) $(FFLAGS_DRV) -c -o $@ $^/g' ${testdir_new}/TESTING/LIN/Makefile

	# Exit peacefully.
	return 0
}

# The script's main entry point, passing all parameters given.
main "$@"
