AC_DEFUN([FLA_REQUIRE_FC],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])

	dnl Save the value of FFLAGS. This will come in useful later in determining
	dnl whether the user provided his own definition of FFLAGS.
	fla_userdef_fflags=$FFLAGS
	
	dnl Find a Fortran compiler.
	dnl If the FC environment variable is not already set, search for the
	dnl compiler defined by fla_requested_fc (which may be empty) and then
	dnl continue searching for the compilers in $fla_f_compiler_list, if
	dnl necessary.
	fla_f_compiler_list="gfortran f77 g77 xlf frt pgf77 cf77 fort77 fl32 af77 xlf90 f90 pgf90 pghpf epcf90 g95 xlf95 f95 fort ifort ifc efc pgfortran pgf95 lf95 ftn nagfor"
	AC_PROG_FC([$fla_requested_fc $fla_f_compiler_list])

	if test "$FC" = "" ; then
		AC_MSG_ERROR([Could not locate any of the following Fortran compilers: $FC $fla_requested_fc $fla_f_compiler_list. Cannot continue without a Fortran compiler.],[1])
	fi

	dnl Ascertain the compiler "vendor".
	dnl FC_VENDOR=$(echo "$FC" | egrep -o 'gfortran|f77|g77|xlf|frt|pgf77|cf77|fort77|fl32|af77|xlf90|f90|pgf90|pghpf|epcf90|g95|xlf95|f95|fort|ifort|ifc|efc|pgfortran|pgf95|lf95|ftn|nagfor' | { read first rest ; echo $first ; })
	dnl AC_MSG_NOTICE([[FC_VENDOR environment variable is set to ${FC_VENDOR}.]])
	
	dnl Substitute the user-defined FFLAGS into the autoconf output files.
	AC_SUBST(fla_userdef_fflags)
	dnl AC_SUBST(FC_VENDOR)
])
