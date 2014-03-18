AC_DEFUN([FLA_REQUIRE_CC],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])

	dnl Save the value of CFLAGS. This will come in useful later in determining
	dnl whether the user provided his own definition of CFLAGS.
	fla_userdef_cflags=$CFLAGS
	
	dnl Find a C compiler. 
	dnl If the CC environment variable is not already set, search for the
	dnl compiler defined by fla_requested_cc (which may be empty) and then
	dnl continue searching for the compilers in $fla_c_compiler_list, if
	dnl necessary. Also, if the C compiler found is not in ANSI mode, then
	dnl try to add an option to make it so. If the GNU gcc was found, then
	dnl GCC shell variable is set to `yes'. 
	AC_PROG_CC([$fla_requested_cc $fla_c_compiler_list])

	if test "$CC" = "" ; then
		AC_MSG_ERROR([Could not locate any of the following C compilers: $CC $fla_requested_cc $fla_c_compiler_list. Cannot continue without a C compiler.],[1])
	fi
	
	dnl Substitute the user-defined CFLAGS into the autoconf output files.
	AC_SUBST(fla_userdef_cflags)
])
