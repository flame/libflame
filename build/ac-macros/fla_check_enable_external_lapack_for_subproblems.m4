AC_DEFUN([FLA_CHECK_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested enabling external LAPACK for small subproblems])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([external-lapack-for-subproblems],
	              AC_HELP_STRING([--enable-external-lapack-for-subproblems],[Enable code that causes most of the computationally-intensive functions within libflame to compute their smallest subproblems by invoking a corresponding (usually unblocked) LAPACK routine. Note that if this option is enabled, lapack2flame MUST be disabled. Also, if this option is enabled, then external-lapack-interfaces MUST also be enabled. Enabling this option is useful when a libflame user wishes to leverage an optimized external implementation of LAPACK to speed up the b-by-b subproblems that arise within libflame's blocked algorithms and algorithms-by-blocks. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_external_lapack_for_subproblems=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_external_lapack_for_subproblems=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_external_lapack_for_subproblems=no
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_external_lapack_for_subproblems" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS,1,
		          [Determines whether to enable external LAPACK for small subproblems.])
		
	elif test "$fla_enable_external_lapack_for_subproblems" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_external_lapack_for_subproblems)

])
