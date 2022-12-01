AC_DEFUN([FLA_CHECK_ENABLE_LAPACK2FLAME],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested the lapack2flame compatibility layer])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([lapack2flame],
	              AS_HELP_STRING([--enable-lapack2flame],[Compile and build into libflame a compatibility layer that maps LAPACK invocations to their corresponding FLAME/C implementations. Note that erroneous input parameters are reported according to libflame conventions, NOT LAPACK conventions. That is, if libflame error checking is disabled, no error checking is performed, and if erroneous input parameters are detected, the library aborts. Also, if this option is enabled, then external-lapack-for-subproblems and lapack2flash MUST be disabled. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_lapack2flame=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_lapack2flame=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_LAPACK2FLAME!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_lapack2flame=no
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_lapack2flame" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_LAPACK2FLAME,1,
		          [Determines whether the LAPACK compatibility layer is included in libflame.])
		
	elif test "$fla_enable_lapack2flame" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_LAPACK2FLAME!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_lapack2flame)

])
