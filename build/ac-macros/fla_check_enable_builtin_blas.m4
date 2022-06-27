AC_DEFUN([FLA_CHECK_ENABLE_BUILTIN_BLAS],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested a built-in BLAS implementation])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([builtin-blas],
	              AS_HELP_STRING([--enable-builtin-blas],[Enable code that provides a built-in implementation of the BLAS. Note that some routines may not be optimized yet. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then

			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_builtin_blas=no

		elif test "$enableval" = "yes" ; then

			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_builtin_blas=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_BUILTIN_BLAS!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_builtin_blas=no
	]
	)

	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_builtin_blas" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_BUILTIN_BLAS,1,
		          [Determines whether a built-in implementation of the BLAS is compiled.])
		
	elif test "$fla_enable_builtin_blas" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else

		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_BUILTIN_BLAS!]])

	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_builtin_blas)

])

