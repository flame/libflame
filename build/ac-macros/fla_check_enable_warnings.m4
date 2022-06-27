AC_DEFUN([FLA_CHECK_ENABLE_WARNINGS],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether to enable compiler warnings])

	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([warnings],
	              AS_HELP_STRING([--enable-warnings],[Use the appropriate flag(s) to request warnings when compiling C and Fortran source code. (Enabled by default.)]),
	[
		
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_compiler_warnings=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_compiler_warnings=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_WARNINGS!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_compiler_warnings=yes
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_compiler_warnings" = "yes" ; then

		dnl Output the result.
		AC_MSG_RESULT([yes])
		
	elif test "$fla_enable_compiler_warnings" = "no" ; then

		dnl Output the result.
		AC_MSG_RESULT([no])

	else

		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_WARNINGS!]])
	fi

	dnl Invoke helper macros.
	FLA_SET_C_WARNING_FLAGS($fla_enable_compiler_warnings)
		
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_compiler_warnings)

])
