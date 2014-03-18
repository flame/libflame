AC_DEFUN([FLA_CHECK_ENABLE_DEBUG],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether to enable compiler debugging symbols])

	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([debug],
	              AC_HELP_STRING([--enable-debug],[Use the appropriate debug flag (usually -g) when compiling C and Fortran source code. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
		
			dnl User provided --enable-<option>=no or --disable-<option>.	
			fla_enable_compiler_debug=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_compiler_debug=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_DEBUG!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_compiler_debug=no
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_compiler_debug" = "yes" ; then

		dnl Output the result.
		AC_MSG_RESULT([yes])
		
	elif test "$fla_enable_compiler_debug" = "no" ; then

		dnl Output the result.
		AC_MSG_RESULT([no])

	else

		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_DEBUG!]])
	fi

	dnl Set the appropriate debug flags
	FLA_SET_C_DEBUG_FLAGS($fla_enable_compiler_debug)
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_compiler_debug)

])
