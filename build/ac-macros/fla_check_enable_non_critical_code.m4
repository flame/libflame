AC_DEFUN([FLA_CHECK_ENABLE_NON_CRITICAL_CODE],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested enabling non-critical code])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([non-critical-code],
	              AS_HELP_STRING([--enable-non-critical-code],[Enable code that provides non-critical functionality. This code has been identified as unnecessary when total library size is of concern. (Enabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_non_critical_code=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_non_critical_code=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_NON_CRITICAL_CODE!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to enable the option.
		fla_enable_non_critical_code=yes
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_non_critical_code" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_NON_CRITICAL_CODE,1,
		          [Determines whether to enable various segments of code identified as providing non-critical functionality.])
		
	elif test "$fla_enable_non_critical_code" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else

		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_NON_CRITICAL_CODE!]])

	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_non_critical_code)

])
