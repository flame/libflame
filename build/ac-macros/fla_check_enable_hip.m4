AC_DEFUN([FLA_CHECK_ENABLE_HIP],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested HIP extensions])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([hip],
	              AC_HELP_STRING([--enable-hip],[Enable code that takes advantage of HIPs when performing certain computations. If enabled, SuperMatrix must also be enabled. Note that this option is experimental. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then

			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_hip=no

		elif test "$enableval" = "yes" ; then

			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_hip=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_HIP!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_hip=no
	]
	)

	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_hip" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_HIP,1,
		          [Determines whether HIP-specific blocks of code should be compiled.])
		
	elif test "$fla_enable_hip" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else

		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_HIP!]])

	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_hip)

])

