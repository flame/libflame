AC_DEFUN([FLA_CHECK_ENABLE_LEGACY_LAPACK],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested the to use legacy lapack with lapack2flame or lapack2flash])

	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([legacy_lapack],
		AS_HELP_STRING([--enable-legacy-lapack],[Compile and build into libflame the legacy f2c files of LAPACK (Version 3.5.0). If this option is not selected, the latest integrated version of LAPACK using the fortran files will be built. This option only works if the --enable-lapack2flash or --enable-lapack2flame are given. Otherwise no LAPACK will be built. (Disabled by default.)]),
    [
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_legacy_lapack=no
		elif test "$enableval" = "yes" ; then
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_legacy_lapack=yes
		else
			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_LEGACY_LAPACK!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_legacy_lapack=no
	]
	)

	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_legacy_lapack" = "yes" ; then
		dnl Output the result.
		AC_MSG_RESULT([yes])

		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_LEGACY_LAPACK,1,
			[Determines whether to use the legacy f2c files or the newer fortran version of lapack with lapack2flame or lapack2flash.])
	elif test "$fla_enable_legacy_lapack" = "no" ; then
		dnl Output the result.
		AC_MSG_RESULT([no])
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_LEGACY_LAPACK!]])
	fi

	dnl Substitute the output variable.
	AC_SUBST(fla_enable_legacy_lapack)

])
