AC_DEFUN([FLA_CHECK_ENABLE_AUTODETECT_F77_LDFLAGS],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested auto-detection of Fortran linker flags])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([autodetect-f77-ldflags],
	              AC_HELP_STRING([--enable-autodetect-f77-ldflags],[Enable automatic detection of any linker flags that may be needed to link against Fortran code. These flags are useful to know about when, for example, linking libflame to a BLAS library that was compiled with the system's Fortran compiler. You may need to disable this option, along with autodetection of Fortran name-mangling, if the environment's Fortran compiler is missing or broken. (Enabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then

			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_autodetect_f77_ldflags=no

		elif test "$enableval" = "yes" ; then

			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_autodetect_f77_ldflags=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_AUTODETECT_F77_LDFLAGS!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to enable the option.
		fla_enable_autodetect_f77_ldflags=yes
	]
	)

	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_autodetect_f77_ldflags" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
	elif test "$fla_enable_autodetect_f77_ldflags" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])

	else

		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_AUTODETECT_F77_LDFLAGS!]])

	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_autodetect_f77_ldflags)

])

