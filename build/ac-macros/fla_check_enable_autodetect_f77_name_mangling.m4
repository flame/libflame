AC_DEFUN([FLA_CHECK_ENABLE_AUTODETECT_F77_NAME_MANGLING],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested auto-detection of Fortran name-mangling])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([autodetect-f77-name-mangling],
	              AC_HELP_STRING([--enable-autodetect-f77-name-mangling],[Enable automatic detection of the (name-mangling) necessary to invoke Fortran routines from C, and C-compiled routines from Fortran. Disabling this option causes a pre-defined default to be used, which may not work in some environments. You may need to disable this option, along with autodetection of Fortran linker flags, if the environment's Fortran compiler is missing or broken. (Enabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then

			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_autodetect_f77_name_mangling=no

		elif test "$enableval" = "yes" ; then

			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_autodetect_f77_name_mangling=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_AUTODETECT_F77_NAME_MANGLING!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to enable the option.
		fla_enable_autodetect_f77_name_mangling=yes
	]
	)

	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_autodetect_f77_name_mangling" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
	elif test "$fla_enable_autodetect_f77_name_mangling" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])

	else

		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_AUTODETECT_F77_NAME_MANGLING!]])

	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_autodetect_f77_name_mangling)

])

