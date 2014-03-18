AC_DEFUN([FLA_CHECK_ENABLE_EXTERNAL_LAPACK_INTERFACES],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested enabling interfaces to external LAPACK routines])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([external-lapack-interfaces],
	              AC_HELP_STRING([--enable-external-lapack-interfaces],[Enable code that allows the user to interface with an external LAPACK implementation via object-based FLAME-like functions. Note that if this option is enabled, an LAPACK library will be required at link-time. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_external_lapack_interfaces=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_external_lapack_interfaces=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_EXTERNAL_LAPACK_INTERFACES!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_external_lapack_interfaces=no
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_external_lapack_interfaces" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES,1,
		          [Determines whether to enable interfaces to external LAPACK routines.])
		
	elif test "$fla_enable_external_lapack_interfaces" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_EXTERNAL_LAPACK_INTERFACES!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_external_lapack_interfaces)

])
