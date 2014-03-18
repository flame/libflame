AC_DEFUN([FLA_CHECK_ENABLE_CBLAS_INTERFACES],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested using the CBLAS interfaces])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([cblas-interfaces],
	              AC_HELP_STRING([--enable-cblas-interfaces],[Enable code that interfaces libflame's external wrapper routines to the BLAS via the CBLAS rather than the traditional Fortran-77 API. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_cblas_interfaces=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_cblas_interfaces=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_CBLAS_INTERFACES!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_cblas_interfaces=no
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_cblas_interfaces" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_CBLAS_INTERFACES,1,
		          [Determines whether to enable CBLAS interfaces instead of Fortran-77 interfaces to the BLAS.])
		
	elif test "$fla_enable_cblas_interfaces" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_CBLAS_INTERFACES!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_cblas_interfaces)

])
