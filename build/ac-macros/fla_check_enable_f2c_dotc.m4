AC_DEFUN([FLA_CHECK_ENABLE_F2C_DOTC],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested enabling f2c based calling convention for complex dotc function])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([f2c-dotc],
	              AC_HELP_STRING([--enable-f2c-dotc],[Enable code that allows the user to invoke f2c based calling convention form complex dotc functio. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_f2c_dotc=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_f2c_dotc=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_F2C_DOTC!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_f2c_dotc=no
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_f2c_dotc" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_F2C_DOTC,1,
		          [Determines whether to enable f2c calling convention for complex dotc function.])
		
	elif test "$fla_enable_f2c_dotc" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_F2C_DOTC!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_f2c_dotc)

])
