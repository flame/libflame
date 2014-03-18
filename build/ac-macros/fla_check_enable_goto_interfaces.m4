AC_DEFUN([FLA_CHECK_ENABLE_GOTO_INTERFACES],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested enabling interfaces to internal/low-level libgoto functionality])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([goto-interfaces],
	              AC_HELP_STRING([--enable-goto-interfaces],[Enable code that interfaces with internal/low-level functionality within GotoBLAS, such as those symbols that may be queried for architecture-dependent blocksize values. When this option is disabled, reasonable static values are used instead. Note that in order to use libflame with a BLAS library other than GotoBLAS, the user must disable this option. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_goto_interfaces=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_goto_interfaces=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_GOTO_INTERFACES!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_goto_interfaces=no
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_goto_interfaces" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_GOTO_INTERFACES,1,
		          [Determines whether to enable interfaces to internal/low-level libgoto symbols.])
		
	elif test "$fla_enable_goto_interfaces" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_GOTO_INTERFACES!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_goto_interfaces)

])
