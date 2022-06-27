AC_DEFUN([FLA_CHECK_ENABLE_SCC],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested SCC extensions])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([scc],
	              AS_HELP_STRING([--enable-scc],[Enable code that takes advantage of the SCC multicore architecture. When using this option, enabling SuperMatrix is recommended, though not strictly required. Note that this option is experimental. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then

			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_scc=no

		elif test "$enableval" = "yes" ; then

			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_scc=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_SCC!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_scc=no
	]
	)

	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_scc" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_SCC,1,
		          [Determines whether SCC-specific blocks of code should be compiled.])
		
	elif test "$fla_enable_scc" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else

		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_SCC!]])

	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_scc)

])

