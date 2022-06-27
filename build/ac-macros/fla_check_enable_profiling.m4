AC_DEFUN([FLA_CHECK_ENABLE_PROFILING],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether to enable compiler profiling symbols])

	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([profiling],
	              AS_HELP_STRING([--enable-profiling],[Use the appropriate profiling flag (usually -pg) when compiling C source code. (Disabled by default.)]),
	[
		
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_compiler_profiling=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_compiler_profiling=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_PROFILING!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_compiler_profiling=no
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_compiler_profiling" = "yes" ; then

		dnl Output the result.
		AC_MSG_RESULT([yes])
		
	elif test "$fla_enable_compiler_profiling" = "no" ; then

		dnl Output the result.
		AC_MSG_RESULT([no])

	else

		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_PROFILING!]])
	fi

	dnl Invoke a helper macro.
	FLA_SET_C_PROF_FLAGS($fla_enable_compiler_profiling)

	dnl Substitute the output variable.
	AC_SUBST(fla_enable_compiler_profiling)

])
