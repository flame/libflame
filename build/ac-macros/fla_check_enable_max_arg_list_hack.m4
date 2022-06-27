AC_DEFUN([FLA_CHECK_ENABLE_MAX_ARG_LIST_HACK],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested hack to circumvent small ARG_MAX])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([max-arg-list-hack],
	              AS_HELP_STRING([--enable-max-arg-list-hack],[Enable makefile code that archives object files from a flat object directory, thus decreasing the potential length of the argument list to ar. Use this option if you get 'Argument list too long' error messages when make tries to archive the library. (Enabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_max_arg_list_hack=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_max_arg_list_hack=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_BUILTIN_LAPACK_FUNCTIONS!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to enable the option.
		fla_enable_max_arg_list_hack=yes
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_max_arg_list_hack" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
	elif test "$fla_enable_max_arg_list_hack" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_BUILTIN_LAPACK_FUNCTIONS!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_max_arg_list_hack)

])
