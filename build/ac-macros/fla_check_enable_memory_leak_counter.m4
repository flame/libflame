AC_DEFUN([FLA_CHECK_ENABLE_MEMORY_LEAK_COUNTER],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested enabling FLA_malloc()/FLA_free() memory leak counter])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([memory-leak-counter],
	              AC_HELP_STRING([--enable-memory-leak-counter],[Enable code that keeps track of the balance between calls to FLA_malloc() and FLA_free(). If enabled, the counter value is output to standard error upon calling FLA_Finalize(). Note that this option determines the default status, which may be changed at runtime. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_memory_leak_counter=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_memory_leak_counter=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_MEMORY_LEAK_COUNTER!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_memory_leak_counter=no
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_memory_leak_counter" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_MEMORY_LEAK_COUNTER,1,
		          [Determines whether to enable the FLA_malloc()/FLA_free() memory counter by default.])
		
	elif test "$fla_enable_memory_leak_counter" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_MEMORY_LEAK_COUNTER!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_memory_leak_counter)

])
