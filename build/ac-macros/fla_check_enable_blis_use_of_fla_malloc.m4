AC_DEFUN([FLA_CHECK_ENABLE_BLIS1_USE_OF_FLA_MALLOC],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested enabling code that defines bl1_malloc() in terms of FLA_malloc()])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([blis-use-of-fla-malloc],
	              AS_HELP_STRING([--enable-blis-use-of-fla-malloc],[Enable code that defines bl1_malloc() in terms of FLA_malloc(). One benefit of this is that BLIS memory allocations can be tracked, along with other libflame memory allocations, if the memory leak counter is enabled. A second benefit is that BLIS memory allocations can be aligned to boundaries if libflame memory alignment is enabled. Note this option may only be set at configure-time. (Enabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_blis_use_of_fla_malloc=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_blis_use_of_fla_malloc=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_BLIS1_USE_OF_FLA_MALLOC!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to enable the option.
		fla_enable_blis_use_of_fla_malloc=yes
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_blis_use_of_fla_malloc" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_BLIS1_USE_OF_FLA_MALLOC,1,
		          [Determines whether to define bl1_malloc() in terms of FLA_malloc().])
		
	elif test "$fla_enable_blis_use_of_fla_malloc" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_BLIS1_USE_OF_FLA_MALLOC!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_blis_use_of_fla_malloc)

])
