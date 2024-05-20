AC_DEFUN([FLA_CHECK_WITH_FC],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested a specific Fortran compiler])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_WITH([fc],
	            AS_HELP_STRING([--with-fc=fc],[Search for and use a Fortran compiler named <fc>. If <fc> is not found, then use the first compiler found from the default search list for the detected build architecture.]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$withval" = "no" ; then
			
			dnl User provided --with-<option>=no or --without-<option>.
			fla_with_specific_fc=without
			fla_requested_fc=''

		elif test "$withval" = "yes" ; then
			
			dnl User provided --with-<option>=yes or --with-<option>.
			fla_with_specific_fc=no
			fla_requested_fc=''
		else
			
			dnl User provided argument value: --with-<option>=value.
			fla_with_specific_fc=yes
			fla_requested_fc=$withval
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_with_specific_fc=no
		fla_requested_fc=''
	]
	)
	
	dnl Now act according to whether a specific value was given.
	if test "$fla_with_specific_fc" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes ($fla_requested_fc)])
		
	elif test "$fla_with_specific_fc" = "without" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
		dnl The user requested --with-fc=no or --without-fc. Scold him.
		AC_MSG_ERROR([[Detected --without-fc. Cannot continue with Fortran compiler disabled!]])
		
	else
		dnl Output the result.
		AC_MSG_RESULT([no])
	fi
	
	dnl Check for FC environment variable, which would override everything else.
	if test "$FC" != "" ; then
		AC_MSG_NOTICE([[FC environment variable is set to $FC, which will override --with-fc option and default search list for Fortran compiler.]])
	fi

])
