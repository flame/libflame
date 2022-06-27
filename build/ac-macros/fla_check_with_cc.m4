AC_DEFUN([FLA_CHECK_WITH_CC],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested a specific C compiler])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_WITH([cc],
	            AS_HELP_STRING([--with-cc=cc],[Search for and use a C compiler named <cc>. If <cc> is not found, then use the first compiler found from the default search list for the detected build architecture.]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$withval" = "no" ; then
			
			dnl User provided --with-<option>=no or --without-<option>.
			fla_with_specific_cc=without
			fla_requested_cc=''

		elif test "$withval" = "yes" ; then
			
			dnl User provided --with-<option>=yes or --with-<option>.
			fla_with_specific_cc=no
			fla_requested_cc=''
		else
			
			dnl User provided argument value: --with-<option>=value.
			fla_with_specific_cc=yes
			fla_requested_cc=$withval
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_with_specific_cc=no
		fla_requested_cc=''
	]
	)
	
	dnl Now act according to whether a specific value was given.
	if test "$fla_with_specific_cc" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes ($fla_requested_cc)])
		
	elif test "$fla_with_specific_cc" = "without" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
		dnl The user requested --with-cc=no or --without-cc. Scold him.
		AC_MSG_ERROR([[Detected --without-cc. Cannot continue with C compiler disabled!]])
		
	else
		dnl Output the result.
		AC_MSG_RESULT([no])
	fi
	
	dnl Check for CC environment variable, which would override everything else.
	if test "$CC" != "" ; then
		AC_MSG_NOTICE([[CC environment variable is set to $CC, which will override --with-cc option and default search list for C compiler.]])
	fi

])
