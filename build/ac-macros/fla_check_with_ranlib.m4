AC_DEFUN([FLA_CHECK_WITH_RANLIB],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested a specific library archive indexer])
	
	dnl Determine whether the user gave the --enable-<option> or
    dnl --disable-<option>. If so, then run the first snippet of code;
    dnl otherwise, run the second code block.
	AC_ARG_WITH([ranlib],
	            AC_HELP_STRING([--with-ranlib=ranlib],[ Search for and use a library archive indexer named <ranlib>. If <ranlib> is not found, then use the first library archive indexer found from the default search list for the detected build architecture. Note: the library archive indexer search list usually consists only of "ranlib".]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$withval" = "no" ; then
			
			dnl User provided --with-<option>=no or --without-<option>.
			fla_with_specific_ranlib=without
			fla_requested_ranlib=''

		elif test "$withval" = "yes" ; then
			
			dnl User provided --with-<option>=yes or --with-<option>.
			fla_with_specific_ranlib=no
			fla_requested_ranlib=''
		else
			
			dnl User provided argument value: --with-<option>=value.
			fla_with_specific_ranlib=yes
			fla_requested_ranlib=$withval
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
        dnl Default behavior is to disable the option.
		fla_with_specific_ranlib=no
		fla_requested_ranlib=''
	]
	)
	
	dnl Now act according to whether a specific value was given.
	if test "$fla_with_specific_ranlib" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes ($fla_requested_ranlib)])
		
	elif test "$fla_with_specific_ranlib" = "without" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
		dnl The user requested --with-ranlib=no or --without-ranlib. Scold him.
		AC_MSG_ERROR([[Detected --without-ranlib. Cannot continue with library archive indexer disabled!]])
		
	else
		dnl Output the result.
		AC_MSG_RESULT([no])
	fi
	
	dnl Check for RANLIB environment variable, which would override everything else.
	if test "$RANLIB" != "" ; then
		AC_MSG_NOTICE([[RANLIB environment variable is set to $RANLIB, which will override --with-ranlib option and default search list for library archive indexer.]])
	fi

])
