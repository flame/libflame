AC_DEFUN([FLA_CHECK_WITH_PYTHON],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested a specific python interpreter])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_WITH([python],
	            AC_HELP_STRING([--with-python=python],[Search for and use a python interpreter named <python>. If <python> is not found, then use the first interpreter found from the default search list.]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$withval" = "no" ; then
			
			dnl User provided --with-<option>=no or --without-<option>.
			fla_with_specific_python=without
			fla_requested_python=''

		elif test "$withval" = "yes" ; then
			
			dnl User provided --with-<option>=yes or --with-<option>.
			fla_with_specific_python=no
			fla_requested_python=''
		else
			
			dnl User provided argument value: --with-<option>=value.
			fla_with_specific_python=yes
			fla_requested_python=$withval
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_with_specific_python=no
		fla_requested_python=''
	]
	)
	
	dnl Now act according to whether a specific value was given.
	if test "$fla_with_specific_python" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes ($fla_requested_python)])
		
	elif test "$fla_with_specific_python" = "without" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
		dnl The user requested --with-python=no or --without-python. Scold him.
		AC_MSG_ERROR([[Detected --without-python. Cannot continue with python interpreter disabled!]])
		
	else
		dnl Output the result.
		AC_MSG_RESULT([no])
	fi
	
	dnl Check for PYTHON environment variable, which would override everything else.
	if test "$PYTHON" != "" ; then
		AC_MSG_NOTICE([[PYTHON environment variable is set to $PYTHON, which will override --with-python option and default search list for python interpreter.]])
	fi

])
