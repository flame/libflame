AC_DEFUN([FLA_CHECK_WITH_AR],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested a specific library archiver])
	
	dnl Determine whether the user gave the --enable-<option> or
    dnl --disable-<option>. If so, then run the first snippet of code;
    dnl otherwise, run the second code block.
	AC_ARG_WITH([ar],
	            AC_HELP_STRING([--with-ar=ar],[ Search for and use a library archiver named <ar>. If <ar> is not found, then use the first library archiver found from the default search list for the detected build architecture. Note: the library archiver search list usually consists only of "ar".]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$withval" = "no" ; then
			
			dnl User provided --with-<option>=no or --without-<option>.
			fla_with_specific_ar=without
			fla_requested_ar=''

		elif test "$withval" = "yes" ; then
			
			dnl User provided --with-<option>=yes or --with-<option>.
			fla_with_specific_ar=no
			fla_requested_ar=''
		else
			
			dnl User provided argument value: --with-<option>=value.
			fla_with_specific_ar=yes
			fla_requested_ar=$withval
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
        dnl Default behavior is to disable the option.
		fla_with_specific_ar=no
		fla_requested_ar=''
	]
	)
	
	dnl Now act according to whether a specific value was given.
	if test "$fla_with_specific_ar" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes ($fla_requested_ar)])
		
	elif test "$fla_with_specific_ar" = "without" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
		dnl The user requested --with-ar=no or --without-ar. Scold him.
		AC_MSG_ERROR([[Detected --without-ar. Cannot continue with library archiver disabled!]])
		
	else
		dnl Output the result.
		AC_MSG_RESULT([no])
	fi
	
	dnl Check for AR environment variable, which would override everything else.
	if test "$AR" != "" ; then
		AC_MSG_NOTICE([[AR environment variable is set to $AR, which will override --with-ar option and default search list for library archiver.]])
	fi

])
