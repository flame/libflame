AC_DEFUN([FLA_CHECK_WITH_EXTRA_CFLAGS],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested using extra C compiler flags])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_WITH([extra-cflags],
	            AC_HELP_STRING([--with-extra-cflags=flagstring],[When compiling C code, use the flags in flagstring in addition to the flags that configure would normally use. This is useful when the user wants some extra flags passed to the compiler but does not want to manually set the CFLAGS environment variable and thus override all of the default compiler flags.]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$withval" = "no" ; then
			
			dnl User provided --with-<option>=no or --without-<option>.
			fla_with_extra_cflags=no
			fla_extra_cflags=''

		else
			
			dnl User provided argument value: --with-<option>=value.
			fla_with_extra_cflags=yes
			fla_extra_cflags=$withval
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_with_extra_cflags=no
		fla_extra_cflags=''
	]
	)
	
	dnl Now act according to whether a specific value was given.
	if test "$fla_with_extra_cflags" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else

		dnl Output the result.
		AC_MSG_RESULT([yes ($fla_extra_cflags)])
		
	fi
	
	dnl Check for CC environment variable, which would override everything else.
	if test "$EXTRA_CFLAGS" != "" ; then
		AC_MSG_NOTICE([[EXTRA_CFLAGS environment variable is set to $EXTRA_CFLAGS, which will override --with-extra-cflags option.]])
        fla_with_extra_cflags=yes
        fla_extra_cflags="$EXTRA_CFLAGS"
	fi

	dnl Declare EXTRA_CFLAGS as precious
	AC_ARG_VAR([EXTRA_CFLAGS],[extra C compiler flags to be used in addition to those determined automatically by configure])
	
	dnl Substitute the output variable.
	AC_SUBST(fla_with_extra_cflags)
	AC_SUBST(fla_extra_cflags)
])
