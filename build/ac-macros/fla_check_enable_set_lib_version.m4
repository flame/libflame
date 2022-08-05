AC_DEFUN([FLA_CHECK_ENABLE_SET_LIB_VERSION],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested setting libflame version])
	
	dnl Determine whether the user gave the --set-lib-version option 
	AC_ARG_ENABLE([set-lib-version],
	              AC_HELP_STRING([--enable-set-lib-version],[Set libFLAME library version]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" != "" ; then
			
			dnl User provided a valid argument (hopefully). Enable.
		        dnl Define the preprocessor macro to enable the option.
		        AC_DEFINE(FLA_LIBFLAME_VERSION,"$enableval",
		                  [Set libFLAME version.])

		        dnl Output the result.
		        AC_MSG_RESULT([yes])

	        else

		       dnl Output the result.
		       AC_MSG_RESULT([no])

		fi
		
	],
	[

	]
	)

])

