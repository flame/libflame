AC_DEFUN([FLA_REQUIRE_RANLIB],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
	
	dnl Declare RANLIB as precious
	AC_ARG_VAR([RANLIB],[library archive indexer])
	
	dnl If RANLIB was not preset externally, then check for an indexer ourselves.
	if test "$RANLIB" = "" ; then

		dnl If $fla_requested_ranlib is set, then check for it first. This variable
		dnl was set in fla_check_with_ranlib. It is usually empty, but could be
		dnl non-empty if the user provided the --with-ranlib option to configure.
		if test "$fla_requested_ranlib" != "" ; then
			
			AC_CHECK_PROG([AR],$fla_requested_ranlib,$fla_requested_ranlib,[no])
		
			if test "$RANLIB" = "no" ; then
				AC_MSG_WARN([Could not locate requested archive indexer ($fla_requested_ranlib)! Continuing search for default archive indexer ($fla_ranlib)],[1])
			fi
		fi
	
		dnl If the previous check for the requested indexer was unsuccessful in
		dnl setting RANLIB, or if a specific indexer was not requested through
		dnl --with-ranlib to begin with, then check for the default indexer
		dnl $fla_ranlib. This variable was set in fla_observe_cpu_type. Most of
		dnl the time it is simply set to "ranlib" but sometimes the default
		dnl indexer is not named "ranlib".
		if test "$RANLIB" = "no" || test "$RANLIB" = "" ; then
			
			AC_CHECK_PROG([RANLIB],$fla_ranlib,$fla_ranlib,[no])

			if test "$RANLIB" = "no" ; then
				AC_MSG_ERROR([Could not locate $fla_ranlib! Cannot continue without a library archive indexer.],[1])
			fi
		fi
	
		dnl Determine the result of the macros and report an error if necessary. The
		dnl previous call to AC_CHECK_PROG also set RANLIB to the value of \$fla_ranlib 
		dnl if the program by that name was found. If it was not found, then we can't
		dnl continue.
	fi
])
