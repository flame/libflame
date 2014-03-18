AC_DEFUN([FLA_REQUIRE_AR],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
	
	dnl Declare AR as precious
	AC_ARG_VAR([AR],[library archiver])
	
	dnl If AR was not preset externally, then check for an archiver ourselves.
	if test "$AR" = "" ; then

		dnl If \$fla_requested_ar is set, then check for it first. This variable
		dnl was set in fla_check_with_ar. It is usually empty, but could be
		dnl non-empty if the user provided the --with-ar option to configure.
		if test "$fla_requested_ar" != "" ; then
			
			AC_CHECK_PROG([AR],$fla_requested_ar,$fla_requested_ar,[no])
		
			if test "$AR" = "no" ; then
				AC_MSG_WARN([Could not locate requested archiver ($fla_requested_ar)! Continuing search for default archiver ($fla_ar)],[1])
			fi
		fi
	
		dnl If the previous check for the requested archiver was unsuccessful in
		dnl setting AR, or if a specific archiver was not requested through
		dnl --with-ar to begin with, then check for the default archiver \$fla_ar.
		dnl This variable was set in fla_observe_cpu_type. Most of the time it is
		dnl simply set to "ar" but sometimes--for example, when cross-compiling
		dnl for the NEC SX-6--the default archiver is not named "ar".
		if test "$AR" = "no" || test "$AR" = "" ; then
			
			AC_CHECK_PROG([AR],$fla_ar,$fla_ar,[no])

			if test "$AR" = "no" ; then
				AC_MSG_ERROR([Could not locate $fla_ar! Cannot continue without a library archiver.],[1])
			fi
		fi
	
		dnl Determine the result of the macros and report an error if necessary. The
		dnl previous call to AC_CHECK_PROG also set AR to the value of \$fla_ar if 
		dnl the program by that name was found. If it was not found, then we can't
		dnl continue.
	fi
])
