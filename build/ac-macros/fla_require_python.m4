AC_DEFUN([FLA_REQUIRE_PYTHON],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
	
	dnl Declare PYTHON as precious
	AC_ARG_VAR([PYTHON],[python interpreter])
	
	dnl If PYTHON was not preset externally, then check for an interpreter.
	if test "$PYTHON" = "" ; then

		dnl If \$fla_requested_python is set, then check for it first. This variable
		dnl was set in fla_check_with_python. It is usually empty, but could be
		dnl non-empty if the user provided the --with-python option to configure.
		if test "$fla_requested_python" != "" ; then
			
			AC_CHECK_PROG([PYTHON],$fla_requested_python,$fla_requested_python,[no])
		
			if test "$PYTHON" = "no" ; then
				AC_MSG_WARN([Could not locate requested python interpreter ($fla_requested_python)! Continuing search for default python interpreter ($fla_python)],[1])
			fi
		fi
	
		dnl If the previous check for the requested interpreter was unsuccessful in
		dnl setting PYTHON, or if a specific interpreter was not requested through
		dnl --with-python to begin with, then walk through the default interpreter
		dnl list in \$fla_python_list.
		if test "$PYTHON" = "no" || test "$PYTHON" = "" ; then
			
			for pyt in $fla_python_list ; do

				AC_CHECK_PROG([PYTHON],$pyt,$pyt,[no])

				if test "$PYTHON" = "no" ; then
					dnl python interpreter $pyt not found. Continue searching
					dnl the list.
					continue;
				else
					break;
				fi
			done

			if test "$PYTHON" = "no" ; then
				AC_MSG_ERROR([Could not locate any of the following python interpreters: $fla_requested_python $fla_python_list. Cannot continue without a python interpreter.],[1])
			fi
		fi
	fi
])
