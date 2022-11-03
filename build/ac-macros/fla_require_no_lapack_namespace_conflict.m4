AC_DEFUN([FLA_REQUIRE_NO_LAPACK_NAMESPACE_CONFLICT],
[
	AC_REQUIRE([FLA_CHECK_ENABLE_LAPACK2FLAME])
	AC_REQUIRE([FLA_CHECK_ENABLE_LAPACK2FLASH])
	AC_REQUIRE([FLA_CHECK_ENABLE_SUPERMATRIX])
	AC_REQUIRE([FLA_CHECK_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS])
	AC_REQUIRE([FLA_CHECK_ENABLE_EXTERNAL_LAPACK_INTERFACES])

	dnl Make sure the user did not request an invalid/dangerous combination of
	dnl LAPACK compatibility/interfacing options.

	dnl Scenarios 1 and 2
	if test "$fla_enable_lapack2flame" = "yes" ; then
		if test "$fla_enable_external_lapack_for_subproblems" = "yes" ; then
			AC_MSG_ERROR([Configuring libflame to enable both lapack2flame and external-lapack-for-subproblems is not allowed. lapack2flame requires that external-lapack-for-subproblems be disabled. Please adjust your configure options accordingly and then re-run configure.],[1])
		fi
	fi

	dnl Scenario 3, 4, and 5
	if test "$fla_enable_lapack2flash" = "yes" ; then
		if test "$fla_enable_lapack2flame" = "yes" ; then
			AC_MSG_ERROR([Configuring libflame to enable both lapack2flame and lapack2flash is not allowed. Please adjust your configure options accordingly and then re-run configure.],[1])
		fi
		if test "$fla_enable_external_lapack_for_subproblems" = "yes" ; then
			AC_MSG_ERROR([Configuring libflame to enable both lapack2flash and external-lapack-for-subproblems is not allowed. lapack2flash requires that external-lapack-for-subproblems be disabled. Please adjust your configure options accordingly and then re-run configure.],[1])
		fi
		if test "$fla_enable_supermatrix" = "no" ; then
			AC_MSG_ERROR([Configuring libflame to enable lapack2flash without also enabling SuperMatrix is not allowed. lapack2flash utilization requires that SuperMatrix be enabled. Please adjust your configure options accordingly and then re-run configure.],[1])
		fi
	fi

	dnl Scenario 6
	if test "$fla_enable_lapack2flame" = "no" ; then
		if test "$fla_enable_external_lapack_for_subproblems" = "yes" ; then
			if test "$fla_enable_external_lapack_interfaces" = "no" ; then
				AC_MSG_ERROR([Configuring libflame to enable external-lapack-for-subproblems without external-lapack-interfaces is not allowed. external-lapack-for-subproblems requires that external-lapack-interfaces be enabled. Please adjust your configure options accordingly and then re-run configure.],[1])
			fi
		fi
	fi
])
