AC_DEFUN([FLA_REQUIRE_MULTITHREADING_ENABLED],
[
	AC_REQUIRE([FLA_CHECK_ENABLE_SUPERMATRIX])
	AC_REQUIRE([FLA_CHECK_ENABLE_MULTITHREADING])
	AC_REQUIRE([FLA_CHECK_ENABLE_HIP])

	dnl Make sure the user did request multithreading if HIP is
	dnl enabled.

	dnl Scenario 1: User wants HIP support but forgot to enable multithreading.
	if test "$fla_enable_hip" = "yes" ; then
		if test "$fla_enable_multithreading" = "no" ; then
			AC_MSG_ERROR([Configuring libflame to enable HIP support without also enabling multithreading is not allowed. HIP utilization requires that multithreading be enabled. Please adjust your configure options accordingly and then re-run configure.],[1])
		fi
	fi
])
