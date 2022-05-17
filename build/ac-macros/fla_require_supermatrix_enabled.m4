AC_DEFUN([FLA_REQUIRE_SUPERMATRIX_ENABLED],
[
	AC_REQUIRE([FLA_CHECK_ENABLE_SUPERMATRIX])
	AC_REQUIRE([FLA_CHECK_ENABLE_MULTITHREADING])
	AC_REQUIRE([FLA_CHECK_ENABLE_GPU])
	AC_REQUIRE([FLA_CHECK_ENABLE_HIP])

	dnl Make sure the user did not request an invalid combination of
	dnl parallelization options.

	dnl Scenario 1: User wants GPU support but forgot to enable SM.
	if test "$fla_enable_gpu" = "yes" ; then
		if test "$fla_enable_supermatrix" = "no" ; then
			AC_MSG_ERROR([Configuring libflame to enable GPU support without also enabling SuperMatrix is not allowed. GPU utilization requires that SuperMatrix be enabled. Please adjust your configure options accordingly and then re-run configure.],[1])
		fi
	fi

	dnl Scenario 2: User wants HIP support but forgot to enable SM.
	if test "$fla_enable_hip" = "yes" ; then
		if test "$fla_enable_supermatrix" = "no" ; then
			AC_MSG_ERROR([Configuring libflame to enable HIP support without also enabling SuperMatrix is not allowed. HIP utilization requires that SuperMatrix be enabled. Please adjust your configure options accordingly and then re-run configure.],[1])
		fi
	fi
])
