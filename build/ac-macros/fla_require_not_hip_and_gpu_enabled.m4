AC_DEFUN([FLA_REQUIRE_NOT_HIP_AND_GPU_ENABLED],
[
	AC_REQUIRE([FLA_CHECK_ENABLE_GPU])
	AC_REQUIRE([FLA_CHECK_ENABLE_HIP])

	dnl Make sure the user did not request an invalid combination of
	dnl acceleration options.

	dnl User enabled both GPU and HIP 
	if test "$fla_enable_gpu" = "yes" ; then
		if test "$fla_enable_hip" = "yes" ; then
			AC_MSG_ERROR([Configuring libflame to enable GPU and HIP support is not allowed. Please adjust your configure options accordingly and then re-run configure.],[1])
		fi
	fi
])
