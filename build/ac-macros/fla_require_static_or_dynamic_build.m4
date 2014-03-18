AC_DEFUN([FLA_REQUIRE_STATIC_OR_DYNAMIC_BUILD],
[
	AC_REQUIRE([FLA_CHECK_ENABLE_STATIC_BUILD])
	AC_REQUIRE([FLA_CHECK_ENABLE_DYNAMIC_BUILD])

	dnl Make sure the user requested the static or the dynamic build (or both).

	dnl Scenario 3
	if test "$fla_enable_static_build" = "no" ; then
		if test "$fla_enable_dynamic_build" = "no" ; then
			AC_MSG_ERROR([Configuring libflame to disable both static and dynamic builds is not allowed. Please enable a static build or a dynamic build (or both) and then re-run configure.],[1])
		fi
	fi
])
