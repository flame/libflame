AC_DEFUN([FLA_CHECK_ENABLE_PORTABLE_TIMER],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested a portable FLA_Clock() timer])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([portable-timer],
	              AC_HELP_STRING([--enable-portable-timer],[Define the FLA_Clock() timer function using clock_gettime(). If that function is not available, then getttimeofday() is used. If neither function is available, FLA_Clock() is will return a static value. (By default, a portable timer is used (if it exists).)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_portable_timer=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_portable_timer=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_PORTABLE_TIMER!]])
		fi
	],
	[
			dnl User did not specify whether to enable or disable the portable
			dnl implementation of FLA_Clock(). Default behavior is to enable
			dnl the portable timer.
			fla_enable_portable_timer=yes
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_portable_timer" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_PORTABLE_TIMER,1,
		          [Determines whether to define a portable FLA_Clock() in terms of clock_gettime() or gettimeofday() from time.h.])
		
		dnl Look for clock_gettime() and gettimeofday().
		AC_CHECK_FUNC([clock_gettime])
		AC_CHECK_FUNC([gettimeofday])

		fla_portable_timer_function=''

		dnl Define the right cpp macro depending on which function was
		dnl present in the environment.
		if   test "$ac_cv_func_clock_gettime" = "yes" ; then

			fla_portable_timer_function='clock_gettime()'
			AC_DEFINE(FLA_PORTABLE_TIMER_IS_CLOCK_GETTIME,1,
			          [Determines whether clock_gettime() was present on the system (via time.h).])

		elif test "$ac_cv_func_clock_gettime" = "no" ; then

			if   test "$ac_cv_func_gettimeofday" = "yes" ; then

				fla_portable_timer_function='gettimeofday()'
				AC_DEFINE(FLA_PORTABLE_TIMER_IS_GETTIMEOFDAY,1,
				          [Determines whether gettimeofday() was present on the system (via time.h).])

			elif test "$ac_cv_func_gettimeofday" = "no" ; then

				AC_MSG_ERROR([[Neither clock_gettime() nor gettimeofday() were found! FLA_Clock() will be broken!]])

				fla_portable_timer_function='none found!'
				AC_DEFINE(FLA_PORTABLE_TIMER_IS_UNKNOWN,1,
				          [Determines whether a timer was found at all.])
			fi
		fi
		
	elif test "$fla_enable_portable_timer" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
		dnl We'll need strchr(), which is used by detect_clocks() in
		dnl FLA_Clock.c. Verify that we have it.
		AC_CHECK_FUNC([strchr],[],
		[
			AC_MSG_ERROR([Failed to find a working version of strchr()! Try enabling the portable timer, which does not need strchr(), and rerun configure.])
		])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_PORTABLE_TIMER!]])
	fi
	
	dnl Substitute the output variables.
	AC_SUBST(fla_enable_portable_timer)
	AC_SUBST(fla_portable_timer_function)

])
