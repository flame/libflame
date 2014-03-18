AC_DEFUN([FLA_REQUIRE_GNU_BASH],
[
	dnl Output a message saying we're going to check for GNU bash.
	AC_MSG_CHECKING([[for GNU bash]])

	dnl Check the cache for the result. If it's not result is not yet cached,
	dnl execute the shell code to set the cache-id (ie: look for bash).
	AC_CACHE_VAL([_cv_gnu_bash_command],
	[
		dnl Initialize the cache-id to null.
		_cv_gnu_bash_command='';
	
		dnl Check that some version of bash is present.
		dnl If bash is not present, then output an error message.
		if ( sh -c "bash --version" 2> /dev/null | grep GNU 2>&1 > /dev/null ); then
			_cv_gnu_bash_command='bash';
		fi
	])
	
	dnl Now that we've checked for GNU bash, let's determine the result of the
	dnl macro and report and error if bash was not found.
	if test "$_cv_gnu_bash_command" = "bash" ; then
		AC_MSG_RESULT([[bash]])
	else
		AC_MSG_RESULT([[not found!]])
		AC_MSG_ERROR([[Could not locate GNU bash! Bailing out...]],[1])
	fi
])
