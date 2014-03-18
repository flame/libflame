AC_DEFUN([FLA_REQUIRE_XARGS],
[
	dnl Declare XARGS as precious
	AC_ARG_VAR([XARGS],[xargs utility])
	
	dnl Check that xargs is present.
	AC_CHECK_PROG([XARGS],[xargs],[xargs],[not found!])
	
	dnl Determine the result of the macro and report an error if necessary.
	if test "$XARGS" != "xargs" ; then
		AC_MSG_WARN([[Could not locate xargs! Note: make-cleaning makefile fragments requires xargs.]],[1])
	fi
])
