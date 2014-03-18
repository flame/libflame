AC_DEFUN([FLA_REQUIRE_FIND],
[
	dnl Declare FIND as precious
	AC_ARG_VAR([FIND],[find utility])
	
	dnl Check that find is present.
	AC_CHECK_PROG([FIND],[find],[find],[not found!])
	
	dnl Determine the result of the macro and report an error if necessary.
	if test "$FIND" != "find" ; then
		AC_MSG_WARN([[Could not locate find! Note: make-cleaning makefile fragments requires find.]],[1])
	fi
])
