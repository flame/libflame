AC_DEFUN([FLA_REQUIRE_CC_VENDOR],
[
	AC_REQUIRE([FLA_REQUIRE_CC])

	dnl Ascertain the compiler "vendor".
	CC_VENDOR=$(echo "$CC" | egrep -o 'icc|gcc|clang|emcc|pnacl|IBM|cc|craycc' | { read first rest ; echo $first ; })

	if test "${CC_VENDOR}" = "" ; then
		AC_MSG_ERROR([configure can't determine the compiler vendor for $CC. Please submit a bug report to the FLAME developers.])
	else
		AC_MSG_NOTICE([[setting the CC_VENDOR environment variable to ${CC_VENDOR}.]])
	fi

	dnl Substitute the user-defined CFLAGS into the autoconf output files.
	AC_SUBST(CC_VENDOR)
])
