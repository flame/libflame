AC_DEFUN([FLA_SET_C_SSE_FLAGS],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
		
	dnl Echo C SSE flags to user.
	AC_MSG_CHECKING([for (guessing) SSE flags for ${CC_VENDOR}])
		
	dnl Set the SSE compiler flags based on which C compiler we're going to use.
	case ${CC_VENDOR} in
		dnl Intel icc.
		icc)
			fla_c_sse_flags='-msse3'
		;;
		dnl GNU gcc.
		gcc)
			fla_c_sse_flags='-msse3'
		;;
		dnl Clang.
		clang)
			fla_c_sse_flags='-msse3'
		;;
		dnl PathScale pathcc
		pathcc)
			fla_c_sse_flags='unknown'
		;;
		dnl PGI pgcc
		pgcc)
			fla_c_sse_flags='-Mvect=sse'
		;;
		dnl NEC sxcc.
		sxcc)
			fla_c_sse_flags='notvalid'
		;;
		dnl IBM xlc
		*xlc*)
			fla_c_sse_flags='notvalid'
		;;
		dnl for all other C compilers.
		*)
			fla_c_sse_flags='unknown'
		;;
	esac
	
	dnl Output the result.
	AC_MSG_RESULT([$fla_c_sse_flags])
	
	dnl Check string in case C SSE flags are not valid.
	if test "$fla_c_sse_flags" = "notvalid" ; then
		
		dnl Tell the user we can't continue because he asked for SSE flags
		dnl for a compiler that (probably) does not support them.
		AC_MSG_ERROR([configure can't continue because the ${CC_VENDOR} compiler (probably) does not support SSE.])
	fi

	dnl Check string in case C SSE flags are unknown.
	if test "$fla_c_sse_flags" = "unknown" ; then
		
		dnl Tell the user we can't continue unless we know what flags
		dnl to pass to the C compiler to enable SSE support.
		AC_MSG_ERROR([configure doesn't know what flag to give ${CC_VENDOR} in order to enable SSE. Please submit a bug report to the FLAME developers.])
	fi

	dnl Output the C SSE flags variable.
	AC_SUBST(fla_c_sse_flags)
])
