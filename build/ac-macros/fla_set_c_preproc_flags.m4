AC_DEFUN([FLA_SET_C_PREPROC_FLAGS],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
	
	dnl Set C preprocessor flags assuming we found...
	case ${CC_VENDOR} in
		dnl GNU gcc.
		gcc)
			dnl Define the _GNU_SOURCE macro.
			AC_DEFINE(_GNU_SOURCE,1,
			          [Enables ANSI C, POSIX.1, POSIX.2, BSD, SVID, X/Open, and GNU extensions to the C language.])
		;;
		dnl Intel cc.
		icc)
		;;
		dnl PathScale pathcc.
		pathcc)
		;;
		dnl PGI pgcc.
		pgcc)
		;;
		dnl NEC sxcc.
		sxcc)
		;;
		dnl IBM xlc.
		*xlc*)
		;;
		dnl ambiguous cc.
		cc)
		;;
		dnl for all other C compilers.
		*)
		;;
	esac
])
