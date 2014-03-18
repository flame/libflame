AC_DEFUN([FLA_SET_C_LANG_FLAGS],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
	
	AC_MSG_CHECKING([for (guessing) appropriate $CC language flags])

	dnl Set C compiler flags assuming we found...
	case $CC in
		dnl GNU gcc.
		gcc)
			fla_c_lang_flags='-std=c99'
		;;
		dnl Intel cc.
		icc)
			fla_c_lang_flags='-std=c99'
		;;
		dnl PathScale pathcc.
		pathcc)
			fla_c_lang_flags='-std=c99'
		;;
		dnl PGI pgcc.
		pgcc)
			fla_c_lang_flags='-c99'
		;;
		dnl NEC sxcc.
		sxcc)
			fla_c_lang_flags=''
		;;
		dnl IBM xlc.
		xlc)
			fla_c_lang_flags='-std=c99'
		;;
		dnl ambiguous cc.
		cc)
			fla_c_lang_flags='-std=c99'
		;;
		dnl for all other C compilers.
		*)
			fla_c_lang_flags=''
		;;
	esac

	dnl Output the result.
	AC_MSG_RESULT([$fla_c_lang_flags])
	
	dnl Substitute the language flags into the autoconf output files
	AC_SUBST(fla_c_lang_flags)

])
