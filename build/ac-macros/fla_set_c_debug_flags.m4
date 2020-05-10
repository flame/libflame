AC_DEFUN([FLA_SET_C_DEBUG_FLAGS],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
	
	AC_MSG_CHECKING([for (guessing) appropriate ${CC_VENDOR} debug flags])

	if test "$1" == "yes" ; then
		
		dnl Set C compiler flags assuming we found...
		case ${CC_VENDOR} in
			dnl GNU gcc.
			gcc)
				fla_c_debug_flags='-g'
			;;
			dnl Intel cc.
			icc)
				fla_c_debug_flags='-g'
			;;
			dnl PathScale pathcc.
			pathcc)
				fla_c_debug_flags='-g'
			;;
			dnl PGI pgcc.
			pgcc)
				fla_c_debug_flags='-g'
			;;
			dnl NEC sxcc.
			sxcc)
				fla_c_debug_flags='-g -C debug'
			;;
			dnl IBM xlc.
			*xlc*)
				fla_c_debug_flags='-g'
			;;
			dnl ambiguous cc.
			cc)
				fla_c_debug_flags='-g'
			;;
			dnl for all other C compilers.
			*)
				fla_c_debug_flags=''
			;;
		esac
	else

		dnl Set C compiler flags assuming we found...
		case ${CC_VENDOR} in
			dnl GNU gcc.
			gcc)
				fla_c_debug_flags='-g0'
			;;
			dnl Intel cc.
			icc)
				fla_c_debug_flags=''
			;;
			dnl PathScale pathcc.
			pathcc)
				fla_c_debug_flags=''
			;;
			dnl PGI pgcc.
			pgcc)
				fla_c_debug_flags=''
			;;
			dnl NEC sxcc.
			sxcc)
				fla_c_debug_flags='-Ng'
			;;
			dnl IBM xlc.
			xlc)
				fla_c_debug_flags=''
			;;
			dnl ambiguous cc.
			cc)
				fla_c_debug_flags=''
			;;
			dnl for all other C compilers.
			*)
				fla_c_debug_flags=''
			;;
		esac
	fi
	
	dnl Output the result.
	AC_MSG_RESULT([$fla_c_debug_flags])
	
	dnl Substitute the debug flags into the autoconf output files
	AC_SUBST(fla_c_debug_flags)

])
