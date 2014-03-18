AC_DEFUN([FLA_SET_C_WARNING_FLAGS],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
	
	AC_MSG_CHECKING([for (guessing) appropriate $CC warning flags])

	if test "$1" == "yes" ; then
		
		dnl Set C compiler flags assuming we found...
		case $CC in
			dnl GNU gcc.
			gcc)
				fla_c_warning_flags='-Wall -Wno-comment'
			;;
			dnl Intel cc.
			icc)
				fla_c_warning_flags='-Wall -wd869,981,1418,1419,1572'
			;;
			dnl PathScale pathcc.
			pathcc)
				fla_c_warning_flags='-Wall'
			;;
			dnl PGI pgcc.
			pgcc)
				fla_c_warning_flags='-Minform=warn'
			;;
			dnl NEC sxcc.
			sxcc)
				fla_c_warning_flags='-w all'
			;;
			dnl IBM xlc.
			xlc)
				fla_c_warning_flags='-qcpluscmt'
			;;
			dnl ambiguous cc.
			cc)
				fla_c_warning_flags=''
			;;
			dnl for all other C compilers.
			*)
				fla_c_warning_flags=''
			;;
		esac
	else

		dnl Set C compiler flags assuming we found...
		case $CC in
			dnl GNU gcc.
			gcc)
				fla_c_warning_flags='-w'
			;;
			dnl Intel cc.
			icc)
				fla_c_warning_flags='-w'
			;;
			dnl PathScale pathcc.
			pathcc)
				fla_c_warning_flags='-w'
			;;
			dnl PGI pgcc.
			pgcc)
				fla_c_warning_flags='-w'
			;;
			dnl NEC sxcc.
			sxcc)
				fla_c_warning_flags='-w none'
			;;
			dnl IBM xlc.
			xlc)
				fla_c_warning_flags='-w -qcpluscmt'
			;;
			dnl ambiguous cc.
			cc)
				fla_c_warning_flags=''
			;;
			dnl for all other C compilers.
			*)
				fla_c_warning_flags=''
			;;
		esac
	fi
	
	dnl Output the result.
	AC_MSG_RESULT([$fla_c_warning_flags])
	
	dnl Substitute the warning flags into the autoconf output files
	AC_SUBST(fla_c_warning_flags)

])
