AC_DEFUN([FLA_SET_C_PROF_FLAGS],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
	
	AC_MSG_CHECKING([for (guessing) appropriate $CC profiling flags])

	if test "$1" == "yes" ; then
		
		dnl Set C compiler flags assuming we found...
		case $CC in
			dnl GNU gcc.
			gcc)
				fla_c_prof_flags='-pg'
			;;
			dnl Intel cc.
			icc)
				fla_c_prof_flags='-p'
			;;
			dnl PathScale pathcc.
			pathcc)
				fla_c_prof_flags='-pg'
			;;
			dnl PGI pgcc.
			pgcc)
				fla_c_prof_flags='-pg'
			;;
			dnl NEC sxcc.
			sxcc)
				fla_c_prof_flags='-p'
			;;
			dnl IBM xlc.
			*xlc*)
				fla_c_prof_flags='-pg'
			;;
			dnl ambiguous cc.
			cc)
				fla_c_prof_flags='-pg'
			;;
			dnl for all other C compilers.
			*)
				fla_c_prof_flags=''
			;;
		esac
	else

		dnl Set C compiler flags assuming we found...
		case $CC in
			dnl GNU gcc.
			gcc)
				fla_c_prof_flags=''
			;;
			dnl Intel cc.
			icc)
				fla_c_prof_flags=''
			;;
			dnl PathScale pathcc.
			pathcc)
				fla_c_prof_flags=''
			;;
			dnl PGI pgcc.
			pgcc)
				fla_c_prof_flags=''
			;;
			dnl NEC sxcc.
			sxcc)
				fla_c_prof_flags='-Np'
			;;
			dnl IBM xlc.
			xlc)
				fla_c_prof_flags=''
			;;
			dnl ambiguous cc.
			cc)
				fla_c_prof_flags=''
			;;
			dnl for all other C compilers.
			*)
				fla_c_prof_flags=''
			;;
		esac
	fi
	
	dnl Output the result.
	AC_MSG_RESULT([$fla_c_prof_flags])
	
	dnl Substitute the prof flags into the autoconf output files
	AC_SUBST(fla_c_prof_flags)

])
