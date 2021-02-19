AC_DEFUN([FLA_SET_C_OPT_FLAGS],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
	
	AC_MSG_CHECKING([for (guessing) appropriate ${CC_VENDOR} optimization flags])

	if test "$1" == "yes" ; then

		dnl Set C compiler flags assuming we found...
		case ${CC_VENDOR} in
			dnl GNU gcc.
			gcc)
				fla_c_opt_flags='-O3'
			;;
			dnl LLVM clang.
			clang)
				fla_c_opt_flags='-O3'
			;;
			dnl Intel cc.
			icc)
				fla_c_opt_flags='-O3'
			;;
			dnl PathScale pathcc.
			pathcc)
				fla_c_opt_flags='-O3'
			;;
			dnl PGI pgcc.
			pgcc)
				fla_c_opt_flags='-O3'
			;;
			dnl NEC sxcc.
			sxcc)
				fla_c_opt_flags='-C vopt -O nooverlap -pvctl,noassume,vwork=stack'
			;;
			dnl IBM xlc.
			*xlc*)
				fla_c_opt_flags='-O3'
			;;
			dnl ambiguous cc.
			cc)
				fla_c_opt_flags='-O'
			;;
			dnl for all other C compilers.
			*)
				fla_c_opt_flags=''
			;;
		esac
	else
		
		dnl Set C compiler flags assuming we found...
		case ${CC_VENDOR} in
			dnl GNU gcc.
			gcc)
				fla_c_opt_flags='-O0'
			;;
			dnl LLVM clang.
			clang)
				fla_c_opt_flags='-O0'
			;;
			dnl Intel cc.
			icc)
				fla_c_opt_flags='-O0'
			;;
			dnl PathScale pathcc.
			pathcc)
				fla_c_opt_flags='-O0'
			;;
			dnl PGI pgcc.
			pgcc)
				fla_c_opt_flags='-O0'
			;;
			dnl NEC sxcc.
			sxcc)
				fla_c_opt_flags='-C noopt'
			;;
			dnl IBM xlc.
			xlc)
				fla_c_opt_flags=''
			;;
			dnl ambiguous cc.
			cc)
				fla_c_opt_flags=''
			;;
			dnl for all other C compilers.
			*)
				fla_c_opt_flags=''
			;;
		esac
	fi

	dnl Output the result.
	AC_MSG_RESULT([$fla_c_opt_flags])
	
	dnl Substitute the optimization flags into the autoconf output files
	AC_SUBST(fla_c_opt_flags)

])
