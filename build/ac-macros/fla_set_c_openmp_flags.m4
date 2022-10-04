AC_DEFUN([FLA_SET_C_OPENMP_FLAGS],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
		
	dnl Echo C OpenMP flags to user.
	AC_MSG_CHECKING([for (guessing) OpenMP flags for ${CC_VENDOR}])
		
	dnl Set the OpenMP compiler flags based on which C compiler we're going to use.
	case ${CC_VENDOR} in
		dnl Intel icc.
		icc)
			fla_c_openmp_flags='-openmp'
		;;
		dnl GNU gcc.
		gcc)
			fla_c_openmp_flags='-fopenmp'
		;;
		dnl PathScale pathcc
		pathcc)
			fla_c_openmp_flags='-mp'
		;;
		dnl PGI pgcc.
		pgcc)
			fla_c_openmp_flags='-mp'
		;;
		dnl NEC sxcc.
		sxcc)
			fla_c_openmp_flags='-P openmp'
		;;
		dnl IBM xlc
		*xlc*)
			fla_c_openmp_flags='-qsmp=omp'
		;;
		dnl Anything ending in 'cc', including craycc.
		*cc)
			fla_c_openmp_flags='-fopenmp'
		;;
		dnl for all other C compilers.
		*)
			fla_c_openmp_flags='unknown'
		;;
	esac
	
	dnl Output the result.
	AC_MSG_RESULT([$fla_c_openmp_flags])
	
	dnl Check string in case C OpenMP compiler is unknown.
	if test "$fla_c_openmp_flags" = "unknown" ; then
		
		dnl Tell the user we can't continue unless we know what flags
		dnl to pass to the C compiler to enable OpenMP.
		AC_MSG_ERROR([configure doesn't know what flag to give ${CC_VENDOR} in order to enable OpenMP support. Please submit a bug report to the FLAME developers.])
	fi

	dnl Output the C OpenMP flags variable.
	AC_SUBST(fla_c_openmp_flags)
])
