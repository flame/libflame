AC_DEFUN([FLA_CHECK_ENABLE_BLAS3_FRONT_END_CNTL_TREES],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested enabling level-3 BLAS front end control trees])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([blas3-front-end-cntl-trees],
	              AS_HELP_STRING([--enable-blas3-front-end-cntl-trees],[Enable code that uses control trees to select a reasonable variant and blocksize when level-3 BLAS front-ends are invoked. When disabled, the front-ends invoke their corresponding external implementations. Note that control trees are always used for LAPACK-like operations. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_blas3_front_end_cntl_trees=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_blas3_front_end_cntl_trees=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in CHECK_ENABLE_BLAS3_FRONT_END_CNTL_TREES!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option
		dnl Default behavior is to disable this feature.
		fla_enable_blas3_front_end_cntl_trees=no
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_blas3_front_end_cntl_trees" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_BLAS3_FRONT_END_CNTL_TREES,1,
		          [Determines whether to use control trees to select a reasonable FLAME variant and blocksize when level-3 BLAS front-ends are invoked.])
		
	elif test "$fla_enable_blas3_front_end_cntl_trees" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_NON_CRITICAL_CODE!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_blas3_front_end_cntl_trees)

])
