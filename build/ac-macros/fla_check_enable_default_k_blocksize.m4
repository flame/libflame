AC_DEFUN([FLA_CHECK_ENABLE_DEFAULT_K_BLOCKSIZE],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested specific default blocksize in k dimension])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([default_k_blocksize],
	              [  --enable-default-k-blocksize=kb],
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" || test "$enableval" = "0" ; then
			
			dnl User provided --enable-<option>=no or --enable-<option>=0.
			AC_MSG_ERROR([[Invalid option to --enable-default-k-blocksize]])

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes.
			AC_MSG_ERROR([[Invalid option to --enable-default-k-blocksize]])
		else
			
			dnl User provided a valid argument (hopefully).
			fla_enable_default_k_blocksize=yes
			fla_default_k_blocksize=$enableval
		fi
		
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_default_k_blocksize=no
	]
	)

	dnl Now act according to whether the option was requested.
	if test "$fla_enable_default_k_blocksize" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Enable option by setting a corresponding preprocessor directive
		dnl to the requested value.
		AC_DEFINE_UNQUOTED(FLA_DEFAULT_K_BLOCKSIZE,$fla_default_k_blocksize,
		                   [Sets the default blocksize in the k dimension.])

		dnl Tell the user we're checking the value given.
		AC_MSG_CHECKING([user-requested k dimension blocksize])
		AC_MSG_RESULT([$fla_default_k_blocksize])
		
		dnl Substitute output variable values.
		AC_SUBST(fla_enable_default_k_blocksize)
		AC_SUBST(fla_default_k_blocksize)

	else

		dnl Output the result.
		AC_MSG_RESULT([no])
	fi
])
