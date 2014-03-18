AC_DEFUN([FLA_CHECK_ENABLE_VECTOR_INTRINSICS],
[
	dnl Initialize some variables.
	fla_enable_vector_intrinsics=no
	fla_vector_intrinsic_type=none
	
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested use of vector intrinsics])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([vector-intrinsics],
	              AC_HELP_STRING([--enable-vector-intrinsics=type],[Enable highly-optimized code that relies upon vector intrinsics to specify certain operations at a very low level. Valid values for type are "sse" and "none". Specifying "none" is the same as disabling the option. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "sse" ; then
			
			dnl Enable with SSE support.
			fla_enable_vector_intrinsics=yes
			fla_vector_intrinsic_type=sse

		elif test "$enableval" = "no" ; then
			
			dnl Disable SSE support.
			fla_enable_vector_intrinsics=no
			fla_vector_intrinsic_type=none

		elif test "$enableval" = "none" ; then
			
			dnl Disable SSE support.
			fla_enable_vector_intrinsics=no
			fla_vector_intrinsic_type=none

		else
			
			dnl Invalid option.
			AC_MSG_ERROR([[Invalid option to --enable-vector-intrinsics. Valid options are "sse" and "none".]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_vector_intrinsics=no
		fla_vector_intrinsic_type=none
	]
	)
	
	dnl Output the result.
	AC_MSG_RESULT([$fla_enable_vector_intrinsics])
		
	dnl Now act according to whether the option was requested.
	if test "$fla_enable_vector_intrinsics" = "yes" ; then
		
		dnl Tell the user we're checking the value given.
		AC_MSG_CHECKING([user-requested vector intrinsic type])
		AC_MSG_RESULT([$fla_vector_intrinsic_type])

		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_VECTOR_INTRINSICS,1,
		          [Determines whether vector intrinsics are used in certain low-level functions.])

		dnl Now we set flags and whatnot related to each vector intrinsic type.
		if test "$fla_vector_intrinsic_type" = "sse" ; then
		
			dnl Encode the C prepropcessor value for SSE vector intrinsics.
			fla_type_val=3
			
			dnl Determine the OpenMP flags to use.
			FLA_SET_C_SSE_FLAGS()
		
			dnl Also make sure we can find pmmintrin.h. If pmmintrin.h is present, do
			dnl nothing. If not, we're in trouble.
			dnl AC_CHECK_HEADER([pmmintrin.h],
			dnl [],
			dnl [
			dnl 	dnl Output an error.
			dnl 	AC_MSG_ERROR([SSE vector intrinsic support requires pmmintrin.h header file!])
			dnl ])
		fi

	else

		dnl Encode the C prepropcessor value for no vector intrinsics.
		fla_type_val=0

	fi

	dnl Define the preprocessor macro VECTOR_INTRINSICS_TYPE to the value
	dnl corresponding to SSE, or none, depending on how fla_type_val was
	dnl set above.
	AC_DEFINE_UNQUOTED(FLA_VECTOR_INTRINSIC_TYPE,$fla_type_val,
	                   [Encodes the type of vector intrinsics requested at configure-time.])

	dnl Substitute the output variables.
	AC_SUBST(fla_enable_vector_intrinsics)
	AC_SUBST(fla_vector_intrinsic_type)
])
