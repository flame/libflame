AC_DEFUN([FLA_CHECK_ENABLE_MULTITHREADING],
[
	dnl Initialize some variables.
	fla_enable_multithreading=no
	fla_multithreading_model=none
	
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested support for multithreading])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([multithreading],
	              AC_HELP_STRING([--enable-multithreading=model],[Enable multithreading support. Valid values for model are "pthreads" and "openmp". Multithreading must be enabled to access SMP parallelized implementations. If enabled, use the environment variable FLA_NUM_THREADS to specify the number of threads to parallelize computations for LAPACK interfaces (FLA external interfaces). (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "openmp" ; then
			
			dnl Enable with OpenMP support.
			fla_enable_multithreading=yes
			fla_multithreading_model=openmp

		elif test "$enableval" = "pthreads" ; then
			
			dnl Enable with POSIX threads support.
			fla_enable_multithreading=yes
			fla_multithreading_model=pthreads

		elif test "$enableval" = "no" ; then
			
			dnl Disable multithreading.
			fla_enable_multithreading=no
			fla_multithreading_model=none

		else
			
			dnl Invalid option.
			AC_MSG_ERROR([[Invalid option to --enable-multithreading. Valid options are "openmp", "pthreads", and "no".]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_multithreading=no
		fla_multithreading_model=none
	]
	)
	
	dnl Output the result.
	AC_MSG_RESULT([$fla_enable_multithreading])
		
	dnl Now act according to whether the option was requested.
	if test "$fla_enable_multithreading" = "yes" ; then
		
		dnl Tell the user we're checking the value given.
		AC_MSG_CHECKING([user-requested multithreading model])
		AC_MSG_RESULT([$fla_multithreading_model])

		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_MULTITHREADING,1,
		          [Determines whether thread-specific blocks of code should be compiled.])

		dnl Now we set flags and whatnot related to each multithreading model.
		if test "$fla_multithreading_model" = "openmp"; then
		
			dnl Encode the C prepropcessor value for OpenMP.
			fla_model_val=1
			
			dnl Determine the OpenMP flags to use.
			FLA_SET_C_OPENMP_FLAGS()
		
			dnl Also make sure we can find omp.h. If omp.h is present, do
			dnl nothing. If not, we're in trouble.
			AC_CHECK_HEADER([omp.h],
			[],
			[
				dnl Output an error.
				AC_MSG_ERROR([OpenMP support requires omp.h header file!])
			])

		elif test "$fla_multithreading_model" = "pthreads" ; then

			dnl Encode the C prepropcessor value for POSIX threads.
			fla_model_val=2
			
			dnl Also make sure we can find pthread.h. If pthread.h is present, do
			dnl nothing. If not, we're in trouble.
			AC_CHECK_HEADER([pthread.h],
			[],
			[
				dnl Output an error.
				AC_MSG_ERROR([POSIX threads support requires pthread.h header file!])
			])

		fi

	else

		dnl Encode the C prepropcessor value for no multithreading.
		fla_model_val=0

	fi

	dnl Define the preprocessor macro MULTITHREADING_MODEL to the value corresponding
	dnl to OpenMP or POSIX threads, or no multithreading, depending on how fla_model_val
	dnl was set above.
	AC_DEFINE_UNQUOTED(FLA_MULTITHREADING_MODEL,$fla_model_val,
	                   [Encodes the type of multithreading chosen at configure-time.])

	dnl Substitute the output variables.
	AC_SUBST(fla_enable_multithreading)
	AC_SUBST(fla_multithreading_model)
])
