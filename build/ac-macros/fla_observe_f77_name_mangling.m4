AC_DEFUN([FLA_OBSERVE_F77_NAME_MANGLING],
[
	dnl This is important. We have to invoke the AC_F77_WRAPPERS macro first so
	dnl it will invoke the name-mangling macro.
	AC_REQUIRE([AC_F77_WRAPPERS])
	
	fla_f77_foobar_unmangled='foobar'
	fla_f77_foobar_mangled=''

	dnl Ask the system to mangle a test name for us using the mangling that
	dnl will be used by the current Fortran environment.
	AC_F77_FUNC($fla_f77_foobar_unmangled,fla_f77_foobar_mangled)

	dnl Substitute the output variables.
	AC_SUBST(fla_f77_foobar_unmangled)
	AC_SUBST(fla_f77_foobar_mangled)
])
