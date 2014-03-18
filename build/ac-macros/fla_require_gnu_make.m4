AC_DEFUN([FLA_REQUIRE_GNU_MAKE],
[
	AC_CACHE_CHECK( for GNU make,_cv_gnu_make_command,
		_cv_gnu_make_command='';
		dnl Search all the common names for GNU make
		for a in "$MAKE" make gmake gnumake ; do
			if test -z "$a" ; then continue ; fi ;
			if ( sh -c "$a --version" 2> /dev/null | grep GNU  2>&1 > /dev/null ) ;  then
				_cv_gnu_make_command=$a;
				break;
			fi
		done;
	);
	dnl If there was a GNU version, then set @fla_gnu_make_found@ to "true"
	if test  "x$_cv_gnu_make_command" != "x"  ; then
		fla_gnu_make_found=yes
	else
		fla_gnu_make_found=no
		AC_MSG_RESULT("not found!");
		AC_MSG_ERROR([[Could not locate GNU make! Bailing out...]],[1])
	fi
	AC_SUBST(fla_gnu_make_found)
])
