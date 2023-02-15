AC_DEFUN([FLA_REQUIRE_REFERENCE_LAPACK],
[
	AC_REQUIRE([FLA_CHECK_ENABLE_LAPACK2FLAME])
	AC_REQUIRE([FLA_CHECK_ENABLE_LAPACK2FLASH])
	AC_REQUIRE([FLA_CHECK_ENABLE_LEGACY_LAPACK])

	AC_MSG_CHECKING([whether the files from the reference lapack tar are needed])
	
	script_name=${0##*/}
	path=${0%/${script_name}}

	dnl Unpack tar for reference lapack if needed.
	dnl if not, clear it out
	if test "$fla_enable_lapack2flame" = "yes" || test "$fla_enable_lapack2flash" = "yes" ; then
		if test "$fla_enable_legacy_lapack" = "yes" ; then
			dnl run script clean
			echo "	Using legacy f2c lapack files, clearing out any fortran files"
			${path}/src/map/common/lapacksrc/scripts/regen-files.sh cleanup
		else
			dnl run script build
			echo "	Using lapack fortran files, unpacking tar"
			${path}/src/map/common/lapacksrc/scripts/regen-files.sh build

		fi
	else
		dnl run script clean
		echo "	Not using reference lapack, clearing out any fortran files"
		${path}/src/map/common/lapacksrc/scripts/regen-files.sh cleanup
	fi

])
