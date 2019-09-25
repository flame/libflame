AC_DEFUN([FLA_OBSERVE_OS_NAME],
[
	os_name=$(uname -s)
	os_vers=$(uname -r)

	dnl os_name_esc=$(echo "${os_name}" | sed 's/\//\\\//g')

	dnl Substitute the OS name and version into autoconf output files
	AC_SUBST(os_name)
	AC_SUBST(os_vers)
])
