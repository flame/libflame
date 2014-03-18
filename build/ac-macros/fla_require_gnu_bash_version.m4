AC_DEFUN([FLA_REQUIRE_GNU_BASH_VERSION],
[
	dnl Require that GNU bash already be present.
	AC_REQUIRE([FLA_REQUIRE_GNU_BASH])

	dnlfla_tmp_file="fla_require_gnu_bash_version_$1.tmp";
	fla_tmp_file=fla_bash.tmp;

	dnl Output a message saying we're going to check that the version of
	dnl GNU bash is equal to or later than the argument.
	AC_MSG_CHECKING([[for GNU bash >= $1]])
		
	dnlbash --version | grep 'version' | sed -e 's;.*version.;;' | sed -e 's;\..*$;;' 2> /dev/null 1> $fla_tmp_file
	bash --version | grep bash | dd skip=18 bs=1 count=3 2> /dev/null 1> $fla_tmp_file
	dnlecho -n "3.0" > $fla_tmp_file
	sh -c "echo ' >= '$1" >> $fla_tmp_file
	sh -c "echo 'halt'" >> $fla_tmp_file
	dnl dnlbroken fla_bash_version=`bash --version | grep 'ver' | sed -e 's;.*version.\([0-9\.][0-9\.]*\).*;\1;'`
	fla_result=`bc -q $fla_tmp_file`

	if test $fla_result = "1" ; then
		_cv_gnu_bash_version=yes;
		AC_MSG_RESULT([[yes]]);
	else
		AC_MSG_RESULT([[no]]);
		AC_MSG_ERROR([[Could not find GNU bash version >= $1! Bailing out...]],[1])
	fi
	rm -f $fla_tmp_file 

])
