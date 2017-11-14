#!/usr/bin/env bash

#
# export-for-release.sh
#
# Field G. Van Zee
#

print_usage()
{
	local script_name
	
	# Get the script name
	script_name=${0##*/}
	
	# Echo usage info
	echo " "
	echo " "$script_name
	echo " "
	echo " Field G. Van Zee"
	echo " "
	echo " Export a copy of the FLAME source tree and package it for release"
	echo " under the version vers_str. Note:"
	echo "  o The host argument must be a UTCS host since the UTCS network"
	echo "    filesystem is where the repository resides."
	echo "  o The release_dir argument should be the absolute path to the"
	echo "    directory where the releases are kept."
	echo " "
	echo " Usage:"
	echo "   ${script_name} [options] -v vstr host release_dir"
	echo " "
	echo " The following arguments are accepted/required:"
	echo " "
	echo "   -r n        revision"
	echo "                 Grab revision n from the repository. If no revision is"
	echo "                 given, the HEAD is exported."
	echo "   -v vstr     version"
	echo "                 Use version vstr to name the revision being exported for"
	echo "                 release. Some examples are \"0.9\" and \"1.0.3b2\"".
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	repos_path=/u/field/repos/flame
	toplevel_dir=trunk
	release_prefix=libflame
	rev_test_filepath=${toplevel_dir}/configure
	linux_revision_filepath=revision
	windows_revision_filepath=windows/revision
	tarball_ownership=':flame'
	tarball_perms=0664
	
	# Process our command line options.
	while getopts ":r:v:" opt; do
		case $opt in
			r  ) rev_num_str=$OPTARG ;;
			v  ) version_str=$OPTARG ;;
			\? ) print_usage
		esac
	done
	shift $(($OPTIND - 1))
	
	# Check the number of arguments after command line option processing.
	if [ $# != "2" ]; then
		print_usage
	fi
	
	# Set the default revision if necessary.
	if [ "$rev_num_str" == "" ]; then
		rev_num_str="HEAD"
	fi
	
	# Make sure that verboseness level is valid
	if [ "$version_str" == "" ]; then
		echo "$0: You must provide a version string with -v."
		exit 1
	fi
	
	# Extract our arguments.
	host=$1
	release_dirpath=$2

	svn_url="svn+ssh://${host}${repos_path}"
	
	# Come up with the revision string if we don't already have it.
	if [ "${rev_num_str}" == "HEAD" ]; then
		
		# Construct the subversion URL to the test file of the HEAD revision.
		svn_url_revtest="${svn_url}/${rev_test_filepath}"
		
		# Acquire the actual revision number for the test file of the HEAD.
		rev_num_str=$(svn info ${svn_url_revtest} | grep 'Revision' \
		                                          | sed 's/Revision: //g')
	fi

	# Construct release directory name.
	export_dirname="${release_prefix}-${version_str}"
	
	# Construct the tarball and zip package names.
	tarball=${export_dirname}.tar.gz
	zipball=${export_dirname}.zip
	
	# Execute the export.
	svn export -r ${rev_num_str} ${svn_url}/${toplevel_dir} ${export_dirname}
	
	# Place the revision number in the top-level Linux and Windows directories.
	echo "${rev_num_str}" > ${export_dirname}/${linux_revision_filepath}
	echo "${rev_num_str}" > ${export_dirname}/${windows_revision_filepath}
	
	# Tar/gzip (for Linux).
	tar czf ${tarball} ${export_dirname}
	
	# Zip (for Windows). We zip from within the directory rather than outside
	# because Windows, when extracting the zip file, likes to place the archive
	# contents in a directory of the same name as the zip archive, which would
	# result in the top-level files residing in a directory structure of the
	# form libflame-rxxxx/libflame-rxxxx, which is not what we want.
	cd ${export_dirname}
	zip -rq ../${zipball} *
	cd ..

	# Remove the old directory.
	rm -rf ${export_dirname}
	
	# Adjust the files' ownership and permissions.
	chown ${tarball_ownership} ${tarball} ${zipball}
	chmod ${tarball_perms}     ${tarball} ${zipball}
	
	# Move the tarball and zip package to the FLAME web directory.
	mv ${tarball} ${release_dirpath}
	mv ${zipball} ${release_dirpath}

	# Exit peacefully
	return 0
}


# The script's main entry point, passing all parameters given.
main "$@"
