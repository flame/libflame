#!/bin/bash

#
# query-revision-file.sh
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
	echo " Query the revision string for the HEAD of the repository."
	echo " "
	echo " Usage:"
	echo "   ${script_name}"
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	# Assume we are being run at the top-level directory.
	host=flame.cs.utexas.edu
	repos_dirpath=/u/field/repos/flame
	test_filepath=trunk/configure

	# Check the number of arguments after command line option processing.
	if [ $# != "0" ]; then
		print_usage
	fi
	
	# Construct the subversion URL to the FLAME repository.
	svn_url_revtest="svn+ssh://${host}/${repos_dirpath}/${test_filepath}"

	# Acquire the revision number for the test file of the HEAD.
	rev_num_str=$(svn info ${svn_url_revtest} | grep 'Revision' \
	                                          | sed 's/Revision: //g')

	# Output the revision number.
	echo "${rev_num_str}"

	# Exit peacefully.
	return 0
}


# The script's main entry point, passing all parameters given.
main "$@"
