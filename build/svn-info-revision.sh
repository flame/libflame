#!/usr/bin/env bash

#
# svn-info-revision.sh
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
	echo " "${script_name}
	echo " "
	echo " Field G. Van Zee"
	echo " "
	echo " Query the subversion revision number for a file using 'svn info'."
	echo " This script is designed for use with doxygen."
	echo " "
	echo " Usage:"
	echo "   ${script_name} filepath"
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	local filepath
	
	# The name of the script
	script_name=${0##*/}

	# Check the number of arguments.
	if [ $# != "1" ]; then
		print_usage
	fi
	
	# Extract arguments.
	filepath="$1"
	
	# Check that old_dir actually exists and is a directory
	if [ ! -f ${filepath} ]; then
		echo "Hmm, ${filepath} does not seem to be a valid regular file."
		exit 1
	fi
	
	# Acquire the revision number for the file given.
	revision=$(svn info ${filepath} | grep 'Revision' | sed 's/Revision: //g')
	
	# Output the revision number.
	echo "r${revision}"

	# Exit peacefully.
	return 0
}

main "$@"
