#!/usr/bin/env bash

#
# lntexinputs.sh
#
# Field G. Van Zee
#
# Creates symbolic links in the current directory of the files passed in as 
# arguments. This is because it seems nearly impossible to set an environment
# variable from within the top-level Makefile. (In this case, TEXINPUTS.)
#

print_usage()
{
	script_name="$0"
	
	echo "Usage: $script_name src_list"
	
	exit 1
}

main()
{
	local src_list
	
	
	# Check the number of arguments after command line option processing.
	if [ $# != "1" ]; then
		print_usage
	fi
	
	
	# Extract arguments.
	src_list="$1"
	
	
	# Soft link each file into the current directory.
	for src_file in $src_list; do
		
		# Create a soft link of the $src_file in the current directory.
		ln -s $src_file .
	done
	
	
	# Exit peacefully.
	return 0
}


main "$@"
