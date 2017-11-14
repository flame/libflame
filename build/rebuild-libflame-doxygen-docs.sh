#!/usr/bin/env bash

#
# rebuild-libflame-book.sh
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
	echo " Rebuild the libflame doxygen documentation."
	echo " Note:"
	echo "  o The toplevel_dir argument should be the absolute path to the"
	echo "    top-level directory where the libflame distribution is located."
	echo "  o The toplevel_dir directory will have \'svn update\' run on it"
	echo "    before the build is performed."
	echo " "
	echo " Usage:"
	echo "   ${script_name} toplevel_dir"
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	# Make sure we have the svn binary.
	export PATH=$PATH:/lusr/bin
	
	# Check the number of arguments.
	if [ $# != "1" ]; then
		print_usage
	fi
	
	# Extract our arguments.
	top_dirpath=$1

	# Change to the source directory.
	cd ${top_dirpath}
	
	# Run 'svn update' to make sure we have the latest source.
	svn update

	# Update the revision file in case it is out-of-date.
	./build/update-check-rev-file.sh
	
	# Update the revision number in the doxygen config file.
	revnum=$(cat revision)
	cat Doxyfile | sed "s/revision_anchor/$revnum/g" > Doxyfile.temp
	mv Doxyfile Doxyfile.orig
	mv Doxyfile.temp Doxyfile

	# Make the doxygen directory so we don't cause doxygen to state that
	# it's creating the directory for us.
	mkdir doxygen
	
	# Run doxygen.
	doxygen

	# Publish the new documentation.
	rm -rf /u/www/users/flame/libflame/docs/html
	mv doxygen/html /u/www/users/flame/libflame/docs
	
	# Restore the original file for next time.
	mv Doxyfile.orig Doxyfile

	# Remove the doxygen docs directory.
	rm -rf doxygen

	# Exit peacefully.
	return 0
}


# The script's main entry point, passing all parameters given.
main "$@"
