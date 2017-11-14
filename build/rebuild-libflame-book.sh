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
	echo " Rebuild the libflame book from LaTeX source."
	echo " "
	echo " Note:"
	echo "  o The src_dir argument should be the absolute path to the"
	echo "    top-level libflame directory."
	echo "  o The docs/libflame subdirectory of src_dir will have \'svn update\'"
	echo "    run on it before the build is performed."
	echo " "
	echo " Usage:"
	echo "   ${script_name} src_dir"
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	# Make sure we have the svn binary.
	export PATH=$PATH:/lusr/bin
	
	# Make sure we have the pdflatex binary.
	export PATH=$PATH:/lusr/tex/bin
	
	# Check the number of arguments.
	if [ $# != "1" ]; then
		print_usage
	fi
	
	# Extract our arguments.
	toplevel_dirpath=$1
	libflame_dirpath=docs/libflame

	# Change to the top-level directory.
	cd ${toplevel_dirpath}
	
	# Run 'git pull' to make sure we have the latest source.
	git pull 2> /dev/null

	# Update the version file in case it is out-of-date.
	./build/update-version-file.sh

	# Change into libflame directory.
	cd ${libflame_dirpath}

	# Clean out the directory and then rebuild the document via make.
	make clean
	make 2> make.err
	
	# Exit peacefully.
	return 0
}


# The script's main entry point, passing all parameters given.
main "$@"
