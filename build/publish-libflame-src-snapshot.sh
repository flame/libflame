#!/bin/bash

#
# publish-libflame-src-snapshot.sh
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
	echo " Export a copy of the libflame source tree and publish it to the"
	echo " FLAME project's webpage. Note:"
	echo "  o The host argument must be a UTCS host since the UTCS network"
	echo "    filesystem is where the repository resides."
	echo "  o The snapshot_dir argument should be the absolute path to the"
	echo "    directory where the subversion snapshots are kept."
	echo " "
	echo " Usage:"
	echo "   ${script_name} [options] host snapshot_dir"
	echo " "
	echo " The following arguments are accepted/required:"
	echo " "
	echo "   -r n        revision"
	echo "                 Grab revision n from the repository. If no revision"
	echo "                 is given, the HEAD is exported."
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	repos_path=/u/field/repos/flame
	toplevel_dirpath=trunk
	export_prefix=libflame
	rev_test_filepath=${toplevel_dirpath}/configure
	linux_revision_filepath=revision
	windows_revision_filepath=windows/revision
	tarball_ownership=':flame'
	tarball_perms=0664
	max_num_snapshots=10  # five tar.gz/zip pairs

	# Make sure we have the svn binary.
	export PATH=$PATH:/lusr/bin
	
	# Process our command line options.
	while getopts ":r:" opt; do
		case $opt in
			r  ) rev_num_str=$OPTARG ;;
			\? ) print_usage
		esac
	done
	shift $(($OPTIND - 1))
	
	# Check the number of arguments after command line option processing.
	if [ $# != "2" ]; then
		print_usage
	fi
	
	# Set the default revision if necessary.
	if [ "${rev_num_str}" == "" ]; then
		rev_num_str="HEAD"
	fi
	
	# Extract our arguments.
	host=$1
	snapshot_dir_path=$2

	# Construct the subversion URL to the FLAME repository.
	svn_url="svn+ssh://${host}${repos_path}"
	
	# Come up with the revision string if we don't already have it.
	if [ "${rev_num_str}" == "HEAD" ]; then
		
		# Construct the subversion URL to the test file of the HEAD revision.
		svn_url_revtest="${svn_url}/${rev_test_filepath}"
		
		# Acquire the actual revision number for the test file of the HEAD.
		rev_num_str=$(svn info ${svn_url_revtest} | grep 'Revision' \
		                                          | sed 's/Revision: //g')
	fi

	# Construct export directory name.
	export_dirname="${export_prefix}-r${rev_num_str}"
	
	# Construct the tarball and zip package names.
	tarball=${export_dirname}.tar.gz
	zipball=${export_dirname}.zip
	
	# Execute the export.
	svn export -r ${rev_num_str} ${svn_url}/${toplevel_dirpath} ${export_dirname} 

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
	mv ${tarball} ${snapshot_dir_path}
	mv ${zipball} ${snapshot_dir_path}

	# Declare and initialize an integer count variable
	declare -i count
	count=0
	
	# Get a reversed list of the current items in the snapshots directory.
	existing_snapshots=$(ls -1r ${snapshot_dir_path})
	
	# Iterate over the list and keep only the most recent max_num_snapshots
	# snapshots, deleting the rest.
	for item in $existing_snapshots; do
		
		if [ $count -ge ${max_num_snapshots} ]; then
			rm -f ${snapshot_dir_path}/${item}
		fi
		
		let count=${count}+1
	done
	
	# Exit peacefully.
	return 0
}


# The script's main entry point, passing all parameters given.
main "$@"
