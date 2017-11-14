#!/usr/bin/env bash

#
# mirror-tree.sh
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
	echo " Recursively descends through the directory given in argument 1 while"
	echo " creating a symmetric directory structure in the new directory specified"
	echo " by argument 2, ignoring regular files along the way."
	echo " "
	echo " Usage:"
	echo "   ${script_name} [-v] existing_dir new_mirror_dir"
	echo " "
	echo " The following options are accepted:"
	echo " "
	echo "   -v    verbose"
	echo "           Echo progress as directories are recursively created."
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}

main()
{
	# Process our command line options. We only respond to the -v flag,
	# in which case we'll echo what we're doing as we go along.
	while getopts ":v" opt; do
		case $opt in
			v  ) verbose_flag=1 ;;
			\? ) print_usage
			     exit 1
		esac
	done
	shift $(($OPTIND - 1))
	
	
	# Check the number of arguments after command line option processing.
	if [ $# != "2" ]; then
		print_usage
	fi
	
	
	# Extract arguments.
	e_dir=$1
	n_dir=$2
	
	
	# If the root new directory does not exist, then create it.
	if [ ! -d $n_dir ]; then
		
		# Be verbose, if -v was one of the command line options.
		if [ -n "$verbose_flag" ]; then
			echo "Creating $n_dir"
		fi
		
		
		# Make the root new directory. Create the parent directories if
		# they do not exist with the -p option.
		mkdir -p $n_dir
	fi
	
	
	# Initialize the recursive variables. We keep a separate variable
	# for the existing and new directories because they have different
	# roots, but they will always change in parallel.
	cur_e_dir=$e_dir
	cur_n_dir=$n_dir
	
	
	# Begin recursion, starting with the contents of the existing
	# directory.
	mirror_tree "$(ls $e_dir)"
	
	
	# Exit peacefully.
	return 0
}

mirror_tree()
{
	# Extract arguments.
	dir_contents="$1"
	
	# Process each item in our argument list (ie: each item in cur_e_dir).
	for thing in ${dir_contents}; do
		
		# Adjust the current existing and new directory paths to
		# include the current instance of thing.
		cur_e_dir="$cur_e_dir/$thing"
		cur_n_dir="$cur_n_dir/$thing"
		
		
		# If the current existing directory exists, then create a
		# corresponding subdirectory in new directory.
		if [ -d ${cur_e_dir} ]; then
			
			# Be verbose, if -v was one of the command line options.
			if [ -n "$verbose_flag" ]; then
				echo "Creating $cur_n_dir"
			fi
			
			
			# Make the new subdirectory, but only if it doesn't
			# already exist.
			if [ ! -d $cur_n_dir ]; then
				mkdir $cur_n_dir
			fi
			
			
			# Continue recursively on the contents of cur_e_dir.
			mirror_tree "$(ls $cur_e_dir)"
		fi
		
		# Delete the end of the path, up to the first / character to
		# prepare for the next "thing" in $@.
		cur_e_dir=${cur_e_dir%/*}
		cur_n_dir=${cur_n_dir%/*}
	done
}

# The script's main entry point, passing all parameters given.
main "$@"

