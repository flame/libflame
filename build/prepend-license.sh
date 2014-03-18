#!/bin/bash

#
# prepend-license.sh
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
	echo " Recusively descends a source tree and prepends comment blocks that contain"
	echo " license info. The type of comment block to prepend will depend on file"
	echo " extension, since different files need different comment character sequences. "
	echo " "
	echo " Usage:"
	echo "   ${script_name} [options] template_dir source_dir"
	echo " "
	echo " The following options are accepted:"
	echo " "
	echo "   -d          dry-run"
	echo "                 Go through all the motions, but don't actually prepend the"
	echo "                 license info to any files."
	echo "   -r          recursive"
	echo "                 Also prepend license info to files of subdirectories."
	echo "   -v [0|1]    verboseness level"
	echo "                 level 0: silent  (no output)"
	echo "                 level 1: default (one line per directory)"
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}







prepend_license()
{
	# Local variable declarations
	local tmpl_dir_path
	local this_dir
	local sub_items
	local item
	local item_path
	local item_path_temp
	local item_suffix
	local tmpl_file_name
	local tmpl_file_path
	
	
	# Extract our arguments to local variables
	tmpl_dir_path=$1
	this_dir=$2
	
	
	# Get a listing of the items in $this_dir
	sub_items=$(ls $this_dir)
	
	
	# Generate a list of the source files we've chosen
	for item in $sub_items; do
		
		# Prepend the directory to the item to get a relative path
		item_path=$this_dir/$item
			
		# Only consider regular files (ie: not directories or other special files).
		if [ -f "$item_path" ]; then
		
			# Generate a temporary file name
			item_path_temp=$item_path.$temp_file_suffix
			
			# Acquire the item's suffix, if it has one
			item_suffix=${item_path##*.}
			
			# For the suffix found above, generate the name of the license template
			# file, if it existed.
			tmpl_file_name=$license_tmpl_root.$item_suffix
			tmpl_file_path=$tmpl_dir_path/$tmpl_file_name
			
			# If the file exists, then we can prepend the file to the current item.
			# Otherwise, we can spit out a warning if verbosity was requested.
			# Note that we avoid (via ! -h) touching symbolic links, though.
			if [   -f $tmpl_file_path ] &&
			   [ ! -h $item_path      ]; then
				
				# Be verbose, if requested, about which file we're looking at,
				# taking into account the dry-run flag.
				if   [ "$verbose_flag" != "0" ] && [ -n "$dry_run_flag" ]; then
					echo ">>> Prepending (dry run) license file to ${item_path}."
				elif [ "$verbose_flag" != "0" ] && [ -z "$dry_run_flag" ]; then
					echo ">>> Prepending license file to ${item_path}."
				fi
				
				# Only prepend if we're not doing a dry-run.
				if [ -z "$dry_run_flag" ]; then
					
					declare -i n_header_lines
					declare -i n_lines
					declare -i diff
					
					#n_header_lines=$(wc -l $tmpl_file_path | cut -d' ' -f1)
					#n_lines=$(wc -l $item_path | cut -d' ' -f1)
					#diff=${n_lines}-${n_header_lines}
					#tail -n ${diff}    $item_path > $item_path_temp
					##tail -n +32 $item_path > $item_path_temp
					#mv  $item_path_temp $item_path
					
					
					# Prepend to a temporary file, and then replace the original with the new file
					# that contains the license comment prepended to the beginning.
					cat $tmpl_file_path $item_path > $item_path_temp
					mv  $item_path_temp $item_path
				fi
			else
				
				# Be verbose, if requested, about the appropriate license file not existing.
				if [ "$verbose_flag" != "0" ]; then
					echo ">>> Could not find appropriate license file for ${item_path}. Skipping!"
				fi
			fi
		fi
	done
	
	
	# Return peacefully.
	return 0
}


prepend_licenses()
{
	# Local variable declarations
	local item sub_items curr_dir this_dir
	
	
	# Extract our argument
	curr_dir=$1
	
	
	# Call our function to prepend the license comments to the files in the current directory.
	prepend_license $tmpl_dir_path $curr_dir
	
	
	# Get a listing of the directories in $directory
	sub_items=$(ls $curr_dir)
	
	# Descend into the contents of curr_dir to prepend the license info to files residing
	# in subdirectories of the current directory.
	for item in $sub_items; do
		
		# If item is a directory, descend into it.
		if [ -d "$curr_dir/$item" ]; then
			
			this_dir=$curr_dir/$item
			prepend_licenses $this_dir
		fi
	done
	
	
	# Return peacefully
	return 0
}




main()
{
	
	# The arguments to this function. They'll get assigned meaningful
	# values after getopts.
	tmpl_dir_path=""
	
	# Flags set by getopts.
	dry_run_flag=""	
	recursive_flag=""
	verbose_flag=""
	
	# The root name of the template license files
	license_tmpl_root="license"

	# The suffix used to create temporary files
	temp_file_suffix="license_temp"
	
	
	# Local variable declarations
	local item sub_items this_dir
	
	
	# Process our command line options.
	while getopts ":drv:" opt; do
		case $opt in
			d  ) dry_run_flag="1" ;;
			r  ) recursive_flag="1" ;;
			v  ) verbose_flag=$OPTARG ;;
			\? ) print_usage
		esac
	done
	shift $(($OPTIND - 1))
	
	
	# Make sure that verboseness level is valid
	if [ "$verbose_flag" != "0" ] && 
	   [ "$verbose_flag" != "1" ]; then
		verbose_flag="1"
	fi
	
	# Check the number of arguments after command line option processing.
	if [ $# != "2" ]; then
		print_usage
	fi
	
	
	# Extract our arguments
	tmpl_dir_path=$1
	src_dir=$2
	
	
	# Strip / from end of directory paths, if there is one.
	tmpl_dir_path=${tmpl_dir_path%/}
	src_dir=${src_dir%/}
	
	
	# Call our function to prepend our license info to the files in the 
	# source directory given.
	prepend_license $tmpl_dir_path $src_dir
	
	
	# If we were asked to act recursively, then continue processing
	# src_dir's contents.
	if [ -n "$recursive_flag" ]; then
		
		# Get a listing of the directories in src_dir
		sub_items=$(ls $src_dir)
		
		# Descend into the contents of src_dir to prepend the license info.
		for item in $sub_items; do
			
			# If item is a directory, descend into it.
			if [ -d "$src_dir/$item" ]; then
				
				this_dir=$src_dir/$item
				prepend_licenses $this_dir
			fi
		done
	fi
	
	
	# Exit peacefully
	return 0
}


# The script's main entry point, passing all parameters given.
main "$@"
