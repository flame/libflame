#!/bin/bash

#
# gen-make-frag.sh
#
# Field G. Van Zee
#


print_usage()
{
	#local script_name
	
	# Get the script name
	#script_name=${0##*/}
	
	# Echo usage info
	echo " "
	echo " "$script_name
	echo " "
	echo " Field G. Van Zee"
	echo " "
	echo " Automatically generates makefile fragments for a given directory. "
	echo " "
	echo " Usage:"
	echo "   ${script_name} [options] template.mk directory"
	echo " "
	echo " The following options are accepted:"
	echo " "
	echo "   -d          dry-run"
	echo "                 Go through all the motions, but don't actually generate any"
	echo "                 makefile fragments."
	echo "   -r          recursive"
	echo "                 Also generate makefile fragments for subdirectories."
	echo "   -h          hide"
	echo "                 Hide the makefile fragments by prepending filenames with '.'."
	echo "   -v [0|1|2]  verboseness level"
	echo "                 level 0: silent  (no output)"
	echo "                 level 1: default (one line per directory)"
	echo "                 level 2: verbose (several lines per directory)."
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}







#
# gen_mkfile()
#
# Creates a single makefile fragment in a user-specified directory and adds
# any local source files found to a top-level Makefile variable.
#
gen_mkfile()
{
	# Local variable declarations
	local mkfile_frag_tmpl_path
	local mkfile_frag_var_name
	local src_file_suffixes
	local this_dir
	local mkfile_frag_tmpl_name 
	local mkfile_name 
	local mkfile_frag_path
	local curr_frag_dir 
	local curr_frag_path
	local local_src_files
	local sub_items
	local item_path
	local item_suffix
	local curr_frag_sub_dirs
	
	
	# Extract our arguments to local variables
	mkfile_frag_tmpl_path=$1
	mkfile_frag_var_name=$2
	src_file_suffixes="$3"
	this_dir=$4
	
	
	# Strip the leading path from the template makefile path to get its
	# simple filename. Hide the output makefile fragment filename, if
	# requested.
	mkfile_frag_tmpl_name=${mkfile_frag_tmpl_path##*/}
	if [ -n "$hide_flag" ]; then
		mkfile_frag_path=$this_dir/.$mkfile_frag_tmpl_name
	else
		mkfile_frag_path=$this_dir/$mkfile_frag_tmpl_name
	fi
	
	
	# Determine the directory in which the fragment will reside.
	curr_frag_path=$this_dir
	curr_frag_dir=${this_dir##*/}
	
	
	# Initialize the local source list to empty
	local_src_files=""
	
	# Get a listing of the items in $this_dir
	sub_items=$(ls $this_dir)
	
	# Generate a list of the source files we've chosen
	for item in $sub_items; do
		
		# Prepend the directory to the item to get a relative path
		item_path=$this_dir/$item
		
		# Acquire the item's suffix, if it has one
		item_suffix=${item_path##*.}
		
		# If the suffix matches, then add it to our list
		#if [ $item_suffix = $src_file_suffix ]; then
		if is_in_list $item_suffix "$src_file_suffixes"
		then
			local_src_files="$local_src_files $item"
		fi
	done
	
	# Delete the leading " " space character in the local source files list.
	local_src_files=${local_src_files##" "}
	
	
	# Initialize the fragment subdirectory list to empty
	curr_frag_sub_dirs=""
	
	# Capture the relative path listing of items in $this_dir.
	sub_items=$(ls $this_dir)
	
	# Determine the fragment's subdirectory names, if any exist
	for item in $sub_items; do
		
		# Prepend the directory to the item to get a relative path
		item_path=$this_dir/$item
		
		# If item is a directory, and it's not in the ignore list, descend into it.
		if [ -d $item_path ] && ! should_ignore $item; then
			curr_frag_sub_dirs=$curr_frag_sub_dirs" "$item
		fi
	done
	
	# Delete the leading " " space character in fragment's subdirectory list.
	curr_frag_sub_dirs=${curr_frag_sub_dirs##" "}
	
	
	# Be verbose, if level 2 was requested.
	if [ "$verbose_flag" = "2" ]; then
		echo "mkf frag tmpl path: $mkfile_frag_tmpl_path"
		echo "mkf frag path:      $mkfile_frag_path"
		echo "curr frag path:     $curr_frag_path"
		echo "curr frag dir:      $curr_frag_dir"
		echo "curr frag sub dirs: $curr_frag_sub_dirs"
		echo "local src files:    $local_src_files"
		echo "src file suffixes:  $src_file_suffixes"
		echo "mkf frag var name:  $mkfile_frag_var_name"
		echo "--------------------------------------------------"
	fi
	
	
	# Copy the template makefile to the directory given, using the new
	# makefile name we just created above.
	if [ -z "$dry_run_flag" ]; then
		cat $mkfile_frag_tmpl_path | sed -e s/"$mkfile_fragment_curr_dir_name_anchor"/"$curr_frag_dir"/g \
		                           | sed -e s/"$mkfile_fragment_sub_dir_names_anchor"/"$curr_frag_sub_dirs"/g \
		                           | sed -e s/"$mkfile_fragment_local_src_files_anchor"/"$local_src_files"/g \
		                           | sed -e s/"$mkfile_fragment_src_var_name_anchor"/"$mkfile_frag_var_name"/g \
		                           > $mkfile_frag_path
	fi
	
	
	# Return peacefully.
	return 0
}


#
# gen_mkfiles
#
# Recursively generates makefile fragments for a directory and all 
# subdirectories. All of the actual work happens in gen_mkfile().
#
gen_mkfiles()
{
	# Local variable declarations
	local item sub_items curr_dir this_dir
	
	
	# Extract our argument
	curr_dir=$1
	
	
	# Append a relevant suffix to the makefile variable name, if necesary
	all_add_src_var_name "$curr_dir"
	
	
	# Be verbose if level 2 was requested
	if   [ "$verbose_flag" = "2" ]; then
		echo ">>>" $script_name $mkfile_frag_tmpl_path ${src_var_name}_$SRC "\"$src_file_suffixes\"" $curr_dir
	elif [ "$verbose_flag" = "1" ]; then
		echo "$script_name: creating makefile fragment in $curr_dir"
	fi
	
	
	# Call our function to generate a makefile in the directory given.
	gen_mkfile $mkfile_frag_tmpl_path "${src_var_name}_$SRC" "$src_file_suffixes" $curr_dir
	
	
	# Get a listing of the directories in $directory
	sub_items=$(ls $curr_dir)
	
	# Descend into the contents of root_dir to generate the subdirectories'
	# makefile fragments.
	for item in $sub_items; do
		
		# If item is a directory, and it's not in the ignore list, descend into it.
		if [ -d "$curr_dir/$item" ] && ! should_ignore $item; then
			this_dir=$curr_dir/$item
			gen_mkfiles $this_dir
		fi
	done
	
	
	# Remove a relevant suffix from the makefile variable name, if necesary
	all_del_src_var_name "$curr_dir"
	
	
	# Return peacefully
	return 0
}



update_src_var_name_lib()
{
	local dir act i name var_suffix
	
	
	# Extract arguments
	act="$1"
	dir="$2"
	
	
	# Strip / from end of directory path, if there is one, and then strip
	# path from directory name.
	dir=${dir%/}
	dir=${dir##*/}
	
	
	# Run through our list
	for i in ${lib_i[@]}; do
		
		# Get the ith name
		name=${lib_name[$i]}
		
		# If the current item matches $dir, then we'll probably have to make
		# a modification of some form to src_var_name.
		if [ "$dir" = "$name" ]; then 
			
			# Get the suffix in uppercase.
			var_suffix=$(echo "$name" | tr '[:lower:]' '[:upper:]')
			
			# Either add or remove the suffix.
			if [ "$act" == "+" ]; then
				
				# This conditional is added so that only the first directory
				# matching an item in the lib_list is appended to the 
				# src_var_name variable. Otherwise we might have source in 
				# base/flamec/wrappers/blas/3/gemm being appended to a 
				# variable named MK_BASE_FLAMEC_BLAS_SRC, which is not what
				# we want.
				if [ "$src_var_name" = "MK" ]; then
					src_var_name=${src_var_name}_$var_suffix
				else
					continue
				fi
			else
				src_var_name=${src_var_name%_$var_suffix}
			fi
			
			# No need to continue iterating.
			break;
		fi
	done
}
update_src_var_name_leaf()
{
	local dir act i name var_suffix
	
	
	# Extract arguments
	act="$1"
	dir="$2"
	
	
	# Strip / from end of directory path, if there is one, and then strip
	# path from directory name.
	dir=${dir%/}
	dir=${dir##*/}
	
	
	# Run through our list
	for i in ${leaf_i[@]}; do
		
		# Get the ith name
		name=${leaf_name[$i]}
		
		# If the current item matches $dir, then we'll have
		# to make a modification of some form.
		if [ "$dir" = "$name" ]; then
			
			# Convert the variable suffix to uppercase.
			var_suffix=$(echo "$name" | tr '[:lower:]' '[:upper:]')
			
			# Get the valid source file suffixes from the leaf array.
			file_suffix_list=${leaf_suffix[$i]}
			
			# Either add or remove the suffix, and also update the
			# source file suffix variable.
			if [ "$act" == "+" ]; then
				src_var_name=${src_var_name}_$var_suffix
				src_file_suffixes="$file_suffix_list"
			else
				src_var_name=${src_var_name%_$var_suffix}
				src_file_suffixes=$no_file_suffix
			fi
			
			# No need to continue iterating.
			break;
		fi
	done
}

init_src_var_name()
{
	local dir="$1"
	
	# Strip off the leading / if there is one
	dir=${dir%%/}
	
	# Convert the / directory separators into spaces to make a list of 
	# directories.
	list=${dir//\// }
	
	# Inspect each item in $list
	for item in $list; do
		
		# Try to initialize the source variable name
		all_add_src_var_name $item
	done
}

all_add_src_var_name()
{
	local dir="$1"
	
	update_src_var_name_lib  "+" "$dir"
	update_src_var_name_leaf "+" "$dir"

}

all_del_src_var_name()
{
	local dir="$1"
	
	update_src_var_name_leaf "-" "$dir"
	update_src_var_name_lib  "-" "$dir"
}

read_mkfile_var_config()
{
	local index lname lsuff
	declare -i count
	
	# Read each line of the file describing the library types that might be
	# built.
	count=0
	for i in $(cat "build/config/lib_list"); do
		
		# Get the index and library name for each line
		#index=${i%%:*}
		#lname=${i##*:}
		lname=${i}
		
		# Save this info into their respective arrays
		lib_i[$count]=$count
		lib_name[$count]=$lname
		
		# Increment the counter
		let count=$count+1
	done
	
	
	# Read each line of the file describing leaf node types
	count=0
	for i in $(cat "build/config/leaf_list"); do
		
		# Get the index, suffix, and directory name for each line
		#index=${i%%:*}
		lname=${i%%:*}
		#lname=${lname#*:}
		lsuff=${i##*:}
		lsuff=${lsuff//,/ }
		
		# Save this info into their respective arrays
		leaf_i[$count]=$count
		leaf_name[$count]=$lname
		leaf_suffix[$count]=$lsuff
		
		# Increment the counter
		let count=$count+1
	done
	
	
	# Read each line of the file describing directories to ignore
	count=0
	for i in $(cat "build/config/ignore_list"); do
	
		# Get the index and name for each line
		#index=${i%%:*}
		lname=${i}
		
		# Save this info into their respective arrays
		ignore_i[$count]=$count
		ignore_name[$count]=$lname
		
		# Increment the counter
		let count=$count+1
	done
}	

main()
{
	# Global array delcarations
	declare -a lib_i
	declare -a lib_name
	declare -a leaf_i
	declare -a leaf_name
	declare -a leaf_suffix
	declare -a ignore_i
	declare -a ignore_name
	
	
	# Define these makefile template "anchors" used in gen_mkfile()
	mkfile_fragment_curr_dir_name_anchor="_mkfile_fragment_curr_dir_name_"
	mkfile_fragment_sub_dir_names_anchor="_mkfile_fragment_sub_dir_names_"
	mkfile_fragment_local_src_files_anchor="_mkfile_fragment_local_src_files_"
	mkfile_fragment_src_var_name_anchor="_mkfile_fragment_src_var_name_"
	
	# The name of the script, stripped of any preceeding path.
	script_name=${0##*/}
	
	# The variable that always holds the string that will be passed to
	# gen_mkfile() as the source variable to insert into the fragment.mk.
	src_var_name='MK'
	
	# The suffix appended to all makefile fragment source variables
	SRC='SRC'
	
	# The placeholder we use to signify that we're not looking for any
	# source files. (Does any file format use .z?) 
	no_file_suffix='z'
	src_file_suffixes='z'
	
	# The arguments to this function. They'll get assigned meaningful
	# values after getopts.
	mkfile_frag_tmpl_path=""
	root_dir=""
	
	# Flags set by getopts.
	dry_run_flag=""	
	hide_flag=""
	recursive_flag=""
	verbose_flag=""
	
	
	# Local variable declarations
	local item sub_items this_dir
	
	
	
	# Read the makefile source variable config files to be used in the
	# makefile fragment generation.
	read_mkfile_var_config
	
	
	# Process our command line options.
	while getopts ":dhrv:" opt; do
		case $opt in
			d  ) dry_run_flag="1" ;;
			h  ) hide_flag="1" ;;
			r  ) recursive_flag="1" ;;
			v  ) verbose_flag=$OPTARG ;;
			\? ) print_usage
		esac
	done
	shift $(($OPTIND - 1))
	
	
	# Make sure that verboseness level is valid
	if [ "$verbose_flag" != "0" ] && 
	   [ "$verbose_flag" != "1" ] && 
	   [ "$verbose_flag" != "2" ]; then
		verbose_flag="1"
	fi
	
	# Check the number of arguments after command line option processing.
	if [ $# != "2" ]; then
		print_usage
	fi
	
	
	# Extract our arguments
	mkfile_frag_tmpl_path=$1
	root_dir=$2
	
	
	# Strip / from end of directory path, if there is one.
	root_dir=${root_dir%/}
	
	
	# Append relevant suffixes to the makefile variable name based on the
	# current root, if any of the directory names match.
	init_src_var_name "$root_dir"
	
	
	# Be verbose if level 2 was requested
	if   [ "$verbose_flag" = "2" ]; then
		echo ">>>" $script_name $mkfile_frag_tmpl_path ${src_var_name}_$SRC "\"$src_file_suffixes\"" $root_dir
	elif [ "$verbose_flag" = "1" ]; then
		echo "$script_name: creating makefile fragment in $root_dir"
	fi
	
	
	# Call our function to generate a makefile in the root directory given.
	gen_mkfile $mkfile_frag_tmpl_path "${src_var_name}_$SRC" "$src_file_suffixes" $root_dir
	
	
	# If we were asked to act recursively, then continue processing
	# root_dir's contents.
	if [ -n "$recursive_flag" ]; then
		
		# Get a listing of the directories in $directory
		sub_items=$(ls $root_dir)
		
		# Descend into the contents of root_dir to generate the makefile
		# fragments.
		for item in $sub_items; do
			
			# If item is a directory, and it's not in the ignore list, descend into it.
			if [ -d "$root_dir/$item" ] && ! should_ignore $item ; then
				
				this_dir=$root_dir/$item
				gen_mkfiles $this_dir
			fi
		done
	fi
	
	
	# Exit peacefully
	return 0
}

should_ignore()
{
	local item name
	
	
	# Extract argument, the item that we may have to ignore.
	item="$1"
	
	
	# Process each index in ignore array
	for i in "${ignore_i[@]}"; do
		
		# Get the ith name
		name=${ignore_name[$i]}
		
		# If the current value of $name matches $item, then we need to
		# signal to calling function that we SHOULD ignore $item.
		# Notice that returning zero value means "success".
		if [ "$item" = "$name" ]; then
			return 0
		fi
	done
	
	
	# If we got this far, then item is not in the ignore list, so we
	# signal that we should NOT ignore item. Notice that returning
	# a non-zero value means "failure".
	return 1
}

is_in_list()
{
	local cur_item the_item item_list
	
	# Extract argument.
	the_item="$1"
	item_list="$2"
	
	# Check each item in the list against the item of interest.
	for cur_item in ${item_list}; do
		
		# If the current item in the list matches the one of interest
		if [ "${cur_item}" = "${the_item}" ]; then
			
			# Return success (ie: item was found).
			return 0
		fi
	done
	
	# If we made it this far, return failure (ie: item not found).
	return 1
}

# The script's main entry point, passing all parameters given.
main "$@"
