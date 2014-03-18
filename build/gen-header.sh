#!/bin/bash

#
# gen-header.sh
#
# Field G. Van Zee
#
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
	echo " Automatically updates the C .h header file for the source code"
	echo " in a given directory. The .h file must already exist--though the"
	echo " contents are not used."
	echo " "
	echo " Usage:"
	echo "   ${script_name} tmpl_dir src_dir"
	echo " "
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	local script_name
	local tmpl_dir
	local src_dir
	local header_filename
	local tmpl_header_filepath
	local src_header_filepath
	local fla_function_regexp
	local fla_prototypes_list
	
	
	# Check the number of arguments after command line option processing.
	if [ $# != "2" ]; then
		print_usage
	fi
	
	
	# Extract arguments.
	script_name=${0##*/}
	tmpl_dir="$1"
	src_dir="$2"
	
	
	# Check existence of the directories provided through the command-line.
	if [ ! -d "${tmpl_dir}" ]; then
		echo "${script_name}: Template directory does not exist (${tmpl_dir})."
		exit 1
	fi
	if [ ! -d "${src_dir}" ]; then
		echo "${script_name}: Source directory does not exist (${src_dir})."
		exit 1
	fi
	
	
	# Determine the template and destination filepaths to the of .h files. (Here we decide that
	# the destination will be the source directory.)
	tmpl_header_filename="FLA_alg.h"
	tmpl_header_filepath="${tmpl_dir}/${tmpl_header_filename}"
	
	# Check existance of the template header file.
	if [ ! -f "${tmpl_header_filepath}" ]; then
		echo "${script_name}: Failed to find template header file (${tmpl_header_filename}) in template directory (${tmpl_dir})."
		exit 1
	fi
	
	
	# Acquire the name of the local .h file. This assumes we have a .h file present in src_dir.
	# Then strip the path to leave just the filename.
	src_header_filepath="$(ls ${src_dir}/*.h 2> /dev/null | head -1)"
	src_header_filename=${src_header_filepath##*/}
	
	
	# We assume we have a .h file present in src_dir. Without a .h file, we can't proceed.
	if [ ! -f "${src_header_filepath}" ]; then
		echo "${script_name}: An FLA_*.h header file must exist in source directory (${src_dir})."
		exit 1
	fi
	
	
	# Acquire a list of all FLAME functions in the local .c files. This assumes the functions
	# begin with 'FLA_'.
	fla_function_regexp='^int FLA_'
	inv_glue_regexp='_glue'
	#fla_prototypes_list=$(grep "${fla_function_regexp}" ${src_dir}/*.c | grep FLA_Obj | grep -v ${inv_glue_regexp} | cut -d: -f2 | sed s/')'/');'/g)
	fla_prototypes_list=$(grep "${fla_function_regexp}" ${src_dir}/*.c | grep FLA_Obj | grep -v ${inv_glue_regexp} | cut -d: -f2 | sed s/')'[[:space:]]*$/');'/g)
	
	
	# First replace the old header file with a fresh template, then create the new header by 
	# appending the prototypes to the new template header file.	
	cp      ${tmpl_header_filepath}     ${src_header_filepath}
	echo    "${fla_prototypes_list}" >> ${src_header_filepath}
	echo -e ""                       >> ${src_header_filepath}
	
	
	# Exit peacefully.
	return 0
}


main "$@"
