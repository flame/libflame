#!/usr/bin/env bash

#
# find-dependents.sh
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
	echo " Finds all dependents of a given routine in a specified source"
	echo " directory and outputs them on standard output. The routine may"
	echo " be specified instead as a quoted string of routines, in which"
	echo " case all dependents will be aggregated indistinguishably."
	echo " "
	echo " Usage:"
	echo "   ${script_name} [-fq] routine src_dir"
	echo " "
	echo " The following options are accepted:"
	echo " "
	echo "   -f    fast"
	echo "           By default, the script builds a hierarchy of dependents, which"
	echo "           may require searching for routines already found by a previous"
	echo "           search. This option allows the script to skip these redundant"
	echo "           searches, which may speed up execution substantially. The final"
	echo "           list of dependents is unaffected by this option."
	echo "   -q    quiet"
	echo "           Suppress the hierarchy of dependents that is output incrementally"
	echo "           as the script executes its recursive search."
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	local script_name
	local the_routine
	local src_dir
	local src_files
	
	# Declare and initialize global varialbes
	fast_flag=""
	quiet_flag=""
	
	# Process our command line options.
	while getopts ":fq" opt; do
		case $opt in
			f  ) fast_flag="1" ;;
			q  ) quiet_flag="1" ;;
			\? ) print_usage
		esac
	done
	shift $(($OPTIND - 1))
	
	# Check the number of arguments after command line option processing.
	if [ $# != "2" ]; then
		print_usage
	fi
	
	# Extract arguments.
	script_name=${0##*/}
	the_routine="$1"
	src_dir="$2"
	
	# Check existence of the directory provided through the command-line.
	if [ ! -d "${src_dir}" ]; then
		echo "${script_name}: Source directory does not exist (${src_dir})."
		exit 1
	fi

	# Echo the current search if verbosity was requested.
	if [ -z "${quiet_flag}" ]; then
		echo "Searching for dependents of:"
	fi
	
	result_files=""	

	for item in ${the_routine}; do
	
		# Translate the routine name to uppercase.
		item=$(echo ${item} | tr '[:lower:]' '[:upper:]')
	
		# Get a list of the files in the source directory.
		src_files="$(ls ${src_dir}/*.f)"
	
		# Call a recursive routine that will find routine dependents until
		# no more are found.
		result_list=""	
		find_routine ${item} "${src_dir}" ""
	
		# Convert the list of routines into a list of their corresponding files.
		for rout in ${result_list}; do
		
			# Translate to lower case, and then add the fortran source file suffix.
			filename=$(echo ${rout} | tr '[:upper:]' '[:lower:]')
			filename="${filename}.f"
			result_files="${result_files} ${filename}"
		done
	
	done
	
	# Output the final list
	echo " "
	echo "Dependents of ${the_routine}:"
	echo "${result_files}"
	
	# Exit peacefully.
	return 0
}

find_routine()
{
	local routine
	local src_dir
	local tabs
	local newtabs
	
	local routine_files
	local routine_list
	local routine_set
	local temp_filename
	local temp_routine
	local prev_rout
	local rout
	
	# Extract arguments.
	routine="$1"
	src_dir="$2"
	tabs="$3"
	
	# Echo the current search if verbosity was requested.
	if [ -z "${quiet_flag}" ]; then
		if [ "${tabs}" == "" ]; then
			#echo -e "${tabs}${routine}"
			echo -e "${routine}"
		else
			echo -e "${tabs}^-- ${routine}"
		fi
	fi

	if [ "${tabs}" == "" ]; then
		newtabs=" "
	else
		newtabs="    "
	fi
	
	# Get a list of all file paths of files that call our routine.
	routine_files=$(grep "${routine}(" ${src_dir}/*.f | grep "CALL" | cut -d':' -f1)
	
	# Initialize the list of routines found to empty.
	routine_list=""
	
	# Process each file that contained the routine invocation.
	for file in ${routine_files}; do
		
		# Strip off the path, then strip the suffix of each file,
		# translate the result to uppercase, and then sort.
		temp_filename=${file##*/}
		temp_routine=${temp_filename%*.f}
		temp_routine=$(echo ${temp_routine} | tr '[:lower:]' '[:upper:]' | sort)
		routine_list="${routine_list} ${temp_routine}"
	done
	
	# Initialize accumulation variables to empty
	routine_set=""
	prev_rout=""
	
	# Add only the first occurance of each routine to our running list.
	for rout in ${routine_list}; do
		
		# Add the item to the list if the previous item was different
		# (ie: built a set).
		if [ "${prev_rout}" != "${rout}" ]; then
			
			# Given the search for routine xyz, exclude the hit we got from xyz.f
			# in our routine set, which represents the set of routines dependent
			# on ${routine} given.
			if [ "${rout}" != "${routine}" ]; then
				routine_set="${routine_set} ${rout}"
			fi
			
			# If it's not already a member of our global result list, add it.
			if ! is_in_list ${rout} "${result_list}"; then
				
				result_list="${result_list} ${rout}"
				
				# If fast flag was given on command line, we call find_routine()
				# here so we only search for routines that have not yet been found.
				if [ -n "${fast_flag}" ] ; then
					find_routine ${rout} ${src_dir} "${newtabs}${tabs}"
				fi
			fi
			
			# Update the previous routine.
			prev_rout=${rout}
		fi
		
	done
	
	# If fast flag was not given on command line, then we always call find_routine()
	# for the sake of building the hierarchy of dependents, even if we've already
	# searched for a routine before.
	if [ -z "${fast_flag}" ] ; then
		for rout in ${routine_set}; do
			
			# Recursively find dependents
			find_routine ${rout} ${src_dir} "${newtabs}${tabs}"
		done
	fi
	
	# Return peacefully.
	return 0
}


is_in_list()
{
	local rout
	local list
	local item
	
	# Extract arguments.
	rout="$1"
	list="$2"
	
	# Check all items in the list.
	for item in ${list}; do
		
		# Return "OK" status code if we find routine in list.
		if [ "${item}" = "${rout}" ]; then
			return 0
		fi
	done
	
	# Return "not-so-OK" status code if routine was never found.
	return 1
}


main "$@"
