#!/bin/bash

#
# runme.sh
# Field G. Van Zee
#

print_usage()
{
	echo "Usage: $0 quoted_dir_list" 
	exit 1
}

main()
{
	local var_dirs
	local input_files
	local output_filepath
	local input_suffix
	local driver
	local output_prefix
	
	# Check the number of arguments
	if [ $# != "1" ]; then
		print_usage
	fi

	# Extract arguments.
	var_dirs="$1"
	
	# Establish some constants.
	driver=$(ls test*.x)
	output_prefix="output"

	
	for var_dir in ${var_dirs}; do
		
		input_files=$(ls ${var_dir}/input*)
		
		for input_filepath in ${input_files}; do
			
			input_filename=${input_filepath##*/}
			input_suffix=${input_filename##input_}
			output_filepath="${var_dir}/${output_prefix}_${input_suffix}.m"
			
			echo ">>> ${driver} < ${input_filepath} > ${output_filepath}"
			./${driver} < ${input_filepath} > ${output_filepath}
		
		done
		
	done
	
	# Exit peacefully.
	return 0
}

main "$@"
