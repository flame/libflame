#!/bin/bash

#
# dvi-ps.sh
#
# Field G. Van Zee
#
# Takes a LaTeX file as input and generates device-independent and/or 
# postscript files in a specified directory.
#
# At least ONE of the output options must be specified for the script to output
# anything. Otherwise, the script is executed by default in "dry-run" mode.
#
# Usage:
#   dvi-ps.sh [options] input_file.tex output_file_root
#
# The following options are accepted:
# 
#   -d            output device-independent (.dvi) file
#   -p            output postscript (.ps) file
#   -v            verbose
#

print_usage()
{
	local script_name
	
	# Get the script name and strip off its path.
	script_name=${0##*/}
	
	echo " "
	echo " $script_name"
	echo " "
	echo " Field G. Van Zee"
	echo " "
	echo " Takes a LaTeX file as input and generates device-independent and/or "
	echo " postscript files in a specified directory."
	echo " "
	echo " At least ONE of the output options must be specified for the script to output"
	echo " anything. Otherwise, the script is executed by default in \"dry-run\" mode."
	echo " "
	echo " Usage:"
	echo "   $script_name [options] input_file.tex output_file_root"
	echo " "
	echo " The following options are accepted:"
	echo " "
	echo "   -d            output device-independent (.dvi) file"
	echo "   -p            output postscript (.ps) file"
	echo "   -v            verbose"
	echo " "

	exit 1
}

gen_dvi_ps()
{
	local script_name	
	local input_file_path_tex
	local output_file_path_base
	local input_file_name_tex
	local output_dir
	local output_file_name_base
	local output_file_name_dvi
	local output_file_name_ps
	local output_file_name_aux
	local output_file_name_log
	
	
	# Extract arguments.
	script_name=$0
	input_file_path_tex=$1
	output_file_path_base=$2
	
	
	# Strip the path from the script name
	script_name="${script_name##*/}"
	
	
	# Determine input filename by stripping the path.
	input_file_name_tex="${input_file_path_tex##*/}"
	
	
	# Determine the output filename root by stripping the path.
	output_file_name_base="${output_file_path_base##*/}"
	
	
	# Determine the output directory by stripping the filename name root.
	output_dir="${output_file_path_base%/*}"
	
	
	# Determine the various output filenames by appending suffixes to the
	# base filename.
	output_file_name_aux="${output_file_name_base}.aux"
	output_file_name_log="${output_file_name_base}.log"
	output_file_name_dvi="${output_file_name_base}.dvi"
	output_file_name_ps="${output_file_name_base}.ps"
	
	
	# If some form of output was requested...
	if [ -n "$output_dvi_flag" ] || [ -n "$output_ps_flag" ]; then
		
		# Run latex; produce a device-independent file. Log
		# output to $script_name.log. Check for failure.
		if latex --interaction nonstopmode \
		      "${input_file_path_tex}" 1>> "${script_name}.log"
		then
			# If verboseness was requested, then echo .dvi file name.
			if [ -n "$verbose_flag" ]; then
				# Echo success.
				echo "Generated ${output_file_name_dvi}..."
			fi
		else
			# Echo an error, remove some temporary files, (but
			# leave the log file) end exit.
			echo "$script_name: latex returned error while processing ${input_file_path_tex}."
			rm -f "${output_file_name_dvi}"
			rm -f "${output_file_name_aux}"
			rm -f "texput.log"
			exit 1
		fi
	fi	
	
	
	# If the dvi file was requested, then move it to the output directory.
	# Otherwise, we can delete it from the current directory.
	if [ -n "$output_dvi_flag" ]; then
		
		# Move the dvi file to the output directory.
		if mv "${output_file_name_dvi}" \
		      "${output_dir}/${output_file_name_dvi}" 2>> "${script_name}.log"
		then
			# If verboseness was requested, then echo .dvi file name.
			if [ -n "$verbose_flag" ]; then
				echo "Moved to  ${output_dir}/${output_file_name_dvi}."
			fi
		else
			# Echo an error and exit.
			echo "$script_name: mv returned error when trying to move ${output_file_name_dvi} to ${output_dir} directory."
			exit 1
		fi
	else
		# Delete the dvi file.
		rm -f "$output_file_name_dvi"
	fi
	
	
	# If the ps file was requested, then we create it from the dvi file
	# in the output directory given.
	if [ -n "$output_ps_flag" ]; then
		
		# Convert the dvi file to a ps file.
		if dvips -f "${output_dir}/${output_file_name_dvi}" -t letter > "${output_dir}/${output_file_name_ps}" 2>> "${script_name}.log"
		then
			# If verboseness was requested, then echo .ps file name.
			if [ -n "$verbose_flag" ]; then
				echo "Generated ${output_dir}/${output_file_name_ps}."
			fi
		else
			# Echo an error, remove bad ps file, and exit.
			echo "$script_name: dvips returned error while processing ${output_dir}/${output_file_name_dvi}."
			rm -f "${output_dir}/${output_file_name_ps}"
			exit 1
		fi
	fi
	
		
	# Clean up temporary files.
	rm -f "${output_file_name_aux}" \
       	      "${output_file_name_log}" \
       	      "texput.log"
	
	
	# Return peacefully.
	return 0
}


main()
{
	
	# Initialize global getopts flags to null.
	output_dvi_flag=""
	output_ps_flag=""
	verbose_flag=""
	
	
	# Process our command line options.
	while getopts ":dpv" opt; do
		case $opt in
			d  ) output_dvi_flag="1" ;;
			p  ) output_ps_flag="1" ;;
			v  ) verbose_flag="1" ;;
			\? ) print_usage
		esac
	done
	shift $(($OPTIND - 1))
	
	
	# Check the number of arguments after command line option processing.
	if [ $# != "2" ]; then
		print_usage
	fi
	
	
	# Generate the dvi and/or ps files from the script's input latex file.
	gen_dvi_ps "$@"
	
	
	# Exit peacefully.
	exit 0
}


# The script's main entry point, passing all parameters given.
main "$@"
