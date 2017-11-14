#!/usr/bin/env bash

#
# print-tree.sh
#
# Field G. Van Zee
#
# Usage:
#   print-tree.sh [-d] dir1 [dir2 ...]
#
# Recursively descends through the directory arguments given and prints to
# standard output all items found (files and directories). Each item printed
# is indented with an additional '\t' tab character according to its depth
# in the directory hierarchy. 
# 


main()
{
	# Process our command line options. We only respond to the -d flag,
	# in which case we'll only echo directories as we go along.
	while getopts ":d" opt; do
		case $opt in
			d  ) dirs_only_flag="1" ;;
			\? ) echo 'Usage: print-tree.sh [-d] dir1 [dir2 ...]'
				exit 1
		esac
	done
	shift $(($OPTIND - 1))
	
	
	# Define our tab character
	one_tab="\t"
	
	
	# Process each command line argument given
	for top_thing in "$@"; do
		
		# Echo the item
		echo $top_thing
		
		
		# If top_thing is a directory, then descend recursively into it.
		if [ -d "$top_thing" ]; then
			this_thing=$top_thing
			print_tree $(command ls $top_thing)
		fi
	done
	
	
	# Exit peacefully.
	return 0
}

print_tree()
{
	# Add a tab character to our tab string
	tab_str=$tab_str$one_tab
	
	
	# Process each item in our argument list (ie: each item in this_thing)
	for thing in "$@"; do
		
		# Append the current item to the current directory
		this_thing=$this_thing/$thing
		
		
		# Echo the tabbed item. If the dirs_only_flag was given,
		# and $thing is a directory, then output. Otherwise, do not
		# output. If the dirs_only_flag was not given, then output
		# regardless of file type.
		if [ -n "$dirs_only_flag" ]; then
			if [ -d $this_thing ]; then
				echo -e "$tab_str$thing"
			fi
		else
			echo -e "$tab_str$thing"
		fi
		
		
		# If this_thing is a directory, then descent recursively into it.
		if [ -d "$this_thing" ]; then
			print_tree $(command ls $this_thing)
		fi
		
		
		# Delete the end of the path, up to the first / character to
		# prepare for the next "file" in $@
		this_thing=${this_thing%/*}
	done
	
	
	# Remove one tab from the end of our tab_str string, given that
	# we're about to return up the function stack one level.
	tab_str=${tab_str%"\t"}
}

# The script's main entry point, passing all parameters given.
main "$@"

