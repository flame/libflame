#!/bin/bash

main()
{
    # local variables
    script_name=${0##*/}

    echo " "
    echo " "
    echo " "${script_name}
    echo " "
    echo " Objective:"
    echo "   The script moves (with -f) routines in the ../netlib/support.list to other directory." 
    echo " "


    file="../netlib/support.list"
    tgt="other"
    one_tab="\t"

    # check the file
    if [ ! -f "${file}" ]; then
        echo "${script_name}: Source directory does not exist (${file})."
        exit 1
    fi
    if [ ! -d "${tgt}" ]; then
        echo "${script_name}: Target directory does not exist (${tgt})."
        exit 1
    fi
    
    # execute the stacked filter
    while read line; do
        cfile=${line%%.*}.c
        echo -ne "   Moving ... ${cfile} to ${tgt}                 "\\r
        mv -f ${cfile} ${tgt}
    done < ${file}

    return 0
}

main "$@"
