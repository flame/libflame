#!/bin/bash

main()
{
    # local variables
    script_name=${0##*/}
    head_file="FLA_lapack2flame_prototypes.h"
    
    echo " "
    echo " "
    echo " "${script_name}
    echo " "
    echo " Objective:"
    echo "   The script generates a C header file collecting function names from the f2c'ed files."
    echo " "
    echo " Output: "${head_file}
    echo " "

    rm -f ${head_file}

    # echo "/* This header is automatically generated from f2c'ed files where ftnlen arguments are removed. */" > ${head_file}

    # dump all c files into a header
    files="$(find . -name "*.c")"
    for file in ${files}; do
        srname=$(basename ${file%.*})
        echo -ne "   Generating header ... ${file}            "\\r
        cat ${file} \
            | sed "s/${srname}[_]*(.*);//g" \
            | sed -n "s/\(${srname}[_]*(.*)\).*/\1;/p" \
            >> ${head_file}
    done

    return 0
}

main "$@"
