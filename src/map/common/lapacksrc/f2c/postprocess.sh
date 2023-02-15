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
    echo "   The script polishes f2c'ed files." 
    echo " "

    # sed filter stack
    declare -a words=( \
        'f90_cycle__()'             'continue' \
        'f90_exit__()'              'break' \
        'include \"f2c.h\"'         'include \"FLA_f2c.h\"' \
        'dotc_'                     'dotc_f2c_' \
        'dotu_'                     'dotu_f2c_' \
        'abs('                      'f2c_abs' \
        '__('                       '_(' \
        )

    # '\\#include \"blaswrap.h\"' ' ' \
    # 's_cat'                     's_kat' \
    #

    src_dir="."
    one_tab="\t"

    # check src file
    if [ ! -d "${src_dir}" ]; then
        echo "${script_name}: Source directory does not exist (${src_dir})."
        exit 1
    fi

    # create sed stack
    files="$(find ${src_dir} -maxdepth 1 -name "*.c")"
    stack="sed \"s/${words[0]}/${words[1]}/g\""
    for (( counter=2; counter<${#words[@]}; counter+=2 )); do
        stack="${stack} | sed \"s/${words[${counter}]}/${words[${counter}+1]}/g\""
    done
    echo " Filter: "
    echo ${stack}
    echo " "
    
    # execute stacked filter
    for file in ${files}; do
        echo -ne "   Replacing ... ${file}              "\\r
        tmp_file=$(echo "${file}.back")

        (   cp -f ${file} ${tmp_file};
            eval "cat ${file} | ${stack}"  > ${tmp_file} ;
            mv ${tmp_file} ${file};
            rm -f ${tmp_file}; ) 
    done

    # apply exceptions
    #
    # for 3.4.2, xpstrf and xpstf2, f90 intrinsic maxloc should be replaced by provided smaxloc and dmaxloc
    files="$(find ${src_dir} -maxdepth 1 -name "[sc]pstrf.c" -o -name "[sc]pstf2.c")"
    for file in ${files}; do
        echo -ne "   Replacing ... ${file}              "\\r
        tmp_file=$(echo "${file}.back")
        (   cp -f ${file} ${tmp_file} ;
            eval "cat ${file} | sed \"s/maxloc\_/smaxloc\_/g\""  > ${tmp_file} ;
            mv ${tmp_file} "${file}" ;
            rm -f ${tmp_file} ) 
    done
    #
    files="$(find ${src_dir} -maxdepth 1 -name "[dz]pstrf.c" -o -name "[dz]pstf2.c")"
    for file in ${files}; do
        echo -ne "   Replacing ... ${file}              "\\r
        tmp_file=$(echo "${file}.back")
        (   cp -f ${file} ${tmp_file} ;
            eval "cat ${file} | sed \"s/maxloc\_/dmaxloc\_/g\""  > ${tmp_file} ;
            mv ${tmp_file} "${file}" ;
            rm -f ${tmp_file} ) 
    done
    # # change double to real strictly...
    # files="$(find ${src_dir} -maxdepth 1 -name "[sc]*.c")"
    # for file in ${files}; do
    #     echo -ne "   Replacing ... ${file}              "\\r
    #     tmp_file=$(echo "${file}.back")
    #     cp -f ${file} ${tmp_file}
    #     eval "cat ${file} | sed \"s/double\ /real\ /g\" | sed \"s/\ sqrt(/\ sqrtf(/g\" | sed \"s/\ sqrtf(doublereal)/\ sqrtf(real)/g\""  > ${tmp_file}
    #     mv ${tmp_file} "${file}"
    #     rm -f ${tmp_file}
    # done

    return 0
}

main "$@"
