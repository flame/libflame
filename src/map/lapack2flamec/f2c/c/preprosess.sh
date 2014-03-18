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
    echo "   The script cleans up fortran files to make them suitable for the f2c conversion." 
    echo " "

    declare -a words=( \
        'CHARACTER(1)'                                        'CHARACTER   ' \
        'RECURSIVE SUBROUTINE'                                'SUBROUTINE ' \
        ' CYCLE'                                              ' CALL F90_CYCLE' \
        ' EXIT'                                               ' CALL F90_EXIT' \
        ', TRANSFER'                                          ' '  \
        'TRANSFER (RWORK(1:2\\*N), (\/ (ZERO, ZERO) \/), N)'  'RWORK(1)' \
        ', MAXLOC'                                            ' ' \
        'MAXLOC( WORK( 1:N ), 1 )'                            'MAXLOC( WORK( 1 ) , N )' \
        'MAXLOC( WORK( (N+J):(2\\*N) ), 1 )'                  'MAXLOC( WORK( (N+J) ) , (N - J + 1) )' \
        )
    ## F90 intrinsic 
    ## F90_XXX are temporary and will be properly replaced later.
    ## TRANSFER is not necessary for C interfaces as it is replaced by a pointer.
    ## MAXLOC format needs to be changed to indicate the proper range of the vector.

##
##        'CHARACTER\\*(\\*)'                                   'CHARACTER' \
##        'INTRINSIC          LEN_TRIM'                         ' ' \
##        '1:LEN_TRIM( SRNAME )'                                '1' 
##        )

    src_dir="../netlib"
    one_tab="\t"

    # check src directory
    if [ ! -d "${src_dir}" ]; then
        echo "${script_name}: Source directory does not exist (${src_dir})."
        exit 1
    fi
    
    # create a sed stack
    files="$(find ${src_dir} -maxdepth 1 -name "*.f")"
    stack="sed \"s/${words[0]}/${words[1]}/g\""
    for (( counter=2; counter<${#words[@]}; counter+=2 )); do
        stack="${stack} | sed \"s/${words[${counter}]}/${words[${counter}+1]}/g\""
    done
    echo " Filter: "
    echo ${stack}
    echo " " 

    # execute the stacked filter
    for file in ${files}; do
        echo -ne "   Replacing ... ${file}                    "\\r
        tmp_file=$(echo "${file}.back")
        (   /bin/cp -f ${file} ${tmp_file} ; 
            eval "cat ${file} | ${stack}"  > ${tmp_file} && /bin/mv ${tmp_file} "${file}" ; 
            /bin/rm -f ${tmp_file} ) 
    done

    return 0
}

main "$@"
