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
    echo "   The script change FLA_ENABLE_LAPACK2FLAME -> FLA_ENABLE_NOT_YET_LAPACK2FLAME."
    echo " "

    one_tab="\t"

    # create sed stack
    files="FLA_gebrd.c FLA_gelqf.c FLA_gelsd.c FLA_geqpf.c FLA_geqrf.c FLA_gesdd.c FLA_gesvd.c FLA_getrf.c FLA_hegst.c FLA_hetrd.c FLA_lapack2flame_util.c FLA_lauum.c FLA_orgbr.c FLA_orglq.c FLA_orgqr.c FLA_orgtr.c FLA_ormbr.c FLA_ormlq.c FLA_ormqr.c FLA_ormtr.c FLA_potrf.c FLA_potri.c FLA_trtri.c"

    # execute stacked filter
    for file in ${files}; do
        echo -ne "   Replacing ... ${file}              "\\r
        tmp_file=$(echo "${file}.back")
        cp -f ${file} ${tmp_file}
        sed "s/FLA_ENABLE_LAPACK2FLAME/FLA_ENABLE_NOT_YET_LAPACK2FLAME/g" ${file} > ${tmp_file}
        # If you want to recover lapack2flame layer uncomment the following and comment the above.
        # sed "s/FLA_ENABLE_NOT_YET_LAPACK2FLAME/FLA_ENABLE_LAPACK2FLAME/g" ${file} > ${tmp_file}
        mv ${tmp_file} ${file}
        rm -f ${tmp_file}
    done

    return 0
}

main "$@"
