#!/bin/bash

print_usage()
{
    local script_name

    # Get the script name
    script_name=${0##*/}

    # Echo usage info
    echo " "
    echo " "${script_name}
    echo " "
    echo " Usage:"
    echo "   ${script_name} old_word new_word src_dir extension"
    echo "  ${script_name} CHARACTER\(1\) CHARACTER netlib f"
    echo " "

    # Exit with non-zero exit status
    exit 1
}

main()
{
    # Check the number of argements
    if [ $# != "4" ]; then
        print_usage
    fi

    # local variables
    script_name=${0##*/}
    old_word="$1"
    new_word="$2"
    src_dir="$3"
    extension="$4"
    one_tab="\t"

    # Check src file
    if [ ! -d "${src_dir}" ]; then
        echo "${script_name}: Source directory does not exist (${src_dir})."
        exit 1
    fi

     files="$(find ${src_dir} -name "*.${extension}")"
#     echo "${files}"

     for file in ${files}; do
         echo "Replacing ... ${file}"
         tmp_file=$(echo "${file}.back")
         /bin/cp -f ${file} ${tmp_file}
         sed "s/${old_word}/${new_word}/g" "${file}" > ${tmp_file} && /bin/mv ${tmp_file} "${file}"
         /bin/rm -f ${tmp_file}
     done
     return 0
}

main "$@"
