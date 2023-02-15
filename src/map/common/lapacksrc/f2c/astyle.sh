#!/bin/bash


main()
{
    script_name=${0##*/}
    echo " "
    echo " "
    echo " "${script_name}
    echo " "
    echo " Objective:"
    echo "   The script reformats  files."
    echo " "
    
    ast=$(which astyle)
    
    if [ -f "$ast" ]
    then
        astyle --style=allman *.c *.h    
        rm -f *.orig
    else
        echo "$ast is not found."
    fi
}

main "$@"
