#!/bin/bash

usage() {
    echo "  "
    echo "  "
    echo " This script has a few uses. It should be run from the scripts/ dir "
    echo " Usage: ./regen-files.sh [Action, Required] [Path to lapack tar, optional] "
    echo "  "
    echo " Arg1: " 
    echo "      'Build' This action is used by the libFlame Makefile at build "
    echo "              time to extract the lapack fortran files "
    echo "      'build_test' This action is used by the libFlame Makefile at build "
    echo "              time to extract the lapack netlib-test for users "
    echo "      'other_installs' This action can be used by a dev to extract netlibs-test "
    echo "              and generate new f2c files in flablas/f2c "
    echo "      'cleanup' This action is used by the libFlame Makefile at build "
    echo "              to clean up any lapack fortran files that were used to build exeicutables "
    echo "      'clean_test' This action is used by the libFlame Makefile at build "
    echo "              to clean up netlibs-test that were made "
    echo "  "
    echo " Arg2: "
    echo "      '../netlib/lapack-X.X.X.tgz' This path should lead to the tar containing " 
    echo "              the current version of lapack libFlame has installed. If no tar is given "
    echo "              this script will look for one in ../netlibs "
    echo "              If this script finds more than 1 tar, this script will error out "
    echo "  "
    echo "  "
}

clean() {
    rm -f $script_dir/../netlib/._*
    rm -f $script_dir/../netlib/*.{f,F,f90,F90}
    rm -f $script_dir/../netlib/*~
    rm -rf $script_dir/../netlib/lapack-*/
    rm -f $script_dir/../fortran/*.{f,F,f90,F90}
    rm -f $script_dir/../fortran/depen_mods/*.{f,F,f90,F90}
}

## Pass an install type and a tar if needed
install_type=$1
tar_file=$2

## Get current script dir
script_dir=$(cd "$(dirname "$0")" && pwd)

if [[ "$install_type" == "" || "$install_type" == help ]]
then
    usage
    exit
fi

if [[ "$tar_file" == "" && "$install_type" != "cleanup" ]]
then
    ## Try to auto select a tar if the user doesn't provide one
    num_of_tars=$(ls -1q $script_dir/../netlib/*.{tar,tgz,tar.gz} 2> /dev/null | wc -l)

    ## Check to make sure there is only one tar file in netlib/
    ## If so, use it
    if [[ "$num_of_tars" == 1 ]]
    then
        tar_file=$(ls $script_dir/../netlib/ | grep ".tar\|.tgz\|.tar.gz")
        echo "No tar given, using $tar_file"
    else
        echo "Could not determine which tar to use."
        usage
        exit
    fi
fi

if [[ "$install_type" == build || "$install_type" == build_test || "$install_type" == other_installs ]]
then
    if [[ "$install_type" == build ]]
    then
        ## Make sure there are no old files
        clean

        ## Untar
        ## It was observed that the tars had only the src files
        ## While the tgz/tar.gz had the full netlibs. However, this
        ## is not a guarantee. When adding a new netlibs tar, the dev
        ## should verify this is still true
        if [[ "$tar_file" == *.tgz || "$tar_file" == *.tar.gz ]]
        then
            cd $script_dir/../netlib
            tar -xzf $tar_file
            cp lapack-*/SRC/*.{f,F,f90,F90} .
            cp lapack-*/INSTALL/*roundup_lwork.{f,F,f90,F90} .
            cp lapack-*/INSTALL/ilaver.{f,F,f90,F90} .
            rm -f ._*.f
            cat support.list | xargs rm -f
            cat install.list | xargs rm -f
        elif [[ "$tar_file" == *.tar ]]
        then
            cd $script_dir/../netlib
            tar -xf $tar_file
            rm -f ._*.f
            cat support.list | xargs rm -f
            cat install.list | xargs rm -f
        else
            echo "Please pass a netlibs tar file as arg 2"
            exit
        fi

        ## Move to the lapack dir
        mkdir $script_dir/../fortran/
        mv $script_dir/../netlib/*.{f,F,f90,F90} $script_dir/../fortran/

        ## Disable xblas-related code so shared libraries link properly.
        cd $script_dir
        ./disable_xblas.sh
    elif [[ "$install_type" == build_test ]]
    then
        ## This section lets a user update netlibs-test and  
        ## Updating other areas in libflame
        ## It was observed that the tars had only the src files
        ## While the tgz/tar.gz had the full netlibs. However, this
        ## is not a guarantee. When adding a new netlibs tar, the dev
        ## should verify this is still true
        if [[ "$tar_file" == *.tgz || "$tar_file" == *.tar.gz ]]
        then
            cd $script_dir/../netlib
            tar -xzf $tar_file

            ## Get the name of the directory
            new_test_dir=$(echo $tar_file | sed "s/\.tgz//" | sed "s/\.tar.gz//" | sed "s/\..\/netlib\///")

            ## Call the script that makes a new netlibs-test
            cd $script_dir/../../../../../netlib-test/
            ./create_new_testdir.sh $script_dir/../netlib/$new_test_dir $new_test_dir
        else
            echo " "
            echo " Skipping adding netlibs-test because it is not in the tar "
            echo " "
            exit
        fi
    elif [[ "$install_type" == other_installs ]]
    then
        ## This section lets a user update netlibs-test and  
        ## Updating other areas in libflame
        ## It was observed that the tars had only the src files
        ## While the tgz/tar.gz had the full netlibs. However, this
        ## is not a guarantee. When adding a new netlibs tar, the dev
        ## should verify this is still true
        if [[ "$tar_file" == *.tgz || "$tar_file" == *.tar.gz ]]
        then
            cd $script_dir/../netlib
            tar -xzf $tar_file

            ## Get the name of the directory
            new_test_dir=$(echo $tar_file | sed "s/\.tgz//" | sed "s/\.tar.gz//" | sed "s/\..\/netlib\///")

            ## Regen the f2c in src/flablas/f2c/
            cd $script_dir/../../../../flablas/f2c/
            rm -f *.{f,F,f90,F90} *.c
            cp $script_dir/../netlib/$new_test_dir/BLAS/SRC/* .
            ./regen-files.sh
        else
            echo " "
            echo " Skipping other installs because they are not in the tar "
            echo " "
            exit
        fi
    else
        echo "Unknown Error"
        exit
    fi

    ## Cleanup
    rm -f $script_dir/../netlib/._*
    rm -f $script_dir/../netlib/*.{f,F,f90,F90}
    rm -rf $script_dir/../netlib/lapack-*/
elif [[ "$install_type" == cleanup ]]
then
    ## Cleanup
    echo "Cleaning up any files from lapack tar."
    clean
elif [[ "$install_type" == clean_test ]]
then
    ## Only clean out the test we made from this script
    echo "Cleaning up any netlibs-tests."
    test_dir=$(echo $tar_file | sed "s/\.tgz//" | sed "s/\.tar.gz//" | sed "s/\..\/netlib\///")
    rm -rf $script_dir/../../../../../netlib-test/$test_dir
else
    echo "Invalid Arg."
    usage
    exit
fi
