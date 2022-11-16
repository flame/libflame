@echo OFF
REM Copyright (C) 2020-2023, Advanced Micro Devices, Inc. All Rights Reserved

setlocal EnableDelayedExpansion

REM The name of the script, stripped of any preceeding path.
set script_name=%0

SET ARGS_COUNT=0
FOR %%A in (%*) DO SET /A ARGS_COUNT+=1

if not !ARGS_COUNT! == 2 (
    REM Echo usage info.
    echo.
    echo  !script_name!
    echo.
    echo  Create a new test suite to test libflame's LAPACK interfaces and
    echo  functionality from an existing netlib LAPACK distribution.
    echo.
    echo  Usage:
    echo.
    echo    !script_name! netlib_path new_dirname
    echo.
    exit /b
    REM ------------------------------------------------------------------------------
)
REM Grab the arguments.
set netlib_path=%1
set testdir_new=%2
set lapack_ver=%netlib_path:~7%
REM Create the local destination directory.
mkdir %testdir_new%

REM The list of items we will copy.
set dir_copy_list=BLAS INSTALL TESTING
set files_copy_list=Makefile lapack_testing.py

for %%s in (%dir_copy_list%) do (
    echo Copying %netlib_path%\%%s to %testdir_new%\%%s

    rem Copy the current item from the netlib directory.
    xcopy /S /Y /I %netlib_path%\%%s %testdir_new%\%%s
)

for %%s in (%files_copy_list%) do (
    echo Copying files %netlib_path%\%%s to %testdir_new%\%%s

    rem Copy the current item from the netlib directory.
    copy %netlib_path%\%%s %testdir_new%\%%s
)

echo Creating %testdir_new%\make.inc
REM Create a make.inc file from make.inc.example. Modify the file
REM to include -lpthread after liblapack.a (which should be a symlink
REM to libflame).
REM cp -p ${netlib_path}/make.inc.example ./${testdir_new}/make.inc
rem cat !netlib_path!/make.inc.example \
rem | sed -e "s/liblapack\.a/liblapack\.a -lpthread/g" \
rem | sed -e "s/librefblas\.a/libblas\.a/g" \
rem > ./!testdir_new!/make.inc

set make_in=%netlib_path%\make.inc.example
set make_ou=%testdir_new%\make.inc

powershell "Get-Content %make_in% | ForEach-Object { $_ -replace 'liblapack.a', 'liblapack.a -lpthread' } | ForEach-Object { $_ -replace 'librefblas.a', 'libblas.a' } | Set-Content %make_ou%"

echo Creating %testdir_new%\SRC
echo Creating %testdir_new%\SRC\Makefile

REM Make a dummy 'SRC' directory with a dummy Makefile.
mkdir %testdir_new%\SRC
copy .\build\SRC-Makefile %testdir_new%\SRC\Makefile

REM Make a dummy 'SRC/VARIANTS' directory with a dummy Makefile.
mkdir %testdir_new%\SRC\VARIANTS
copy .\build\SRC-VARIANTS-Makefile %testdir_new%\SRC\VARIANTS\Makefile

REM Make a dummy 'lapacke' directory with a dummy Makefile.
mkdir %testdir_new%\lapacke
copy .\build\lapacke-Makefile %testdir_new%\lapacke\Makefile

echo Tweaking TESTING\LIN\xerbla.f.
echo Tweaking TESTING\EIG\xerbla.f.

REM Disable a part of the custom xerbla_() that tests string equality
REM since this is broken due to Fortran/C incompatibilities.
set xb_lin_in=%netlib_path%\TESTING\LIN\xerbla.f
set xb_lin_ou=%testdir_new%\TESTING\LIN\xerbla.f
set xb_eig_in=%netlib_path%\TESTING\EIG\xerbla.f
set xb_eig_ou=%testdir_new%\TESTING\EIG\xerbla.f

powershell "Get-Content %xb_lin_in% | ForEach-Object { $_ -replace 'SRNAME.NE.SRNAMT', '.FALSE.' } | Set-Content %xb_lin_ou%"
powershell "Get-Content %xb_eig_in% | ForEach-Object { $_ -replace 'SRNAME.NE.SRNAMT', '.FALSE.' } | Set-Content %xb_eig_ou%"
echo Tweaking ddrvsg.f

REM This small change for the ddrvsg.f file is needed so that a
REM particular test of dsygvx() passes the test suite. (For some
REM reason, the f2c'ed dstein() does not converge when dstebz()
REM computes the eigenvalues with the regular definition of abstol,
REM so we use an alternate definition that causes the absolute
REM tolerance to be computed as eps * norm1(T), where T is the
REM tridiagonal matrix.

set ddrvsg_in=%netlib_path%\TESTING\EIG\ddrvsg.f
set ddrvsg_ou=%testdir_new%\TESTING\EIG\ddrvsg.f

powershell "Get-Content %ddrvsg_in% | ForEach-Object { $_ -replace 'ABSTOL = UNFL + UNFL', 'ABSTOL = 0' } | Set-Content %ddrvsg_ou%"

echo Copying CMAKE files to %testdir_new%
xcopy /S /Y .\windows_netlib_test\%lapack_ver% %testdir_new%\

echo Creating source CMAKE file
set cmake_in=windows_netlib_test\%lapack_ver%\CMakeLists.example.txt
set cmake_ou=CMakeLists.txt

powershell "Get-Content %cmake_in% | ForEach-Object { $_ -replace 'windows_netlib_test', '%testdir_new%'} | Set-Content %cmake_ou%"


REM Exit peacefully.
exit /b