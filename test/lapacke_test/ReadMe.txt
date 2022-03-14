
Procedure to build Gtest based lapacke test-suite and Run the test-suite:
========================================================================

1) Install Gtest (Google test) by below Commands:
   Ex: cd /home/
  •	wget https://github.com/google/googletest/archive/refs/tags/release-1.11.0.tar.gz
  •	tar -xf release-1.11.0.tar.gz
  •	cd googletest-release-1.11.0/
  •	mkdir build;
  •	cd build;
  •	cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<Gtest_Install_folder_path>
  •	make -j
  •	make install


2) Build Netlib LAPACK (latest version 3.10.0) from the sources to generate (libblas.so, liblapacke.so).
  •	wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.0.tar.gz
  •	tar -xf v3.10.0.tar.gz
  •	cd lapack-3.10.0/
  •	mkdir build
  •	cd build
  •	cmake -DCMAKE_INSTALL_PREFIX=../LAPACK3p10-Install -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=ON -DLAPACKE=ON ../
  •	make -j
  •	make install

   The 'Netlib refrence libraries' (libblas.so, liblapacke.so, liblapack.so) are generated in the "<lapack-3.10.0>/LAPACK3p10-Install/lib" folder.

3) Edit the "make.inc" file with the following:
  •	Assign the variables 'BLASLIB', 'LIBFLAME_INSTALL_DIR' with the aocl-blis, aocl-libflame install directory paths respectively, to link the gtest application with AOCL libraries.
  •	Assign the variable 'GTESTDIR' with the <Gtest_Install_folder_path> (as mentioned in step-1).
  • Assign the variable 'NETLIB_INSTALL_DIR' to the folder path of Netlib reference libraries (as mentioned in step-2).

  NOTE: Default build settings is for 'gcc'.

  LLVM compiler Build settings: Below mofications needed in the 'make.inc' file.
  •CC = clang
  •CXX = clang++
  •FC = flang


5) Steps to run specific API tests:
  Edit the 'GT_LPKE.inc' file to include only the corresponding test case source files from line-2 onwards.

  EX: To include 3 test cases 'getrf', 'potrf'  API, the 'GT_LPKE.inc' file contents should be as below:
  -------------------------------------------------------------
  GTEST_SRC_FILE_LIST = LIN/lapacke_gtest_lin_common.o \
  LIN/lapacke_gtest_getrf.o \
  EIG/ORTHO/lapacke_gtest_geqrf.o \
  LIN/lapacke_gtest_potrf.o \
  -------------------------------------------------------------

5) In the terminal, Go to the '<lapacke_test>/' folder and run "make clean; make -j" command to generate the lapacke gtest executable "xlapacke_test_main" in the '<lapacke_test>/' folder.

NOTE: To enable Lapacke test suite automation, the below option is also available to override the 'Netlib refrence libraries' path variable.
make -j  "EXTFLAGS= -DREF_LPKE_LIB="<Netlib_lapacke_library_folder_path>/liblapacke.so.3.10.0" -DREF_BLAS_LIB="<Netlib_blas_library_folder_path>/libblas.so.3.10.0""

6) run ./xlapacke_test_main to run the lapacke test suite.
