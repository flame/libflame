
Procedure to build Google Bench based libflame benchmark test-suite and Run the test-suite:
=============================================================================================

1) Install Google Bench by compiling from the sources. Refer the below link for installation steps:
   https://www.intellect.ind.in/pages/installing-google-benchmark-library-on-linux-from-source-code.html

2) export LD_LIBRARY_PATH, PATH variables with the google bench library, install folder paths respectively.

3) Update the blis, libflame library paths variables 'BLASLIB', 'LAPACKLIB' variables in the 'make.inc' file.

4) The test suite has below 3 mode support:
   i)  Discrete mode (0): Runs the matrix sizes for given discrete sizes.
       Ex: Effective Parameters in the .dat file for matrix sizes:
            Mode 0  // Enter:- 0: discrete, 1: Combinational, 2: Range with step increments of Matrix sizes
            Num_Tests 2
            M     100 40
            N     60 30
        In this mode the tests run for the below (M,N) combinations:
        (100, 60) , (40,30)

   ii) Combinational mode (1): Runs the matrix sizes given for m, n combinations.
       Ex: Effective Parameters in the .dat file for matrix sizes:
            Mode 1  // Enter:- 0: discrete, 1: Combinational, 2: Range with step increments of Matrix sizes
            Num_Tests 2
            M     100 40
            N     60 30
        In this mode the tests run for the below (M,N) combinations:
        (100, 60) , (100, 30), (40, 60), (40,30)
       
   iii) Range mode (1): Runs the matrix sizes for given ranges & step sizes.
       Ex: Effective Parameters in the .dat file for matrix sizes:
            Mode 2  // Enter:- 0: discrete, 1: Combinational, 2: Range with step increments of Matrix sizes
            Num_Ranges 2
            MRange_start 10 100
            MRange_end   30 300
            MRange_step_size 10 100
            NRange_start 5 120
            NRange_end   20 300
            NRange_step_size 15 130
            
        In this mode the tests run for the below (M,N) combinations:
            Range-1: 
            M     10 20 30
            N     5  20
            the below (M,N) combinations run for range-1:
            (10, 5) , (10, 20), (20, 5), (20,20),  (30, 5), (30,20)
            
            Range-2: 
            M     100 200 300
            N     120 250
            the below (M,N) combinations run for range-1:
            (100, 120) , (100, 250), (200, 120) , (200, 250), (300, 120) , (300, 250)
            

4) Build Procedure:  
   
   i)  Provide Blis, libflame library paths (for the BLASLIB, LAPACKLIB variables) in the make.inc file.
   ii) In the terminal, Go to the '<libflame>/performance_test/' folder and run "make clean; make -j" command to 
       generate the libflame's test suite application "xFlame_gbench".
   iii) Execute "./xFlame_gbench" to run the test suite.


