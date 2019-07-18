/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "omp.h"
#include "FLAME.h"
#include "LU_prototypes.h"

#define ENABLE_OMP_THREADS 1

int main(int argc, char *argv[])
{
  int 
    n, nfirst, nlast, ninc, i, irep,
    nrepeats, nb_alg;

  double
    dtime,
    dtime_best,
    gflops,
    max_gflops,
    diff,
    d_n;

  FLA_Obj
    A[3], Aref, Aold, delta;
  
  /* Initialize FLAME */
  FLA_Init( );

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf( "%% number of repeats:" );
  //scanf( "%d", &nrepeats );
  nrepeats = 3;
  printf( "%% %d\n", nrepeats );

  /* Enter the max GFLOPS attainable 
     This is used to set the y-axis range for the graphs. Here is how
     you figure out what to enter (on Linux machines):
     1) more /proc/cpuinfo   (this lists the contents of this file).
     2) read through this and figure out the clock rate of the machine (in GHz).
     3) Find out (from an expert of from the web) the number of floating point
        instructions that can be performed per core per clock cycle.
     4) Figure out if you are using "multithreaded BLAS" which automatically
        parallelize calls to the Basic Linear Algebra Subprograms.  If so,
        check how many cores are available.
     5) Multiply 2) x 3) x 4) and enter this in response to the below.
  */

  printf( "%% enter max GFLOPS:" );
  //scanf( "%lf", &max_gflops );
  max_gflops = 6.8;
  printf( "%% %lf\n", max_gflops );

  /* Enter the algorithmic block size */
  printf( "%% enter nb_alg:" );
  //scanf( "%d", &nb_alg );
  nb_alg = 128;
  printf( "%% %d\n", nb_alg );

  /* Turn on parameter checking */
  FLA_Check_error_level_set( FLA_FULL_ERROR_CHECKING );

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed */
  printf( "%% enter nfirst, nlast, ninc:" );
  //scanf( "%d%d%d", &nfirst, &nlast, &ninc );
  nfirst = 50; nlast = 500; ninc = 50;
  printf( "%% %d %d %d\n", nfirst, nlast, ninc );

  i = 1;
  for ( n=nfirst; n<= nlast; n+=ninc ){
   
    /* Allocate space for the matrices */
    FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &A[0] );

    FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &A[1] );

    FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &A[2] );
    FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &Aref );
    FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &Aold );
    FLA_Obj_create( FLA_DOUBLE, 1, 1, 0, 0, &delta );

    /* Generate random matrix A and save in Aold */
    FLA_Random_matrix( Aold );

    /* Add something large to the diagonal to make sure it isn't ill-conditionsed */
    d_n = ( double ) n;
    *( ( double * ) FLA_Obj_base_buffer( delta ) ) = d_n;
    FLA_Shift_diag( FLA_NO_CONJUGATE, delta, Aold );
    
    /* Set gflops = billions of floating point operations that will be performed */
    gflops = 2.0/3.0 * n * n * n * 1.0e-09;

    /* Time the reference implementation */

    for ( irep=0; irep<nrepeats; irep++ ){
      FLA_Copy( Aold, Aref );

      dtime = FLA_Clock();

      REF_LU( Aref );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    printf( "data_REF( %d, 1:2 ) = [ %d %le ];\n", i, n,
            gflops / dtime_best );
    fflush( stdout );

  FLA_Finalize( );

    /* Time FLA_LU */
#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for ( irep=0; irep<nrepeats; irep++ ){
  FLA_Init( );
  /* Turn on parameter checking */
  FLA_Check_error_level_set( FLA_FULL_ERROR_CHECKING );
      FLA_Copy( Aold, A[irep] );

      dtime = FLA_Clock();

      FLA_LU_nopiv( A[irep] );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );

  FLA_Finalize( );
    }

 //   sleep(5);
    printf("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n");
//    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A[2] );

//    int ij;
//    for(ij = 0; ij < n; ij++)
//	printf("%d %d %lf %lf %lf %lf\n", n, ij, buff_A[4 * ij + 0], buff_A[4 * ij + 1], buff_A[4 * ij + 2], buff_A[4 * ij + 3]);

    printf( "data_FLAME( %d, 1:2 ) = [ %d %le ];\n", i, n,
            gflops / dtime_best );
    fflush( stdout );


    /* Time the your implementations */


    /* Variant 1 unblocked */

#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for ( irep=0; irep<nrepeats; irep++ ){

  FLA_Init( );
      FLA_Copy( Aold, A[irep] );
    
      dtime = FLA_Clock();

      LU_unb_var1( A[irep] );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
  FLA_Finalize( );
    }    

    diff = FLA_Max_elemwise_diff( A[0], Aref );
    diff += FLA_Max_elemwise_diff( A[1], Aref );
    diff += FLA_Max_elemwise_diff( A[2], Aref );
    diff = diff / 3.0;

    printf( "data_unb_var1( %d, 1:3 ) = [ %d %le  %le];\n", i, n,
            gflops / dtime_best, diff );
    fflush( stdout );

    /* Variant 1 blocked */

#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for ( irep=0; irep<nrepeats; irep++ ){
  FLA_Init( );
      FLA_Copy( Aold, A[irep] );
    
      dtime = FLA_Clock();

      LU_blk_var1( A[irep], nb_alg );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
  FLA_Finalize( );
    }

    diff = FLA_Max_elemwise_diff( A[0], Aref );
    diff += FLA_Max_elemwise_diff( A[1], Aref );
    diff += FLA_Max_elemwise_diff( A[2], Aref );
    diff = diff / 3.0;

    printf( "data_blk_var1( %d, 1:3 ) = [ %d %le  %le];\n", i, n,
            gflops / dtime_best, diff );
    fflush( stdout );


    /* Variant 2 unblocked */

#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for ( irep=0; irep<nrepeats; irep++ ){

  FLA_Init( );
      FLA_Copy( Aold, A[irep] );
    
      dtime = FLA_Clock();

      LU_unb_var2( A[irep] );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
  FLA_Finalize( );
    }    

    diff = FLA_Max_elemwise_diff( A[0], Aref );
    diff += FLA_Max_elemwise_diff( A[1], Aref );
    diff += FLA_Max_elemwise_diff( A[2], Aref );
    diff = diff / 3.0;

    printf( "data_unb_var2( %d, 1:3 ) = [ %d %le  %le];\n", i, n,
            gflops / dtime_best, diff );
    fflush( stdout );

    /* Variant 2 blocked */

#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for ( irep=0; irep<nrepeats; irep++ ){
  FLA_Init( );
      FLA_Copy( Aold, A[irep] );
    
      dtime = FLA_Clock();

      LU_blk_var2( A[irep], nb_alg );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
  FLA_Finalize( );
    }

    diff = FLA_Max_elemwise_diff( A[0], Aref );
    diff += FLA_Max_elemwise_diff( A[1], Aref );
    diff += FLA_Max_elemwise_diff( A[2], Aref );
    diff = diff / 3.0;

    printf( "data_blk_var2( %d, 1:3 ) = [ %d %le  %le];\n", i, n,
            gflops / dtime_best, diff );
    fflush( stdout );

    /* Variant 3 unblocked */

#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for ( irep=0; irep<nrepeats; irep++ ){

  FLA_Init( );
      FLA_Copy( Aold, A[irep] );
    
      dtime = FLA_Clock();

      LU_unb_var3( A[irep] );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
  FLA_Finalize( );
    }    

    diff = FLA_Max_elemwise_diff( A[0], Aref );
    diff += FLA_Max_elemwise_diff( A[1], Aref );
    diff += FLA_Max_elemwise_diff( A[2], Aref );
    diff = diff / 3.0;

    printf( "data_unb_var3( %d, 1:3 ) = [ %d %le  %le];\n", i, n,
            gflops / dtime_best, diff );
    fflush( stdout );

    /* Variant 3 blocked */

#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for ( irep=0; irep<nrepeats; irep++ ){
  FLA_Init( );
      FLA_Copy( Aold, A[irep] );
    
      dtime = FLA_Clock();

      LU_blk_var3( A[irep], nb_alg );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
  FLA_Finalize( );
    }

    diff = FLA_Max_elemwise_diff( A[0], Aref );
    diff += FLA_Max_elemwise_diff( A[1], Aref );
    diff += FLA_Max_elemwise_diff( A[2], Aref );
    diff = diff / 3.0;

    printf( "data_blk_var3( %d, 1:3 ) = [ %d %le  %le];\n", i, n,
            gflops / dtime_best, diff );
    fflush( stdout );

    /* Variant 4 unblocked */

#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for ( irep=0; irep<nrepeats; irep++ ){

  FLA_Init( );
      FLA_Copy( Aold, A[irep] );
    
      dtime = FLA_Clock();

      LU_unb_var4( A[irep] );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
  FLA_Finalize( );
    }    

    diff = FLA_Max_elemwise_diff( A[0], Aref );
    diff += FLA_Max_elemwise_diff( A[1], Aref );
    diff += FLA_Max_elemwise_diff( A[2], Aref );
    diff = diff / 3.0;

    printf( "data_unb_var4( %d, 1:3 ) = [ %d %le  %le];\n", i, n,
            gflops / dtime_best, diff );
    fflush( stdout );

    /* Variant 4 blocked */

#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for ( irep=0; irep<nrepeats; irep++ ){
  FLA_Init( );
      FLA_Copy( Aold, A[irep] );
    
      dtime = FLA_Clock();

      LU_blk_var4( A[irep], nb_alg );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
  FLA_Finalize( );
    }

    diff = FLA_Max_elemwise_diff( A[0], Aref );
    diff += FLA_Max_elemwise_diff( A[1], Aref );
    diff += FLA_Max_elemwise_diff( A[2], Aref );
    diff = diff / 3.0;

    printf( "data_blk_var4( %d, 1:3 ) = [ %d %le  %le];\n", i, n,
            gflops / dtime_best, diff );
    fflush( stdout );

    /* Variant 5 unblocked */

#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for ( irep=0; irep<nrepeats; irep++ ){

  FLA_Init( );
      FLA_Copy( Aold, A[irep] );
    
      dtime = FLA_Clock();

      LU_unb_var5( A[irep] );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
  FLA_Finalize( );
    }    

    diff = FLA_Max_elemwise_diff( A[0], Aref );
    diff += FLA_Max_elemwise_diff( A[1], Aref );
    diff += FLA_Max_elemwise_diff( A[2], Aref );
    diff = diff / 3.0;

    printf( "data_unb_var5( %d, 1:3 ) = [ %d %le  %le];\n", i, n,
            gflops / dtime_best, diff );
    fflush( stdout );

    /* Variant 5 blocked */

#if ENABLE_OMP_THREADS
#pragma omp parallel
#pragma omp for
#endif
    for ( irep=0; irep<nrepeats; irep++ ){
  FLA_Init( );
      FLA_Copy( Aold, A[irep] );
    
      dtime = FLA_Clock();

      LU_blk_var5( A[irep], nb_alg );

      dtime = FLA_Clock() - dtime;

      if ( irep == 0 ) 
        dtime_best = dtime;
      else
        dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
  FLA_Finalize( );
    }

    diff = FLA_Max_elemwise_diff( A[0], Aref );
    diff += FLA_Max_elemwise_diff( A[1], Aref );
    diff += FLA_Max_elemwise_diff( A[2], Aref );
    diff = diff / 3.0;

    printf( "data_blk_var5( %d, 1:3 ) = [ %d %le  %le];\n", i, n,
            gflops / dtime_best, diff );
    fflush( stdout );

    FLA_Obj_free( &A[0] );

    FLA_Obj_free( &A[1] );

    FLA_Obj_free( &A[2] );
    FLA_Obj_free( &Aold );
    FLA_Obj_free( &Aref );
    FLA_Obj_free( &delta );

    printf( "\n" );

    i++;
  }

  /* Print the MATLAB commands to plot the data */

  /* Delete all existing figures */
  printf( "close all\n" );

  /* Plot the performance of FLAME */
  printf( "plot( data_FLAME( :,1 ), data_FLAME( :, 2 ), 'k--' ); \n" );

  /* Indicate that you want to add to the existing plot */
  printf( "hold on\n" );

  /* Plot the performance of the reference implementation */
  printf( "plot( data_REF( :,1 ), data_REF( :, 2 ), 'k-' ); \n" );

  /* Plot the performance of your implementations */

  printf( "plot( data_unb_var1( :,1 ), data_unb_var1( :, 2 ), 'r-.' ); \n" );
  printf( "plot( data_unb_var2( :,1 ), data_unb_var2( :, 2 ), 'g-.' ); \n" );
  printf( "plot( data_unb_var3( :,1 ), data_unb_var3( :, 2 ), 'b-.' ); \n" );
  printf( "plot( data_unb_var4( :,1 ), data_unb_var4( :, 2 ), 'm-.' ); \n" );
  printf( "plot( data_unb_var5( :,1 ), data_unb_var5( :, 2 ), 'c-.' ); \n" );
  printf( "plot( data_blk_var1( :,1 ), data_blk_var1( :, 2 ), 'r--' ); \n" );
  printf( "plot( data_blk_var2( :,1 ), data_blk_var2( :, 2 ), 'g--' ); \n" );
  printf( "plot( data_blk_var3( :,1 ), data_blk_var3( :, 2 ), 'b--' ); \n" );
  printf( "plot( data_blk_var4( :,1 ), data_blk_var4( :, 2 ), 'm--' ); \n" );
  printf( "plot( data_blk_var5( :,1 ), data_blk_var5( :, 2 ), 'c--' ); \n" );

  printf( "hold on \n");

  printf( "xlabel( 'matrix dimension m=n' );\n");
  printf( "ylabel( 'GFLOPS/sec.' );\n");
  printf( "axis( [ 0 %d 0 %3.1f ] ); \n", nlast, max_gflops );
  printf( "legend( 'FLA LU nopiv', ...\n");
  printf( "        'Simple loops', ...\n");
  printf( "        'unb var1', ...\n");
  printf( "        'unb var2', ...\n");
  printf( "        'unb var3', ...\n");
  printf( "        'unb var4', ...\n");
  printf( "        'unb var5', ...\n");
  printf( "        'blk var1', ...\n");
  printf( "        'blk var2', ...\n");
  printf( "        'blk var3', ...\n");
  printf( "        'blk var4', ...\n");
  printf( "        'blk var5', 2);\n");

  printf( "print -r100 -depsc LU.eps\n");

  FLA_Finalize( );

  return 0;
}
