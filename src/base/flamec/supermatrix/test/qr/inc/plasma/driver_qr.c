/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "plasma.h"
#include "FLAME.h" 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 


#define OUTPUT_PATH "./results"
#define OUTPUT_FILE "qr_plasma"


integer main( integer argc, char *argv[] ) 
{ 
   integer
      i, j,
      size,
      n_threads,
      n_repeats,
      n_trials,
      nb_alg,
      increment,
      begin;

   FLA_Datatype
      datatype = FLA_DOUBLE;

   FLA_Obj 
      A;
   
   double 
      b_norm_value = 0.0,
      dtime, 
      *dtimes,
      *flops,
      *T;

   char
      output_file_m[100];
   
   FILE
      *fpp;

   fprintf( stdout, "%c Enter number of repeats: ", '%' );
   scanf( "%d", &n_repeats );
   fprintf( stdout, "%c %d\n", '%', n_repeats );

   fprintf( stdout, "%c Enter blocksize: ", '%' );
   scanf( "%d", &nb_alg );
   fprintf( stdout, "%c %d\n", '%', nb_alg );

   fprintf( stdout, "%c Enter problem size parameters: first, inc, num: ", '%' );
   scanf( "%d%d%d", &begin, &increment, &n_trials );
   fprintf( stdout, "%c %d %d %d\n", '%', begin, increment, n_trials );

   fprintf( stdout, "%c Enter number of threads: ", '%' );
   scanf( "%d", &n_threads );
   fprintf( stdout, "%c %d\n\n", '%', n_threads );

   sprintf( output_file_m, "%s/%s_output.m", OUTPUT_PATH, OUTPUT_FILE );
   fpp = fopen( output_file_m, "a" );

   fprintf( fpp, "%%\n" );
   fprintf( fpp, "%% | Matrix Size |    PLASMA   |\n" );
   fprintf( fpp, "%% |    n x n    |    GFlops   |\n" );
   fprintf( fpp, "%% -----------------------------\n" );

   FLA_Init();
   PLASMA_Init( n_threads );

   PLASMA_Disable( PLASMA_AUTOTUNING );
   PLASMA_Set( PLASMA_TILE_SIZE, nb_alg );
   PLASMA_Set( PLASMA_INNER_BLOCK_SIZE, nb_alg / 4 );

   dtimes = ( double * ) FLA_malloc( n_repeats * sizeof( double ) );
   flops  = ( double * ) FLA_malloc( n_trials  * sizeof( double ) );

   fprintf( fpp, "%s = [\n", OUTPUT_FILE );

   for ( i = 0; i < n_trials; i++ )
   {
      size = begin + i * increment;
      
      FLA_Obj_create( datatype, size, size, 0, 0, &A );
      
      for ( j = 0; j < n_repeats; j++ )
      {
         FLA_Random_matrix( A );

         PLASMA_Alloc_Workspace_dgeqrf( size, size, &T );

         dtime = FLA_Clock();

         PLASMA_dgeqrf( size, size, FLA_Obj_buffer_at_view( A ), size, T );

         dtime = FLA_Clock() - dtime;
         dtimes[j] = dtime;
         
         free( T );
      }
      
      dtime = dtimes[0];
      for ( j = 1; j < n_repeats; j++ )
         dtime = fla_min( dtime, dtimes[j] );
      flops[i] = 4.0 / 3.0 * size * size * size / dtime / 1e9;
      
      fprintf( fpp, "   %d   %6.3f\n", size, flops[i] );

      printf( "Time: %e  |  GFlops: %6.3f\n",
              dtime, flops[i] );
      printf( "Matrix size: %d x %d  |  nb_alg: %d\n",
              size, size, nb_alg );
      printf( "Norm of difference: %le\n\n", b_norm_value );

      FLA_Obj_free( &A );
   }

   fprintf( fpp, "];\n" );
   
   fflush( fpp );
   fclose( fpp );
   
   FLA_free( dtimes );
   FLA_free( flops );

   PLASMA_Finalize();
   FLA_Finalize(); 
   
   return 0; 
}
