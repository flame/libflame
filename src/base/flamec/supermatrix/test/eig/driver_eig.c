/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h" 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 


#define OUTPUT_PATH "./results"
#define OUTPUT_FILE "eig"


integer main( integer argc, char *argv[] ) 
{ 
   integer
      i, j,
      n_threads,
      n_repeats,
      n_trials,
      increment,
      begin,
      sorting,
      caching,
      work_stealing,
      data_affinity;

   dim_t
      size,
      nb_alg;

   FLA_Datatype
      datatype = FLA_DOUBLE;

   FLA_Inv    
      inv = FLA_NO_INVERSE;

   FLA_Uplo
      uplo = FLA_LOWER_TRIANGULAR;

   FLA_Obj 
      A, B, x, b, b_norm,
      AH, BH;
   
   double 
      length,
      b_norm_value = 0.0,
      dtime, 
      *dtimes,
      *flops;

#ifndef FLA_ENABLE_WINDOWS_BUILD
   char
      output_file_m[100];
   
   FILE
      *fpp;
#endif

   fprintf( stdout, "%c Enter number of repeats: ", '%' );
   scanf( "%d", &n_repeats );
   fprintf( stdout, "%c %d\n", '%', n_repeats );

   fprintf( stdout, "%c Enter blocksize: ", '%' );
   scanf( "%u", &nb_alg );
   fprintf( stdout, "%c %u\n", '%', nb_alg );

   fprintf( stdout, "%c Enter problem size parameters: first, inc, num: ", '%' );
   scanf( "%d%d%d", &begin, &increment, &n_trials );
   fprintf( stdout, "%c %d %d %d\n", '%', begin, increment, n_trials );

   fprintf( stdout, "%c Enter number of threads: ", '%' );
   scanf( "%d", &n_threads );
   fprintf( stdout, "%c %d\n", '%', n_threads );

   fprintf( stdout, "%c Enter SuperMatrix parameters: sorting, caching, work stealing, data affinity: ", '%' );
   scanf( "%d%d%d%d", &sorting, &caching, &work_stealing, &data_affinity );
   fprintf( stdout, "%c %s %s %s %s\n\n", '%', ( sorting ? "TRUE" : "FALSE" ), ( caching ? "TRUE" : "FALSE" ), ( work_stealing ? "TRUE" : "FALSE" ), ( data_affinity ? ( data_affinity == 1 ? "FLASH_QUEUE_AFFINITY_2D_BLOCK_CYCLIC" : "FLASH_QUEUE_AFFINITY_OTHER" ) : "FLASH_QUEUE_AFFINITY_NONE" ) );

#ifdef FLA_ENABLE_WINDOWS_BUILD
   fprintf( stdout, "%s_%u = [\n", OUTPUT_FILE, nb_alg );
#else
   sprintf( output_file_m, "%s/%s_output.m", OUTPUT_PATH, OUTPUT_FILE );
   fpp = fopen( output_file_m, "a" );

   fprintf( fpp, "%%\n" );
   fprintf( fpp, "%% | Matrix Size |    FLASH    |\n" );
   fprintf( fpp, "%% |    n x n    |    GFlops   |\n" );
   fprintf( fpp, "%% -----------------------------\n" );
   fprintf( fpp, "%s_%u = [\n", OUTPUT_FILE, nb_alg );
#endif

   FLA_Init();

   dtimes = ( double * ) FLA_malloc( n_repeats * sizeof( double ) );
   flops  = ( double * ) FLA_malloc( n_trials  * sizeof( double ) );
   
   FLASH_Queue_set_num_threads( n_threads );
   FLASH_Queue_set_sorting( sorting );
   FLASH_Queue_set_caching( caching );
   FLASH_Queue_set_work_stealing( work_stealing );
   FLASH_Queue_set_data_affinity( data_affinity );

   for ( i = 0; i < n_trials; i++ )
   {
      size = begin + i * increment;
      
      FLA_Obj_create( datatype, size, size, 0, 0, &A ); 
      FLA_Obj_create( datatype, size, size, 0, 0, &B ); 
      FLA_Obj_create( datatype, size, 1,    0, 0, &x ); 
      FLA_Obj_create( datatype, size, 1,    0, 0, &b ); 
      FLA_Obj_create( datatype, 1,    1,    0, 0, &b_norm ); 
      
      for ( j = 0; j < n_repeats; j++ )
      {
         FLA_Random_matrix( A );
         FLA_Random_matrix( B );
         FLA_Random_matrix( x );
         FLA_Random_matrix( b );

         FLA_Symmetrize( uplo, A );
         FLA_Symmetrize( uplo, B );

         length = ( double ) FLA_Obj_length( B );
         FLA_Add_to_diag( &length, B );
         FLA_Symv_external( uplo, FLA_ONE, B, x, FLA_ZERO, b );

         FLASH_Obj_create_hier_copy_of_flat( A, 1, &nb_alg, &AH );  
         FLASH_Obj_create_hier_copy_of_flat( B, 1, &nb_alg, &BH );  

         FLASH_Chol( uplo, BH );
         
         dtime = FLA_Clock();
         
         FLASH_Eig_gest( inv, uplo, AH, BH );
         
         dtime = FLA_Clock() - dtime;
         dtimes[j] = dtime;
         
         FLASH_Obj_free( &AH );
         FLASH_Obj_free( &BH );
      }
      
      dtime = dtimes[0];
      for ( j = 1; j < n_repeats; j++ )
         dtime = fla_min( dtime, dtimes[j] );
      flops[i] = 1.0 * size * size * size / dtime / 1e9;

#ifdef FLA_ENABLE_WINDOWS_BUILD      
      fprintf( stdout, "   %d   %6.3f   %le\n", size, flops[i], b_norm_value );
#else
      fprintf( fpp, "   %d   %6.3f\n", size, flops[i] );
      
      fprintf( stdout, "Time: %e  |  GFlops: %6.3f\n", dtime, flops[i] );
      fprintf( stdout, "Matrix size: %u x %u  |  nb_alg: %u\n", 
               size, size, nb_alg ); 
      fprintf( stdout, "Norm of difference: %le\n\n", b_norm_value ); 
#endif
 
      FLA_Obj_free( &A ); 
      FLA_Obj_free( &B ); 
      FLA_Obj_free( &x ); 
      FLA_Obj_free( &b ); 
      FLA_Obj_free( &b_norm ); 
   }

#ifdef FLA_ENABLE_WINDOWS_BUILD
   fprintf( stdout, "];\n\n" );
#else
   fprintf( fpp, "];\n" );
   
   fflush( fpp );
   fclose( fpp );
#endif

   FLA_free( dtimes );
   FLA_free( flops );

   FLA_Finalize(); 
   
   return 0; 
}
