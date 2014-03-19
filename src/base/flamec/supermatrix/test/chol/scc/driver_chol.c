/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h" 
//extern "C" {
#include "RCCE.h"
//}

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

//Change the user name in the path below to yours
#define OUTPUT_PATH "/shared/bamarker/results"
#define OUTPUT_FILE "chol"

void Synch_all() {
  RCCE_barrier( &RCCE_COMM_WORLD );
}


//extern "C" {
int RCCE_APP( int argc, char* argv[] )
{ 
   int
      i, j,
      rank,
      n_threads,
      n_repeats,
      n_trials,
      increment,
      begin,
      sorting = 0,
      caching = 0,
      work_stealing = 0,
      data_affinity = 0;

   dim_t
      size,
      nb_alg;

   FLA_Datatype
      datatype = FLA_DOUBLE;

   FLA_Obj 
      A, x, b, b_norm,
      AH;
   
   double 
      length,
      b_norm_value,
      dtime, 
      *dtimes,
      flops;

   char
      output_file_m[100];

   FILE
      *fpp;

   RCCE_init( &argc, &argv );
   RCCE_comm_rank( RCCE_COMM_WORLD, &rank );
   RCCE_comm_size( RCCE_COMM_WORLD, &n_threads );

   if ( argc < 6 ) 
   {
      if ( rank == 0 ) 
         fprintf( stdout, "args: %d  |  Usage: <number of repeats> <block size> <first problem size> <problem size increment> <number of trials>\n", argc );
      
      return 0;
   }
   
   n_repeats = atoi( argv[1] );
   nb_alg    = atoi( argv[2] );
   begin     = atoi( argv[3] );
   increment = atoi( argv[4] );
   n_trials  = atoi( argv[5] );

   if ( rank == 0 )
   {
      sprintf( output_file_m, "%s/%s_output.m", OUTPUT_PATH, OUTPUT_FILE );
      fpp = fopen( output_file_m, "a" );

      fprintf( fpp, "%%\n" );
      fprintf( fpp, "%% | Matrix Size |    FLASH    |\n" );
      fprintf( fpp, "%% |    n x n    |    GFLOPS   |\n" );
      fprintf( fpp, "%% -----------------------------\n" );
      fprintf( fpp, "%s_%u = [\n", OUTPUT_FILE, nb_alg );
   }

#ifdef FLA_ENABLE_SCC   
   FLA_Init();
   
   FLASH_Queue_set_num_threads( n_threads );
   FLASH_Queue_set_sorting( sorting );
   FLASH_Queue_set_caching( caching );
   FLASH_Queue_set_work_stealing( work_stealing );
   FLASH_Queue_set_data_affinity( data_affinity );
   
   dtimes = ( double * ) FLA_malloc( n_repeats * sizeof( double ) );
   
   for ( i = 0; i < n_trials; i++ )
   {
      size = begin + i * increment;
      
      FLA_Obj_create( datatype, size, size, 0, 0, &A ); 
      FLA_Obj_create( datatype, size, 1,    0, 0, &x ); 
      FLA_Obj_create( datatype, size, 1,    0, 0, &b ); 
      FLA_Obj_create( datatype, 1,    1,    0, 0, &b_norm ); 
      
      for ( j = 0; j < n_repeats; j++ )
      {
         if ( rank == 0 )
         {
            FLA_Random_matrix( A );
            FLA_Random_matrix( x );
            FLA_Random_matrix( b );
            
            length = ( double ) FLA_Obj_length( A );
            FLA_Add_to_diag( &length, A );
            FLA_Symv_external( FLA_LOWER_TRIANGULAR, FLA_ONE, A, x, FLA_ZERO, b );
         }
         
         RCCE_barrier( &RCCE_COMM_WORLD );
         
         FLASH_Obj_create_hier_copy_of_flat( A, 1, &nb_alg, &AH );  

         RCCE_barrier( &RCCE_COMM_WORLD );
         
         dtime = RCCE_wtime();
         
         FLASH_Chol( FLA_LOWER_TRIANGULAR, AH );
         
	 RCCE_barrier( &RCCE_COMM_WORLD );
         	 
         dtime = RCCE_wtime() - dtime;
         dtimes[j] = dtime;

         FLASH_Obj_flatten( AH, A );         

         RCCE_barrier( &RCCE_COMM_WORLD );

         FLASH_Obj_free( &AH );
      }
      
      dtime = dtimes[0];
      for ( j = 1; j < n_repeats; j++ )
         dtime = min( dtime, dtimes[j] );
      flops = 1.0 / 3.0 * size * size * size / dtime / 1e9;
      
      if ( rank == 0 ) 
      {
         FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, 
                            FLA_NONUNIT_DIAG, A, b );
         FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, 
                            FLA_NONUNIT_DIAG, A, b ); 
         FLA_Axpy_external( FLA_MINUS_ONE, x, b );
         FLA_Nrm2_external( b, b_norm );
         FLA_Obj_extract_real_scalar( b_norm, &b_norm_value );

         fprintf( fpp, "   %d   %6.3f\n", size, flops );
         
         fprintf( stdout, "\nTime: %e  |  GFLOPS: %6.3f\n", dtime, flops );
         fprintf( stdout, "Matrix size: %u x %u  |  Block size: %u\n", 
                  size, size, nb_alg ); 
         fprintf( stdout, "Norm of difference: %le\n", b_norm_value ); 
      }
      
      FLA_Obj_free( &A ); 
      FLA_Obj_free( &x ); 
      FLA_Obj_free( &b ); 
      FLA_Obj_free( &b_norm ); 
   }

   FLA_free( dtimes );
   
   FLA_Finalize(); 
#else
   fprintf( stdout, "libflame is not configured for SCC.\n" );
#endif

   if ( rank == 0 )
   {
      fprintf( fpp, "];\n" );

      fflush( fpp );
      fclose( fpp );
   }
   
   RCCE_finalize();
   
   return 0; 
}
//}
