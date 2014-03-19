/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

void FLA_Task_partitioning_init( void );
void FLA_Task_partitioning_set( int n_threads, int tag0_val, int tag1_val, int tag2_val, int tag3_val, int tag4_val );

int  FLA_Task_compute_blocksize( int tag, FLA_Obj A, FLA_Obj A_proc, FLA_Quadrant from );

int  FLA_task_get_num_partitions( int n_threads, int tag );

int  FLA_task_determine_matrix_size( FLA_Obj A, FLA_Quadrant from );
int  FLA_task_determine_relative_blocksize( int A_size, int A_proc_size, int n_part );
int  FLA_task_determine_absolute_blocksize( int A_size, int A_proc_size, int nb_alg );

