
#include "FLAME.h"
#include "FLA_queue_omp.h"
#include "FLA_task_partitioning.h"

#if FLA_MULTITHREADING_MODEL == FLA_OPENMP

#define NUM_THREADS_MAX 128
#define NUM_TAGS_MAX      5

int factors[NUM_THREADS_MAX+1][NUM_TAGS_MAX];

void FLA_Task_partitioning_init()
{
	int i;
	
	for( i = 0; i < NUM_THREADS_MAX; ++i )
		FLA_Task_partitioning_set( i, -1, -1, -1, -1, -1 );
}

void FLA_Task_partitioning_set( int n_threads, int tag0_val, int tag1_val, int tag2_val, int tag3_val, int tag4_val )
{
	factors[n_threads][0] = tag0_val;
	factors[n_threads][1] = tag1_val;
	factors[n_threads][2] = tag2_val;
	factors[n_threads][3] = tag3_val;
	factors[n_threads][4] = tag4_val;
}

int FLA_task_get_num_partitions( int n_threads, int tag )
{
	return factors[n_threads][tag];
}

int FLA_Task_compute_blocksize( int tag, FLA_Obj A, FLA_Obj A_proc, FLA_Quadrant from )
{
	int n_threads = FLA_Queue_get_num_threads();
	int A_size, A_proc_size;
	int n_part;
	int b;
	
	// Determine the sizes of the matrix partitions.
	A_size      = FLA_task_determine_matrix_size( A, from );
	A_proc_size = FLA_task_determine_matrix_size( A_proc, from );
	
	// Determine the raw blocksize value.
	n_part      = FLA_task_get_num_partitions( n_threads, tag );
	
	// Determine the blocksize based on the sign of the value from
	// _get_num_partitions().
	if( n_part > 0 )
	{
		b = FLA_task_determine_absolute_blocksize( A_size,
	                                               A_proc_size,
	                                               n_part );
	}
	else if( n_part < 0 )
	{
	    b = FLA_task_determine_relative_blocksize( A_size,
	                                               A_proc_size,
	                                               abs(n_part) );
	}
	else
	{
		FLA_Print_message( "Detected blocksize of 0!", __FILE__, __LINE__ );
        FLA_Abort();
	}

	return b;
}

int FLA_task_determine_matrix_size( FLA_Obj A, FLA_Quadrant from )
{
	int r_val = 0;
	
	// Determine the size of the matrix dimension along which we are moving.
	switch( from )
	{
		case FLA_TOP:
		case FLA_BOTTOM:
		{
			r_val = FLA_Obj_length( A );
			break;
		}
		case FLA_LEFT:
		case FLA_RIGHT:
		{
			r_val = FLA_Obj_width( A );
			break;
		}
		case FLA_TL:
		case FLA_TR:
		case FLA_BL:
		case FLA_BR:
		{
			// If A happens to be the full object, we need to use min_dim() here
			// because the matrix might be rectangular. If A is the processed
			// partition, it is very probably square, and min_dim() doesn't hurt.
			r_val = FLA_Obj_min_dim( A );
			break;
		}
		default:
			FLA_Print_message( "Unexpected default in switch statement!", __FILE__, __LINE__ );
			FLA_Abort();
	}

	return r_val;
}


int FLA_task_determine_relative_blocksize( int A_size, int A_proc_size, int n_part )
{
	int b, i, z;

    // Return early if the size is zero (scenario (A)).
    if( A_size == 0 ) return 0;

    // Compute the base blocksize (according to (1) above).
    b = A_size / n_part;

    // If the base blocksize is zero (because A_size < n_part), then override
    // n_part and proceed as if only one partition was requested (scenario (B)).
	if( b == 0 )
    {
        n_part = 1;
        b      = A_size;
    }

    // Compute partition index i.
    i = A_proc_size / b;

    // Compute the index z below which we will use the base blocksize.
    z = n_part - (A_size % n_part);

    // If the current partition index i is at least z, then increment the
    // base blocksize (according to (2) above).
    if( z <= i ) b++;

    // Return the blocksize.
	return b;
}

int FLA_task_determine_absolute_blocksize( int A_size, int A_proc_size, int nb_alg )
{
    int A_unproc_size = A_size - A_proc_size;
    int b;

    b = min( A_unproc_size, nb_alg );

    return b;
}

#endif
