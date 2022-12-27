/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


#if FLA_MULTITHREADING_MODEL == FLA_OPENMP

#define N_OMP_LOCKS 64

omp_lock_t* fla_omp_lock;
int         n_omp_stages;

FLA_Error FLA_omp_init()
{
	int i;
	
	// Allocate storage for the locks.
	fla_omp_lock = (omp_lock_t*)malloc( N_OMP_LOCKS * sizeof(omp_lock_t) );

	// Initialize the locks.
	for( i = 0; i < N_OMP_LOCKS; ++i )
		omp_init_lock( &fla_omp_lock[i] );
	
	// Set the default number of lockable "stages" to one (ie: the entire
	// matrix is locked).
	FLA_omp_set_num_stages( 1 );

	return FLA_SUCCESS;
}

void FLA_omp_finalize()
{
	// Free the locks.
	free( fla_omp_lock );
}

int FLA_omp_get_num_threads()
{
	return omp_get_num_threads();
}

void FLA_omp_set_num_threads( int n_threads )
{
	omp_set_num_threads( n_threads );
}

int FLA_omp_get_num_stages()
{
	return n_omp_stages;
}

void FLA_omp_set_num_stages( int n_stages )
{
	n_omp_stages = n_stages;
}

int FLA_omp_compute_stage_width( FLA_Obj A )
{
	int n_stages = FLA_omp_get_num_stages();
	int stage_width;
	
	// Compute the desired width of one lockable panel partition.
	if( n_stages == 1 )
		stage_width = FLA_Obj_width(A);
	else if( FLA_Obj_width(A) % n_stages == 0 )
		stage_width = fla_max( FLA_Obj_width(A)/n_stages, 1 );
	else
		stage_width = fla_max( FLA_Obj_width(A)/n_stages + 1, 1 );
	
	return stage_width;
}

#endif

