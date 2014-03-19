/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#if FLA_MULTITHREADING_MODEL == FLA_OPENMP

#include <omp.h>

extern omp_lock_t* fla_omp_lock;

FLA_Error FLA_omp_init( void );
void      FLA_omp_finalize( void );
int       FLA_omp_get_num_threads( void );
void      FLA_omp_set_num_threads( int n_threads );
int       FLA_omp_get_num_stages( void );
void      FLA_omp_set_num_stages( int n_stages );
int       FLA_omp_compute_stage_width( FLA_Obj A );

#endif
