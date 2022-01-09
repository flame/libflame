/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

integer fla_scomp_f( const void* a, const void* b );
integer fla_scomp_b( const void* a, const void* b );
integer fla_dcomp_f( const void* a, const void* b );
integer fla_dcomp_b( const void* a, const void* b );

FLA_Error FLA_Sort( FLA_Direct direct, FLA_Obj x )
{
	FLA_Datatype datatype;
	FLA_Obj      x_use;
	dim_t        m_x;
	dim_t        inc_x;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Sort_check( direct, x );

	datatype = FLA_Obj_datatype( x );

	m_x      = FLA_Obj_vector_dim( x );
	inc_x    = FLA_Obj_vector_inc( x );

	// If the vector does not have unit stride, copy it to a temporary vector
	// that does have unit stride.
	if ( inc_x != 1 )
	{
		FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, x, &x_use );
		inc_x = FLA_Obj_vector_inc( x_use );
	}
	else
	{
		x_use = x;
	}

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float* x_p = ( float* ) FLA_FLOAT_PTR( x_use );

			if ( direct == FLA_FORWARD )
				FLA_Sort_f_ops( m_x,
				                x_p, inc_x );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_b_ops( m_x,
				                x_p, inc_x );

			break;
		}

		case FLA_DOUBLE:
		{
			double* x_p = ( double* ) FLA_DOUBLE_PTR( x_use );

			if ( direct == FLA_FORWARD )
				FLA_Sort_f_opd( m_x,
				                x_p, inc_x );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_b_opd( m_x,
				                x_p, inc_x );

			break;
		}

	}

	if ( inc_x != 1 )
	{
		FLA_Copy( x_use, x );
		FLA_Obj_free( &x_use );
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Sort_f_ops( integer     m_x,
                          float*  x, integer inc_x )
{
	qsort( x,
	       m_x,
	       sizeof( float ),
	       fla_scomp_f );

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_b_ops( integer     m_x,
                          float*  x, integer inc_x )
{
	qsort( x,
	       m_x,
	       sizeof( float ),
	       fla_scomp_b );

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_f_opd( integer     m_x,
                          double* x, integer inc_x )
{
	qsort( x,
	       m_x,
	       sizeof( double ),
	       fla_dcomp_f );

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_b_opd( integer     m_x,
                          double* x, integer inc_x )
{
	qsort( x,
	       m_x,
	       sizeof( double ),
	       fla_dcomp_b );

	return FLA_SUCCESS;
}





integer fla_scomp_f( const void* a, const void* b )
{
	float*  da = ( float*  ) a;
	float*  db = ( float*  ) b;
	integer     r_val;

	if      ( *da < *db ) r_val = -1;
	else if ( *da > *db ) r_val =  1;
	else                  r_val =  0;

	return r_val;
}

integer fla_scomp_b( const void* a, const void* b )
{
	float*  da = ( float*  ) a;
	float*  db = ( float*  ) b;
	integer     r_val;

	if      ( *da < *db ) r_val =  1;
	else if ( *da > *db ) r_val = -1;
	else                  r_val =  0;

	return r_val;
}

integer fla_dcomp_f( const void* a, const void* b )
{
	double* da = ( double* ) a;
	double* db = ( double* ) b;
	integer     r_val;

	if      ( *da < *db ) r_val = -1;
	else if ( *da > *db ) r_val =  1;
	else                  r_val =  0;

	return r_val;
}

integer fla_dcomp_b( const void* a, const void* b )
{
	double* da = ( double* ) a;
	double* db = ( double* ) b;
	integer     r_val;

	if      ( *da < *db ) r_val =  1;
	else if ( *da > *db ) r_val = -1;
	else                  r_val =  0;

	return r_val;
}

