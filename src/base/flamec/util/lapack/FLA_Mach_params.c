/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Mach_params( FLA_Machval machval, FLA_Obj val )
{
	FLA_Datatype datatype;

	datatype = FLA_Obj_datatype( val );

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Mach_params_check( machval, val );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*  val_p = ( float* ) FLA_FLOAT_PTR( val );

			*val_p = FLA_Mach_params_ops( machval );

			break;
		}

		case FLA_DOUBLE:
		{
			double* val_p = ( double* ) FLA_DOUBLE_PTR( val );

			*val_p = FLA_Mach_params_opd( machval );

			break;
		}
	}

	return FLA_SUCCESS;
}


float FLA_Mach_params_ops( FLA_Machval machval )
{
	// No cache: the previous version filled a static table on first use,
	// guarded by an unsynchronized first_time flag, so concurrent first
	// callers could read a partially filled table. fla_slamch() is now a
	// stateless switch over <float.h> constants, so computing on demand is
	// both thread-safe and cheap.
	char lapack_machval;
	float  eps;

	if ( machval == FLA_MACH_EPS2 )
	{
		// epsilon^2 has no LAPACK character code.
		FLA_Param_map_flame_to_netlib_machval( FLA_MACH_EPS, &lapack_machval );
		eps = fla_slamch( &lapack_machval, 1 );
		return eps * eps;
	}

	FLA_Param_map_flame_to_netlib_machval( machval, &lapack_machval );
	return fla_slamch( &lapack_machval, 1 );
}

double FLA_Mach_params_opd( FLA_Machval machval )
{
	// No cache: the previous version filled a static table on first use,
	// guarded by an unsynchronized first_time flag, so concurrent first
	// callers could read a partially filled table. fla_dlamch() is now a
	// stateless switch over <float.h> constants, so computing on demand is
	// both thread-safe and cheap.
	char lapack_machval;
	double  eps;

	if ( machval == FLA_MACH_EPS2 )
	{
		// epsilon^2 has no LAPACK character code.
		FLA_Param_map_flame_to_netlib_machval( FLA_MACH_EPS, &lapack_machval );
		eps = fla_dlamch( &lapack_machval, 1 );
		return eps * eps;
	}

	FLA_Param_map_flame_to_netlib_machval( machval, &lapack_machval );
	return fla_dlamch( &lapack_machval, 1 );
}

