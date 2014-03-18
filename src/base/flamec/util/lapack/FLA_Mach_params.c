
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
	static int    first_time = TRUE;
	static float  vals[FLA_MACH_N_VALS];

	if ( first_time )
	{
		char lapack_machval;
		int  i;

		for( i = 0; i < FLA_MACH_N_VALS - 1; ++i )
		{
			FLA_Param_map_flame_to_netlib_machval( FLA_MACH_START + i, &lapack_machval );
//printf( "querying %d %c\n", FLA_MACH_START + i, lapack_machval );
			vals[i] = fla_slamch( &lapack_machval, 1 );
//printf( "got back  %34.29e\n", vals[i] );
		}

		// Store epsilon^2 in the last element.
		vals[i] = vals[0] * vals[0];

		first_time = FALSE;
	}

	return vals[ machval - FLA_MACH_START ];
}

double FLA_Mach_params_opd( FLA_Machval machval )
{
	static int    first_time = TRUE;
	static double vals[FLA_MACH_N_VALS];

	if ( first_time )
	{
		char lapack_machval;
		int  i;

		for( i = 0; i < FLA_MACH_N_VALS - 1; ++i )
		{
			FLA_Param_map_flame_to_netlib_machval( FLA_MACH_START + i, &lapack_machval );
//printf( "querying %d %c\n", FLA_MACH_START + i, lapack_machval );
			vals[i] = fla_dlamch( &lapack_machval, 1 );
//printf( "got back  %34.29e\n", vals[i] );
		}

		// Store epsilon^2 in the last element.
		vals[i] = vals[0] * vals[0];

		first_time = FALSE;
	}

	return vals[ machval - FLA_MACH_START ];
}

