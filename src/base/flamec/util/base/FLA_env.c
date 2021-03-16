/*
    Copyright (c) 2021 Advanced Micro Devices, Inc.  All rights reserved.
    Mar 16, 2021
*/

#include "FLAME.h"

int FLA_env_get_var( const char* env, int fallback )
{
	int r_val;
	char*  str;

	// Query the environment variable and store the result in str.
	str = getenv( env );

	// Set the return value based on the string obtained from getenv().
	if ( str != NULL )
	{
		// If there was no error, convert the string to an integer and
		// prepare to return that integer.
		r_val = ( int )strtol( str, NULL, 10 );
	}
	else
	{
		// If there was an error, use the "fallback" as the return value.
		r_val = fallback;
	}

	return r_val;
}

int FLASH_get_num_threads( int fallback )
{
    int omp_threads, fla_threads;

    omp_threads = FLA_env_get_var( "OMP_NUM_THREADS", -1 );
    fla_threads = FLA_env_get_var( "FLA_NUM_THREADS", -1 );

    if( ( omp_threads == -1 ) && ( fla_threads == -1 ) )
    {
        fla_threads = fallback;
    }
    else if( fla_threads == -1 )
    {
        fla_threads = omp_threads;
    }

    if( fla_threads <= 0 )
    {
        fla_threads = fallback;
    }

    return fla_threads;
}

