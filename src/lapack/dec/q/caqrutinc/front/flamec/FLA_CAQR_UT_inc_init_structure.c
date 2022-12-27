/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

void FLA_CAQR_UT_inc_init_structure( dim_t p, dim_t nb_part, FLA_Obj R )
{
	dim_t    m, n;
	dim_t    rs, cs;
	dim_t    i, j, ip;
	FLA_Obj* buff_R;

	m      = FLA_Obj_length( R );
	n      = FLA_Obj_width( R );
	rs     = FLA_Obj_row_stride( R );
	cs     = FLA_Obj_col_stride( R );
	buff_R = FLA_Obj_buffer_at_view( R );

	// Fill in R by row panels.
	for ( ip = 0; ip < p; ++ip )
	{
		FLA_Obj* buff_R1 = buff_R + (ip*nb_part)*rs;

		integer  m_behind   = ip*nb_part;
		integer  m_ahead    = m - m_behind;

		integer  m_cur      = fla_min( nb_part, m_ahead );
		integer  n_cur      = n;

		// Iterate across columns for the current panel.
		for ( j = 0; j < n_cur; ++j )
		{
			FLA_Obj* rho = buff_R1 + j*cs;

			// Mark the above-diagonal blocks as full.
			for ( i = 0; i < j; ++i )
			{
				rho->base->uplo = FLA_FULL_MATRIX;
				rho += rs;
			}

			// Mark the diagonal block as triangular.
			rho->base->uplo = FLA_UPPER_TRIANGULAR;
			rho += rs;
			
			// Mark the below-diagonal blocks as zero.
			for ( i = j + 1; i < m_cur; ++i )
			{
				rho->base->uplo = FLA_ZERO_MATRIX;
				rho += rs;
			}
		}
	}
}

