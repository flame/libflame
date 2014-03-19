/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_shift_U( FLA_Uplo uplo, FLA_Obj A )
{
	FLA_Datatype datatype;
	int          m_A;
	int          rs_A, cs_A;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Tridiag_UT_shift_U_check( uplo, A );

	datatype = FLA_Obj_datatype( A );

        // Play with swapping of cs rs; we do not need "u" version.
        if ( uplo == FLA_LOWER_TRIANGULAR ) 
          {
            m_A      = FLA_Obj_length( A );
            rs_A     = FLA_Obj_row_stride( A );
            cs_A     = FLA_Obj_col_stride( A );
          }
        else
          {
            m_A      = FLA_Obj_width( A );
            cs_A     = FLA_Obj_row_stride( A );
            rs_A     = FLA_Obj_col_stride( A );
          }

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_A = ( float* ) FLA_FLOAT_PTR( A );
                        FLA_Tridiag_UT_shift_U_l_ops( m_A,
                                                      buff_A, rs_A, cs_A );
			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_A = ( double* ) FLA_DOUBLE_PTR( A );
                        FLA_Tridiag_UT_shift_U_l_opd( m_A,
                                                      buff_A, rs_A, cs_A );
			break;
		}

		case FLA_COMPLEX:
		{
			scomplex* buff_A = ( scomplex* ) FLA_COMPLEX_PTR( A );
                        FLA_Tridiag_UT_shift_U_l_opc( m_A,
                                                      buff_A, rs_A, cs_A );
			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			dcomplex* buff_A = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );
                        FLA_Tridiag_UT_shift_U_l_opz( m_A,
                                                      buff_A, rs_A, cs_A );
			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Tridiag_UT_shift_U_l_ops( int       m_A,
                                        float*    buff_A, int rs_A, int cs_A )
{
	float*  a00  = buff_A;
	float*  a10  = buff_A + rs_A;
	float   zero = bl1_s0();
	float   one  = bl1_s1();
	int     j;

	for ( j = m_A - 1; j > 0; --j )
	{
		float*    alpha01  = buff_A + (j  )*cs_A + (0  )*rs_A;
		float*    alpha11  = buff_A + (j  )*cs_A + (j  )*rs_A;
		float*    a20      = buff_A + (j-1)*cs_A + (j+1)*rs_A;
		float*    a21      = buff_A + (j  )*cs_A + (j+1)*rs_A;

		int       m_ahead  = m_A - j - 1;

		*alpha01 = zero;
		*alpha11 = one;
		bl1_scopyv( BLIS1_NO_CONJUGATE,
		            m_ahead,
		            a20, rs_A,
		            a21, rs_A );
	}

	*a00 = one;
	bl1_ssetv( m_A - 1,
	           &zero,
	           a10, rs_A );

	return FLA_SUCCESS;
}

FLA_Error FLA_Tridiag_UT_shift_U_l_opd( int       m_A,
                                        double*   buff_A, int rs_A, int cs_A )
{
	double* a00  = buff_A;
	double* a10  = buff_A + rs_A;
	double  zero = bl1_d0();
	double  one  = bl1_d1();
	int     j;

	for ( j = m_A - 1; j > 0; --j )
	{
		double*   alpha01  = buff_A + (j  )*cs_A + (0  )*rs_A;
		double*   alpha11  = buff_A + (j  )*cs_A + (j  )*rs_A;
		double*   a20      = buff_A + (j-1)*cs_A + (j+1)*rs_A;
		double*   a21      = buff_A + (j  )*cs_A + (j+1)*rs_A;

		int       m_ahead  = m_A - j - 1;

		*alpha01 = zero;
		*alpha11 = one;
		bl1_dcopyv( BLIS1_NO_CONJUGATE,
		            m_ahead,
		            a20, rs_A,
		            a21, rs_A );
	}

	*a00 = one;
	bl1_dsetv( m_A - 1,
	           &zero,
	           a10, rs_A );

	return FLA_SUCCESS;
}

FLA_Error FLA_Tridiag_UT_shift_U_l_opc( int       m_A,
                                        scomplex* buff_A, int rs_A, int cs_A )
{
	scomplex* a00  = buff_A;
	scomplex* a10  = buff_A + rs_A;
	scomplex  zero = bl1_c0();
	scomplex  one  = bl1_c1();
	int       j;

	for ( j = m_A - 1; j > 0; --j )
	{
		scomplex* alpha01  = buff_A + (j  )*cs_A + (0  )*rs_A;
		scomplex* alpha11  = buff_A + (j  )*cs_A + (j  )*rs_A;
		scomplex* a20      = buff_A + (j-1)*cs_A + (j+1)*rs_A;
		scomplex* a21      = buff_A + (j  )*cs_A + (j+1)*rs_A;

		int       m_ahead  = m_A - j - 1;

		*alpha01 = zero;
		*alpha11 = one;
		bl1_ccopyv( BLIS1_NO_CONJUGATE,
		            m_ahead,
		            a20, rs_A,
		            a21, rs_A );
	}

	*a00 = one;
	bl1_csetv( m_A - 1,
	           &zero,
	           a10, rs_A );

	return FLA_SUCCESS;
}

FLA_Error FLA_Tridiag_UT_shift_U_l_opz( int       m_A,
                                        dcomplex* buff_A, int rs_A, int cs_A )
{
	dcomplex* a00  = buff_A;
	dcomplex* a10  = buff_A + rs_A;
	dcomplex  zero = bl1_z0();
	dcomplex  one  = bl1_z1();
	int       j;

	for ( j = m_A - 1; j > 0; --j )
	{
		dcomplex* alpha01  = buff_A + (j  )*cs_A + (0  )*rs_A;
		dcomplex* alpha11  = buff_A + (j  )*cs_A + (j  )*rs_A;
		dcomplex* a20      = buff_A + (j-1)*cs_A + (j+1)*rs_A;
		dcomplex* a21      = buff_A + (j  )*cs_A + (j+1)*rs_A;

		int       m_ahead  = m_A - j - 1;

		*alpha01 = zero;
		*alpha11 = one;
		bl1_zcopyv( BLIS1_NO_CONJUGATE,
		            m_ahead,
		            a20, rs_A,
		            a21, rs_A );
	}

	*a00 = one;
	bl1_zsetv( m_A - 1,
	           &zero,
	           a10, rs_A );

	return FLA_SUCCESS;
}

