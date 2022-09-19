/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_HIP

#include "hip/hip_runtime_api.h"
#include "rocblas/rocblas.h"
#include "rocsolver/rocsolver.h"

FLA_Error FLA_Bidiag_form_U_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip )
{
  FLA_Datatype datatype;
  int          m_A, n_A, k_A;
  int          cs_A;
  rocblas_storev blas_vect = rocblas_column_wise;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Bidiag_form_U_check( A, t );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  if ( blas_vect == rocblas_column_wise ) k_A = FLA_Obj_vector_dim( t );
  else                    k_A = FLA_Obj_vector_dim( t ) + 1;

  void* A_vecs = NULL;
  void* t_scals = NULL;
  if ( FLASH_Queue_get_malloc_managed_enabled_hip() )
  {
    A_vecs = FLA_Obj_buffer_at_view( A );
    t_scals = FLA_Obj_buffer_at_view( t );
  }
  else
  {
    A_vecs = A_hip;
    t_scals = t_hip;
  }


  switch( datatype ){

    case FLA_FLOAT:
    {
      float*    buff_A    = ( float    * ) A_vecs;
      float*    buff_t    = ( float    * ) t_scals;

      rocsolver_sorgbr( handle,
                        blas_vect,
                        m_A,
                        n_A,
                        k_A,
                        buff_A, cs_A,
                        buff_t );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A    = ( double   * ) A_vecs;
      double*   buff_t    = ( double   * ) t_scals;

      rocsolver_dorgbr( handle,
                        blas_vect,
                        m_A,
                        n_A,
                        k_A,
                        buff_A, cs_A,
                        buff_t );

      break;
    }

    case FLA_COMPLEX:
    {
      rocblas_float_complex* buff_A    = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( A );
      rocblas_float_complex* buff_t    = ( rocblas_float_complex * ) FLA_COMPLEX_PTR( t );

      rocsolver_cungbr( handle,
                        blas_vect,
                        m_A,
                        n_A,
                        k_A,
                        buff_A, cs_A,
                        buff_t );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      rocblas_double_complex *buff_A    = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      rocblas_double_complex *buff_t    = ( rocblas_double_complex * ) FLA_DOUBLE_COMPLEX_PTR( t );

      rocsolver_zungbr( handle,
                        blas_vect,
                        m_A,
                        n_A,
                        k_A,
                        buff_A, cs_A,
                        buff_t );

      break;
    }

  }


  return FLA_SUCCESS;
}

#endif
