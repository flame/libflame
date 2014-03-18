
#include "FLAME.h"

#define FLA_SORT_BSVD_EXT_DEFINE_OBJ_VARIABLES( A, apply, dtype, m, measure, rs, cs ) \
  {                                                                     \
    if ( apply == TRUE ) {                                              \
      dtype = FLA_Obj_datatype( A );                                    \
      m     = measure ( A );                                            \
      rs    = FLA_Obj_row_stride( A );                                  \
      cs    = FLA_Obj_col_stride( A );                                  \
    } else {                                                            \
      m     = 0;                                                        \
      cs    = 0;                                                        \
      rs    = 0;                                                        \
    }                                                                   \
  }

#define FLA_SORT_BSVD_EXT_COMPARE_FORWARD( a, b ) ( a < b )
#define FLA_SORT_BSVD_EXT_COMPARE_BACKWARD( a, b ) ( a > b )
#define FLA_SORT_BSVD_EXT_BODY( direct, vector_swap )                   \
  {                                                                     \
    for ( ii = 1; ii < m_s; ++ii ) {                                    \
      i = ii - 1;                                                       \
      k = i;                                                            \
      p = s[ i*inc_s ];                                                 \
                                                                        \
      for ( j = ii; j < m_s; ++j ) {                                    \
        if ( FLA_SORT_BSVD_EXT_COMPARE_ ## direct (s[ j*inc_s ], p) ) { \
          k = j;                                                        \
          p = s[ j*inc_s ];                                             \
        }                                                               \
      }                                                                 \
                                                                        \
      if ( k != i ) {                                                   \
        /* printf("bsvd sort: exchange %d %d\n", i,k);*/                \
        s[ k*inc_s ] = s[ i ];                                          \
        s[ i       ] = p;                                               \
                                                                        \
        if ( U != NULL )                                                \
          vector_swap( m_U,                                             \
                       U + i*cs_U, rs_U,                                \
                       U + k*cs_U, rs_U );                              \
        if ( V != NULL )                                                \
          vector_swap( m_V,                                             \
                       V + i*cs_V, rs_V,                                \
                       V + k*cs_V, rs_V );                              \
        if ( C != NULL )                                                \
          vector_swap( n_C,                                             \
                       C + i*rs_C, cs_C,                                \
                       C + k*rs_C, cs_C );                              \
      }                                                                 \
    }                                                                   \
  }


// According to the sorted order of a given vector s,
// U and V are reordered in columns while C is reordered
// in rows when they need to be applied.
FLA_Error FLA_Sort_bsvd_ext( FLA_Direct direct, FLA_Obj s,
                             FLA_Bool apply_U, FLA_Obj U,
                             FLA_Bool apply_V, FLA_Obj V,
                             FLA_Bool apply_C, FLA_Obj C )
{
    FLA_Datatype datatype;
    dim_t        m_U, rs_U, cs_U;
    dim_t        m_V, rs_V, cs_V;
    dim_t        n_C, rs_C, cs_C;
    dim_t        m_s, inc_s;

    //if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    //	FLA_Sort_bsvd_check( direct, s,
    //                         apply_U, U,
    //                         apply_V, V,
    //                         apply_C, C );

    // Sort singular values only; quick sort
    if ( apply_U == FALSE && apply_V == FALSE  )
        return FLA_Sort( direct, s );

    // s dimensions must be provided.
    m_s      = FLA_Obj_vector_dim( s );
    inc_s    = FLA_Obj_vector_inc( s );

    // Datatype of U, V and C must be consistent and must be defined from one of them.
    FLA_SORT_BSVD_EXT_DEFINE_OBJ_VARIABLES( U, apply_U, datatype, m_U, FLA_Obj_length, rs_U, cs_U );
    FLA_SORT_BSVD_EXT_DEFINE_OBJ_VARIABLES( V, apply_V, datatype, m_V, FLA_Obj_length, rs_V, cs_V );
    FLA_SORT_BSVD_EXT_DEFINE_OBJ_VARIABLES( C, apply_C, datatype, n_C, FLA_Obj_width,  rs_C, cs_C );
    
    switch ( datatype )
    {
    case FLA_FLOAT:
    {
        float* s_p = ( float* ) FLA_FLOAT_PTR( s );
        float* U_p = ( apply_U == TRUE ? ( float* ) FLA_FLOAT_PTR( U ) : NULL );
        float* V_p = ( apply_V == TRUE ? ( float* ) FLA_FLOAT_PTR( V ) : NULL );
        float* C_p = ( apply_C == TRUE ? ( float* ) FLA_FLOAT_PTR( C ) : NULL );

        if ( direct == FLA_FORWARD )
            FLA_Sort_bsvd_ext_f_ops( m_s, s_p, inc_s,
                                     m_U, U_p, rs_U, cs_U,
                                     m_V, V_p, rs_V, cs_V,
                                     n_C, C_p, rs_C, cs_C );
        else // if ( direct == FLA_BACKWARD )
            FLA_Sort_bsvd_ext_b_ops( m_s, s_p, inc_s,
                                     m_U, U_p, rs_U, cs_U,
                                     m_V, V_p, rs_V, cs_V,
                                     n_C, C_p, rs_C, cs_C );
        break;
    }
    case FLA_DOUBLE:
    {
        double* s_p = ( double* ) FLA_DOUBLE_PTR( s );
        double* U_p = ( apply_U == TRUE ? ( double* ) FLA_DOUBLE_PTR( U ) : NULL );
        double* V_p = ( apply_V == TRUE ? ( double* ) FLA_DOUBLE_PTR( V ) : NULL );
        double* C_p = ( apply_C == TRUE ? ( double* ) FLA_DOUBLE_PTR( C ) : NULL );

        if ( direct == FLA_FORWARD )
            FLA_Sort_bsvd_ext_f_opd( m_s, s_p, inc_s,
                                     m_U, U_p, rs_U, cs_U,
                                     m_V, V_p, rs_V, cs_V,
                                     n_C, C_p, rs_C, cs_C );
        else // if ( direct == FLA_BACKWARD )
            FLA_Sort_bsvd_ext_b_opd( m_s, s_p, inc_s,
                                     m_U, U_p, rs_U, cs_U,
                                     m_V, V_p, rs_V, cs_V,
                                     n_C, C_p, rs_C, cs_C );
        break;
    }
    case FLA_COMPLEX:
    {
        float*    s_p = ( float*    ) FLA_FLOAT_PTR( s );
        scomplex* U_p = ( apply_U == TRUE ? ( scomplex* ) FLA_COMPLEX_PTR( U ) : NULL );
        scomplex* V_p = ( apply_V == TRUE ? ( scomplex* ) FLA_COMPLEX_PTR( V ) : NULL );
        scomplex* C_p = ( apply_C == TRUE ? ( scomplex* ) FLA_COMPLEX_PTR( C ) : NULL );

        if ( direct == FLA_FORWARD )
            FLA_Sort_bsvd_ext_f_opc( m_s, s_p, inc_s,
                                     m_U, U_p, rs_U, cs_U,
                                     m_V, V_p, rs_V, cs_V,
                                     n_C, C_p, rs_C, cs_C );
        else // if ( direct == FLA_BACKWARD )
            FLA_Sort_bsvd_ext_b_opc( m_s, s_p, inc_s,
                                     m_U, U_p, rs_U, cs_U,
                                     m_V, V_p, rs_V, cs_V,
                                     n_C, C_p, rs_C, cs_C );
        break;
    }
    case FLA_DOUBLE_COMPLEX:
    {
        double*   s_p = ( double*   ) FLA_DOUBLE_PTR( s );
        dcomplex* U_p = ( apply_U == TRUE ? ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( U ) : NULL );
        dcomplex* V_p = ( apply_V == TRUE ? ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( V ) : NULL );
        dcomplex* C_p = ( apply_C == TRUE ? ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( C ) : NULL );

        if ( direct == FLA_FORWARD )
            FLA_Sort_bsvd_ext_f_opz( m_s, s_p, inc_s,
                                     m_U, U_p, rs_U, cs_U,
                                     m_V, V_p, rs_V, cs_V,
                                     n_C, C_p, rs_C, cs_C );
        else // if ( direct == FLA_BACKWARD )
            FLA_Sort_bsvd_ext_b_opz( m_s, s_p, inc_s,
                                     m_U, U_p, rs_U, cs_U,
                                     m_V, V_p, rs_V, cs_V,
                                     n_C, C_p, rs_C, cs_C );
        break;
    }
    }
    return FLA_SUCCESS;
}

// single
FLA_Error FLA_Sort_bsvd_ext_f_ops( int m_s, float* s, int inc_s,
                                   int m_U, float* U, int rs_U, int cs_U,
                                   int m_V, float* V, int rs_V, int cs_V,
                                   int n_C, float* C, int rs_C, int cs_C )
{
    int    i, ii, j, k;
    float  p;
    FLA_SORT_BSVD_EXT_BODY( FORWARD, bl1_sswapv );
    return FLA_SUCCESS;
}
FLA_Error FLA_Sort_bsvd_ext_b_ops( int m_s, float* s, int inc_s,
                                   int m_U, float* U, int rs_U, int cs_U,
                                   int m_V, float* V, int rs_V, int cs_V,
                                   int n_C, float* C, int rs_C, int cs_C )
{
    int    i, ii, j, k;
    float  p;
    FLA_SORT_BSVD_EXT_BODY( BACKWARD, bl1_sswapv );
    return FLA_SUCCESS;
}

// double
FLA_Error FLA_Sort_bsvd_ext_f_opd( int m_s, double* s, int inc_s,
                                   int m_U, double* U, int rs_U, int cs_U,
                                   int m_V, double* V, int rs_V, int cs_V,
                                   int n_C, double* C, int rs_C, int cs_C )
{
    int    i, ii, j, k;
    float  p;
    FLA_SORT_BSVD_EXT_BODY( FORWARD, bl1_dswapv );
    return FLA_SUCCESS;
}
FLA_Error FLA_Sort_bsvd_ext_b_opd( int m_s, double* s, int inc_s,
                                   int m_U, double* U, int rs_U, int cs_U,
                                   int m_V, double* V, int rs_V, int cs_V,
                                   int n_C, double* C, int rs_C, int cs_C )
{
    int    i, ii, j, k;
    double p;
    FLA_SORT_BSVD_EXT_BODY( BACKWARD, bl1_dswapv );
    return FLA_SUCCESS;
}

// scomplex
FLA_Error FLA_Sort_bsvd_ext_f_opc( int m_s, float*    s, int inc_s,
                                   int m_U, scomplex* U, int rs_U, int cs_U,
                                   int m_V, scomplex* V, int rs_V, int cs_V,
                                   int n_C, scomplex* C, int rs_C, int cs_C )
{
    int    i, ii, j, k;
    float  p;
    FLA_SORT_BSVD_EXT_BODY( FORWARD, bl1_cswapv );
    return FLA_SUCCESS;
}
FLA_Error FLA_Sort_bsvd_ext_b_opc( int m_s, float*    s, int inc_s,
                                   int m_U, scomplex* U, int rs_U, int cs_U,
                                   int m_V, scomplex* V, int rs_V, int cs_V,
                                   int n_C, scomplex* C, int rs_C, int cs_C )
{
    int    i, ii, j, k;
    float  p;
    FLA_SORT_BSVD_EXT_BODY( BACKWARD, bl1_cswapv );
    return FLA_SUCCESS;
}

// dcomplex
FLA_Error FLA_Sort_bsvd_ext_f_opz( int m_s, double*   s, int inc_s,
                                   int m_U, dcomplex* U, int rs_U, int cs_U,
                                   int m_V, dcomplex* V, int rs_V, int cs_V,
                                   int n_C, dcomplex* C, int rs_C, int cs_C )
{
    int    i, ii, j, k;
    double p;
    FLA_SORT_BSVD_EXT_BODY( FORWARD, bl1_zswapv );
    return FLA_SUCCESS;
}
FLA_Error FLA_Sort_bsvd_ext_b_opz( int m_s, double*   s, int inc_s,
                                   int m_U, dcomplex* U, int rs_U, int cs_U,
                                   int m_V, dcomplex* V, int rs_V, int cs_V,
                                   int n_C, dcomplex* C, int rs_C, int cs_C )
{
    int    i, ii, j, k;
    double p;
    FLA_SORT_BSVD_EXT_BODY( BACKWARD, bl1_zswapv );
    return FLA_SUCCESS;
}
