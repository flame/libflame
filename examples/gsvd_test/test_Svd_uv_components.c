
#include "FLAME.h"

#define FLA_ALG_REFERENCE     0
#define FLA_ALG_UNBLOCKED     1
#define FLA_ALG_UNB_OPT       2
#define FLA_ALG_BLOCKED       3


void time_Svd_uv_components(
               int variant, int type, int n_repeats, int m, int n, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
               double* dtime, double* diff1, double* diff2, double* gflops,
               double* dtime_bred, double* gflops_bred,
               double* dtime_bsvd, double* gflops_bsvd,
               double* dtime_appq, double* gflops_appq,
               double* dtime_qrfa, double* gflops_qrfa,
               double* dtime_gemm, double* gflops_gemm, int* k_perf );


int main(int argc, char *argv[])
{
  int 
    m_input, n_input,
    m, n,
    p_first, p_last, p_inc,
    p,
    min_m_n,
    k_accum,
    b_alg,
    n_iter_max,
    dist_type,
    variant,
    n_repeats,
    i,
    k_perf,
    first_var = 1,
    last_var  = 2;
  
  double
    dtime,
    gflops,
    diff1,
    diff2;
  double
    dtime_tred,
    dtime_tevd,
    dtime_appq,
    dtime_qrfa,
    dtime_gemm,
    gflops_tred,
    gflops_tevd,
    gflops_appq,
    gflops_qrfa,
    gflops_gemm;

  FLA_Datatype datatype, dt_real;

  FLA_Obj
    A, s, U, V, alpha, shift, nc;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter n_iter_max:", '%' );
  scanf( "%d", &n_iter_max );
  fprintf( stdout, "%c %d\n", '%', n_iter_max );

  fprintf( stdout, "%c enter number of Givens rotations to accumulate:", '%' );
  scanf( "%d", &k_accum );
  fprintf( stdout, "%c %d\n", '%', k_accum );

  fprintf( stdout, "%c enter blocking size:", '%' );
  scanf( "%d", &b_alg );
  fprintf( stdout, "%c %d\n", '%', b_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d %d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );

  fprintf( stdout, "%c enter distribution type (0=u, 1=i, 2=g, 3=l, 4=r, 5=c): ", '%' );
  scanf( "%d", &dist_type );
  fprintf( stdout, "%c %d\n", '%', dist_type );


  fprintf( stdout, "\n" );



  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;
    n = n_input;

    if( m < 0 ) m = p / abs(m_input);
    if( n < 0 ) n = p / abs(n_input);

    min_m_n = min( m, n );

    //datatype = FLA_FLOAT;
    //datatype = FLA_DOUBLE;
    //datatype = FLA_COMPLEX;
    datatype = FLA_DOUBLE_COMPLEX;

    FLA_Obj_create( datatype, m,       n, 0, 0, &A );
    FLA_Obj_create( datatype, m,       m, 0, 0, &U );
    FLA_Obj_create( datatype, n,       n, 0, 0, &V );

    dt_real = FLA_Obj_datatype_proj_to_real( A );

    FLA_Obj_create( dt_real,  min_m_n, 1, 0, 0, &s );
    FLA_Obj_create( dt_real,  1,       1, 0, 0, &alpha );
    FLA_Obj_create( dt_real,  1,       1, 0, 0, &shift );
    FLA_Obj_create( FLA_INT,  1,       1, 0, 0, &nc );

    FLA_Random_unitary_matrix( U );
    FLA_Random_unitary_matrix( V );

    if ( dist_type == 0 )
    {
      // Linear
      *FLA_DOUBLE_PTR( shift ) = 0.0;
      *FLA_DOUBLE_PTR( alpha ) = 1.0;
      fprintf( stdout, "%c using linear dist.\n", '%' );
      fprintf( stdout, "%c delta = %9.3e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_linear_dist( shift, alpha, s );
    }
    else if ( dist_type == 1 )
    {
      // Inverse
      *FLA_DOUBLE_PTR( alpha ) = 1.0;
      fprintf( stdout, "%c using inverse dist.\n", '%' );
      fprintf( stdout, "%c alpha = %9.3e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_inverse_dist( alpha, s );
    }
    else if ( dist_type == 2 )
    {
      // Geometric
      *FLA_DOUBLE_PTR( alpha ) = 1.0 / (double)min_m_n;
      fprintf( stdout, "%c using geometric dist.\n", '%' );
      fprintf( stdout, "%c alpha = %10.4e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_geometric_dist( alpha, s );
    }
    else if ( dist_type == 3 )
    {
      // Logarithmic
      *FLA_DOUBLE_PTR( alpha ) = 1.20;
      fprintf( stdout, "%c using logarithmic dist.\n", '%' );
      fprintf( stdout, "%c alpha = %9.3e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_logarithmic_dist( alpha, s );
    }
    else if ( dist_type == 4 )
    {
      // Random
      *FLA_DOUBLE_PTR( shift ) = 0.0;
      *FLA_DOUBLE_PTR( alpha ) = (double)min_m_n;
      fprintf( stdout, "%c using random dist.\n", '%' );
      fprintf( stdout, "%c shift = %13.8e\n", '%', *FLA_DOUBLE_PTR( shift ) );
      fprintf( stdout, "%c alpha = %9.3e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_random_dist( shift, alpha, s );
    }
    else if ( dist_type == 5 )
    {
      // Cluster
      *FLA_INT_PTR( nc )       = 10;
      *FLA_DOUBLE_PTR( alpha ) = 1.0e-9;
      fprintf( stdout, "%c using cluster dist.\n", '%' );
      fprintf( stdout, "%c num clusters  = %d\n", '%', *FLA_INT_PTR( nc ) );
      fprintf( stdout, "%c cluster width = %9.3e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_cluster_dist( nc, alpha, s );
    }

    {
      FLA_Obj UL, UR;
      FLA_Obj VL, VR;

      FLA_Part_1x2( U,   &UL, &UR,   min_m_n, FLA_LEFT );
      FLA_Part_1x2( V,   &VL, &VR,   min_m_n, FLA_LEFT );

      FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, s, UL );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                FLA_ONE, UL, VL, FLA_ZERO, A );
    }

    FLA_Set( FLA_ZERO, s );
    FLA_Set( FLA_ZERO, U );
    FLA_Set( FLA_ZERO, V );

    fprintf( stdout, "%c                               total----------                       reduction------    bi svd---------   form/app QUV---   QR fact--------   gemm-----------\n", '%' );
    fprintf( stdout, "%c                         %4s  gflops dtime      resid    |I-QQ'|    gflops dtime       gflops dtime      gflops dtime      gflops dtime      gflops dtime     niter\n", '%', "p" );

    time_Svd_uv_components( -2, FLA_ALG_REFERENCE, n_repeats,
                            m, n, n_iter_max, k_accum, b_alg,
                            A, s, U, V, &dtime, &diff1, &diff2, &gflops,
                            &dtime_tred, &gflops_tred,
                            &dtime_tevd, &gflops_tevd,
                            &dtime_appq, &gflops_appq,
                            &dtime_qrfa, &gflops_qrfa,
                            &dtime_gemm, &gflops_gemm, &k_perf );
    if ( dtime_tred == 1.0 ) gflops_tred = gflops_tevd = gflops_appq = gflops_qrfa = gflops_gemm = 
                              dtime_tred =  dtime_tevd =  dtime_appq =  dtime_qrfa =  dtime_gemm = 0.0;

    fprintf( stdout, "data_refq( %2d, 1:16 ) = [ %4d %6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e  %6.3lf %9.2e  %6.3lf %9.2e %5d ];\n",
             i, p, gflops, dtime, diff1, diff2,
             gflops_tred, dtime_tred,
             gflops_tevd, dtime_tevd,
             gflops_appq, dtime_appq,
             gflops_qrfa, dtime_qrfa,
             gflops_gemm, dtime_gemm, k_perf );
    fflush( stdout );

    time_Svd_uv_components( -3, FLA_ALG_REFERENCE, n_repeats,
                            m, n, n_iter_max, k_accum, b_alg,
                            A, s, U, V, &dtime, &diff1, &diff2, &gflops,
                            &dtime_tred, &gflops_tred,
                            &dtime_tevd, &gflops_tevd,
                            &dtime_appq, &gflops_appq,
                            &dtime_qrfa, &gflops_qrfa,
                            &dtime_gemm, &gflops_gemm, &k_perf );
    if ( dtime_tred == 1.0 ) gflops_tred = gflops_tevd = gflops_appq = gflops_qrfa = gflops_gemm = 
                              dtime_tred =  dtime_tevd =  dtime_appq =  dtime_qrfa =  dtime_gemm = 0.0;

    fprintf( stdout, "data_refd( %2d, 1:16 ) = [ %4d %6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e  %6.3lf %9.2e  %6.3lf %9.2e %5d ];\n",
             i, p, gflops, dtime, diff1, diff2,
             gflops_tred, dtime_tred,
             gflops_tevd, dtime_tevd,
             gflops_appq, dtime_appq,
             gflops_qrfa, dtime_qrfa,
             gflops_gemm, dtime_gemm, k_perf );
    fflush( stdout );

    time_Svd_uv_components( 0, FLA_ALG_REFERENCE, n_repeats,
                            m, n, n_iter_max, k_accum, b_alg,
                            A, s, U, V, &dtime, &diff1, &diff2, &gflops,
                            &dtime_tred, &gflops_tred,
                            &dtime_tevd, &gflops_tevd,
                            &dtime_appq, &gflops_appq,
                            &dtime_qrfa, &gflops_qrfa,
                            &dtime_gemm, &gflops_gemm, &k_perf );
    if ( dtime_tred == 1.0 ) gflops_tred = gflops_tevd = gflops_appq = gflops_qrfa = gflops_gemm = 
                              dtime_tred =  dtime_tevd =  dtime_appq =  dtime_qrfa =  dtime_gemm = 0.0;

    fprintf( stdout, "data_REFq( %2d, 1:16 ) = [ %4d %6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e  %6.3lf %9.2e  %6.3lf %9.2e %5d ];\n",
             i, p, gflops, dtime, diff1, diff2,
             gflops_tred, dtime_tred,
             gflops_tevd, dtime_tevd,
             gflops_appq, dtime_appq,
             gflops_qrfa, dtime_qrfa,
             gflops_gemm, dtime_gemm, k_perf );
    fflush( stdout );

    time_Svd_uv_components( -1, FLA_ALG_REFERENCE, n_repeats,
                            m, n, n_iter_max, k_accum, b_alg,
                            A, s, U, V, &dtime, &diff1, &diff2, &gflops,
                            &dtime_tred, &gflops_tred,
                            &dtime_tevd, &gflops_tevd,
                            &dtime_appq, &gflops_appq,
                            &dtime_qrfa, &gflops_qrfa,
                            &dtime_gemm, &gflops_gemm, &k_perf );
    if ( dtime_tred == 1.0 ) gflops_tred = gflops_tevd = gflops_appq = gflops_qrfa = gflops_gemm = 
                              dtime_tred =  dtime_tevd =  dtime_appq =  dtime_qrfa =  dtime_gemm = 0.0;

    fprintf( stdout, "data_REFd( %2d, 1:16 ) = [ %4d %6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e  %6.3lf %9.2e  %6.3lf %9.2e %5d ];\n",
             i, p, gflops, dtime, diff1, diff2,
             gflops_tred, dtime_tred,
             gflops_tevd, dtime_tevd,
             gflops_appq, dtime_appq,
             gflops_qrfa, dtime_qrfa,
             gflops_gemm, dtime_gemm, k_perf );
    fflush( stdout );


    for ( variant = first_var; variant <= last_var; variant++ ){
      
      fprintf( stdout, "data_var%d( %2d, 1:16 ) = [ %4d ", variant, i, p );
      fflush( stdout );

      time_Svd_uv_components( variant, FLA_ALG_UNBLOCKED, n_repeats,
                              m, n, n_iter_max, k_accum, b_alg,
                              A, s, U, V, &dtime, &diff1, &diff2, &gflops,
                              &dtime_tred, &gflops_tred,
                              &dtime_tevd, &gflops_tevd,
                              &dtime_appq, &gflops_appq,
                              &dtime_qrfa, &gflops_qrfa,
                              &dtime_gemm, &gflops_gemm, &k_perf );

      fprintf( stdout, "%6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e  %6.3lf %9.2e  %6.3lf %9.2e %5d ",
               gflops, dtime, diff1, diff2,
               gflops_tred, dtime_tred,
               gflops_tevd, dtime_tevd,
               gflops_appq, dtime_appq,
               gflops_qrfa, dtime_qrfa,
               gflops_gemm, dtime_gemm, k_perf );

      fprintf( stdout, "];\n" );
      fflush( stdout );
    }

    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &U );
    FLA_Obj_free( &V );
    FLA_Obj_free( &s );
    FLA_Obj_free( &alpha );
    FLA_Obj_free( &shift );
    FLA_Obj_free( &nc );
  }

  FLA_Finalize( );

  return 0;
}

