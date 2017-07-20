/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE     0
#define FLA_ALG_UNBLOCKED     1
#define FLA_ALG_UNB_OPT       2
#define FLA_ALG_BLOCKED       3

void fill_eigenvalues( FLA_Obj l );

void time_Hevd_ln(
               int variant, int type, int n_repeats, int m, int b_alg,
               FLA_Obj A, FLA_Obj l,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    m_input,
    m,
    p_first, p_last, p_inc,
    p,
    b_alg,
    variant,
    n_repeats,
    i,
    datatype,
    n_variants = 1;
  
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];

  double max_gflops=6.0;
  double safemin;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, l, Q, T, W;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c Enter blocking size:", '%' );
  scanf( "%d", &b_alg );
  fprintf( stdout, "%c %d\n", '%', b_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m (-1 means bind to problem size): ", '%' );
  scanf( "%d", &m_input );
  fprintf( stdout, "%c %d\n", '%', m_input );


  fprintf( stdout, "\n" );


  if     ( m_input >  0 ) {
    sprintf( m_dim_desc, "m = %d", m_input );
    sprintf( m_dim_tag,  "m%dc", m_input);
  }
  else if( m_input <  -1 ) {
    sprintf( m_dim_desc, "m = p/%d", -m_input );
    sprintf( m_dim_tag,  "m%dp", -m_input );
  }
  else if( m_input == -1 ) {
    sprintf( m_dim_desc, "m = p" );
    sprintf( m_dim_tag,  "m%dp", 1 );
  }
/*
char ch = 's';
safemin = dlamch_( &ch );
printf( "safemin = %23.15e\n", safemin );
ch = 'e';
double eps = dlamch_( &ch );
printf( "eps dla = %23.15e\n", eps );
printf( "eps fla = %23.15e\n", FLA_EPSILON_D );
*/

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;

    if( m < 0 ) m = p / f2c_abs(m_input);

    //datatype = FLA_FLOAT;
    //datatype = FLA_DOUBLE;
    //datatype = FLA_COMPLEX;
    datatype = FLA_DOUBLE_COMPLEX;

    FLA_Obj_create( datatype,                           m,      m, 0, 0, &A );
    FLA_Obj_create( datatype,                           m,      m, 0, 0, &Q );
    FLA_Obj_create( datatype,                           32,     m, 0, 0, &T );
    FLA_Obj_create( datatype,                           32,     m, 0, 0, &W );
    FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), m,      1, 0, 0, &l );

    //FLA_Random_herm_matrix( FLA_LOWER_TRIANGULAR, A );
    //FLA_Random_spd_matrix( FLA_LOWER_TRIANGULAR, A );

    FLA_Random_matrix( A );
    FLA_Obj_set_to_identity( Q );
    FLA_QR_UT( A, T );
    FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, T, W, Q );
    fill_eigenvalues( l );
    //FLA_Obj_show( "eig", l, "%9.2e ", "" );
    FLA_Apply_diag_matrix( FLA_LEFT, FLA_NO_CONJUGATE, l, Q );
    FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, A, T, W, Q );
    FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, Q );
    FLA_Copy( Q, A );

    time_Hevd_ln( 0, FLA_ALG_REFERENCE, n_repeats, m, b_alg,
                  A, l, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REFs( %d, 1:2 ) = [ %d %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );

    time_Hevd_ln( -1, FLA_ALG_REFERENCE, n_repeats, m, b_alg,
                  A, l, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REFd( %d, 1:2 ) = [ %d %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );


    for ( variant = 1; variant <= n_variants; variant++ ){
      
      fprintf( stdout, "data_var%d( %d, 1:9 ) = [ %d ", variant, i, p );
      fflush( stdout );

      time_Hevd_ln( variant, FLA_ALG_UNBLOCKED, n_repeats, m, b_alg,
                    A, l, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      //time_Hevd_ln( variant, FLA_ALG_UNB_OPT, n_repeats, m, b_alg,
      //              A, l, &dtime, &diff, &gflops );

      //fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      //fflush( stdout );

      fprintf( stdout, "];\n" );
      fflush( stdout );
    }

    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &T );
    FLA_Obj_free( &W );
    FLA_Obj_free( &Q );
    FLA_Obj_free( &l );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), '-' ); \n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 1; i <= n_variants; i++ ) {
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 2 ), '%c:%c' ); \n",
            i, i, colors[ i-1 ], ticks[ i-1 ] );
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 4 ), '%c-.%c' ); \n",
            i, i, colors[ i-1 ], ticks[ i-1 ] );
  }

  fprintf( stdout, "legend( ... \n" );
  fprintf( stdout, "'Reference', ... \n" );

  for ( i = 1; i < n_variants; i++ )
    fprintf( stdout, "'unb\\_var%d', 'blk\\_var%d', ... \n", i, i );
  fprintf( stdout, "'unb\\_var%d', 'blk\\_var%d' ); \n", i, i );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME Hevd_ln performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc tridiag_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

void fill_eigenvalues( FLA_Obj l )
{
  FLA_Obj lT,              l0,
          lB,              lambda1,
                           l2;
  FLA_Obj alpha;

  FLA_Obj_create( FLA_Obj_datatype( l ), 1, 1, 0, 0, &alpha );
  FLA_Copy( FLA_ONE, alpha );

  FLA_Part_2x1( l,    &lT, 
                      &lB,            0, FLA_TOP );

  while ( FLA_Obj_length( lT ) < FLA_Obj_length( l ) ){

    FLA_Repart_2x1_to_3x1( lT,                &l0, 
                        /* ** */            /* ******* */
                                              &lambda1, 
                           lB,                &l2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Copy( alpha, lambda1 );
    FLA_Mult_add( FLA_ONE, FLA_ONE, alpha );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &lT,                l0, 
                                                  lambda1, 
                            /* ** */           /* ******* */
                              &lB,                l2,     FLA_TOP );

  }

  FLA_Obj_free( &alpha );
}
