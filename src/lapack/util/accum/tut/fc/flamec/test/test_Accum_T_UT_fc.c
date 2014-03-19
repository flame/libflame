/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_UNBLOCKED 1
#define FLA_ALG_UNB_OPT   2

void time_Accum_T_UT_fc(
                 int variant, int type, int n_repeats, int m, int n,
                 FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj W, FLA_Obj b, FLA_Obj b_ref,
                 double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    m_input, n_input,
    m, n, min_m_n,
    p_first, p_last, p_inc,
    p,
    variant,
    n_repeats,
    i,
    n_variants = 1;
  
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char n_dim_desc[14];
  char m_dim_tag[10];
  char n_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, t, T, W, b, b_ref;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d %d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );


  fprintf( stdout, "\nclear all;\n\n" );


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
  if     ( n_input >  0 ) {
    sprintf( n_dim_desc, "n = %d", n_input );
    sprintf( n_dim_tag,  "n%dc", n_input);
  }
  else if( n_input <  -1 ) {
    sprintf( n_dim_desc, "n = p/%d", -n_input );
    sprintf( n_dim_tag,  "n%dp", -n_input );
  }
  else if( n_input == -1 ) {
    sprintf( n_dim_desc, "n = p" );
    sprintf( n_dim_tag,  "n%dp", 1 );
  }


  //datatype = FLA_FLOAT;
  //datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  datatype = FLA_DOUBLE_COMPLEX;

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;
    n = n_input;

    if( m < 0 ) m = p / abs(m_input);
    if( n < 0 ) n = p / abs(n_input);

    min_m_n = min( m, n );

    for ( variant = 0; variant < n_variants; variant++ ){
      
/*
      FLA_Obj_create( datatype, m, n, 0, 0, &A );
      FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &t );
      FLA_Obj_create( datatype, n, n, 0, 0, &T );
      FLA_Obj_create( datatype, m, 1, 0, 0, &W );
      FLA_Obj_create( datatype, m, 1, 0, 0, &b );
      FLA_Obj_create( datatype, m, 1, 0, 0, &b_ref );
*/
      FLA_Obj_create( datatype, m, n, n, 1, &A );
      FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &t );
      FLA_Obj_create( datatype, n, n, n, 1, &T );
      FLA_Obj_create( datatype, m, 1, 1, 1, &W );
      FLA_Obj_create( datatype, m, 1, 1, 1, &b );
      FLA_Obj_create( datatype, m, 1, 1, 1, &b_ref );

      FLA_Random_matrix( A );
      FLA_Random_matrix( b );

      fprintf( stdout, "data_qr_ut( %d, 1:7 ) = [ %d  ", i, p );
      fflush( stdout );

      time_Accum_T_UT_fc( variant, FLA_ALG_REFERENCE, n_repeats, m, n,
	                      A, t, T, W, b, b_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Accum_T_UT_fc( variant, FLA_ALG_UNBLOCKED, n_repeats, m, n,
	                      A, t, T, W, b, b_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Accum_T_UT_fc( variant, FLA_ALG_UNB_OPT, n_repeats, m, n,
	                      A, t, T, W, b, b_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &A );
      FLA_Obj_free( &t );
      FLA_Obj_free( &T );
      FLA_Obj_free( &W );
      FLA_Obj_free( &b );
      FLA_Obj_free( &b_ref );
    }

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_variants; i++ ) {
    fprintf( stdout, "plot( data_qr_ut( :,1 ), data_qr_ut( :, 2 ), '%c:%c' ); \n",
            colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_qr_ut( :,1 ), data_qr_ut( :, 4 ), '%c-.%c' ); \n",
            colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_variants; i++ )
    fprintf( stdout, "'ref\\_qr\\_ut', 'fla\\_qr\\_ut', ... \n" );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME Accum_T_UT_fc front-end performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc qr_ut_front_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

