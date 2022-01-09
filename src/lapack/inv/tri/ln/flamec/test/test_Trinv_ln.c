/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define N_VARIANTS        4

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_BLOCKED   1
#define FLA_ALG_UNBLOCKED 2
#define FLA_ALG_UNB_OPT   3


void time_Trinv_ln(
                  integer variant, integer type, integer n_repeats, integer m, integer nb_alg,
                  FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj norm,
                  double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  integer 
    datatype,
    m_input,
    m,
    p_first, p_last, p_inc,
    p,
    nb_alg,
    variant,
    n_repeats,
    i, j,
    n_variants = N_VARIANTS;
  
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];

  double max_gflops=6.0;
  

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, b, b_orig, norm;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c Enter blocking size:", '%' );
  scanf( "%d", &nb_alg );
  fprintf( stdout, "%c %d\n", '%', nb_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m (-1 means bind to problem size): ", '%' );
  scanf( "%d", &m_input );
  fprintf( stdout, "%c %d\n", '%', m_input );


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

  //datatype = FLA_FLOAT;
  //datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  datatype = FLA_DOUBLE_COMPLEX;

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;

    if( m < 0 ) m = p / f2c_abs(m_input);

    FLA_Obj_create( datatype, m, m, 0, 0, &A );
    FLA_Obj_create( datatype, m, 1, 0, 0, &b );
    FLA_Obj_create( datatype, m, 1, 0, 0, &b_orig );
/*
    FLA_Obj_create( datatype, m, m, m, 1, &A );
    FLA_Obj_create( datatype, m, 1, 1, 1, &b );
    FLA_Obj_create( datatype, m, 1, 1, 1, &b_orig );
*/

    if ( FLA_Obj_is_single_precision( A ) )
      FLA_Obj_create( FLA_FLOAT, 1, 1, 0, 0, &norm );
    else
      FLA_Obj_create( FLA_DOUBLE, 1, 1, 0, 0, &norm );

    FLA_Random_tri_matrix( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
    FLA_Random_matrix( b );
    FLA_Copy_external( b, b_orig );

/*
    time_Trinv_ln( 0, FLA_ALG_REFERENCE, n_repeats, m, nb_alg,
                   A, b, b_orig, norm, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REF( %d, 1:2 ) = [ %d  %6.3lf ]; \n", i, p, gflops );
    fflush( stdout );
*/

    for ( variant = 1; variant <= n_variants; variant++ ){
      
      fprintf( stdout, "data_var%d( %d, 1:7 ) = [ %d  ", variant, i, p );
      fflush( stdout );

      time_Trinv_ln( variant, FLA_ALG_UNBLOCKED, n_repeats, m, nb_alg,
                     A, b, b_orig, norm, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Trinv_ln( variant, FLA_ALG_UNB_OPT, n_repeats, m, nb_alg,
                     A, b, b_orig, norm, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Trinv_ln( variant, FLA_ALG_BLOCKED, n_repeats, m, nb_alg,
                     A, b, b_orig, norm, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );
    }

    FLA_Obj_free( &A );
    FLA_Obj_free( &b );
    FLA_Obj_free( &b_orig );
    FLA_Obj_free( &norm );

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), '-' ); \n" );

  for ( i = 1; i <= n_variants; i++ ){
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 2 ), '%c:%c' ); \n", 
             variant, variant, colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );
  fprintf( stdout, "'Reference', ... \n" );

  for ( i = 1; i <= n_variants; i++ )
    fprintf( stdout, "'FLAME var%d', ... \n", i );

  fprintf( stdout, "'Location', 'SouthWest' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME trinv\\_u performance (%s)' );\n", 
           m_dim_desc );
  fprintf( stdout, "print -depsc trinv_l_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );
}

