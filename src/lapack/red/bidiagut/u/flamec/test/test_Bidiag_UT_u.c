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
#define FLA_ALG_UNB_OPT_FUSED 3
#define FLA_ALG_BLOCKED       4
#define FLA_ALG_BLOCKED_FUSED 5


void time_Bidiag_UT_u(
               integer variant, integer type, integer n_repeats, integer m, integer n, integer nb_alg,
               FLA_Obj A, FLA_Obj TU, FLA_Obj TV, FLA_Obj TTU, FLA_Obj TTV, FLA_Obj tu, FLA_Obj tv,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  integer 
    m_input, n_input,
    m, n,
    min_m_n,
    p_first, p_last, p_inc,
    p,
    nb_alg,
    variant,
    n_repeats,
    i,
    datatype,
    n_variants = 5;
  
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];
  char n_dim_desc[14];
  char n_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, TU, TV, TTU, TTV, tu, tv;
  

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

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d %d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );


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

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;
    n = n_input;

    if( m < 0 ) m = p / f2c_abs(m_input);
    if( n < 0 ) n = p / f2c_abs(n_input);

    min_m_n = fla_min( m, n );

    //datatype = FLA_FLOAT;
    //datatype = FLA_DOUBLE;
    //datatype = FLA_COMPLEX;
    datatype = FLA_DOUBLE_COMPLEX;

    FLA_Obj_create( datatype, m,       n,       0, 0, &A );
    FLA_Obj_create( datatype, min_m_n, min_m_n, 0, 0, &TU );
    FLA_Obj_create( datatype, min_m_n, min_m_n, 0, 0, &TV );
    FLA_Obj_create( datatype, nb_alg,  min_m_n, 0, 0, &TTU );
    FLA_Obj_create( datatype, nb_alg,  min_m_n, 0, 0, &TTV );
    FLA_Obj_create( datatype, min_m_n,       1, 0, 0, &tu );
    FLA_Obj_create( datatype, min_m_n,       1, 0, 0, &tv );

    FLA_Random_matrix( A );

    time_Bidiag_UT_u( 0, FLA_ALG_REFERENCE, n_repeats, m, n, nb_alg,
                      A, TU, TV, TTU, TTV, tu, tv, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REFb( %d, 1:3 ) = [ %d %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );

    for ( variant = 1; variant <= n_variants; variant++ ){
      
      fprintf( stdout, "data_var%d( %d, 1:11 ) = [ %d ", variant, i, p );
      fflush( stdout );

      time_Bidiag_UT_u( variant, FLA_ALG_UNBLOCKED, n_repeats, m, n, nb_alg,
                        A, TU, TV, TTU, TTV, tu, tv, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Bidiag_UT_u( variant, FLA_ALG_UNB_OPT, n_repeats, m, n, nb_alg,
                        A, TU, TV, TTU, TTV, tu, tv, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Bidiag_UT_u( variant, FLA_ALG_UNB_OPT_FUSED, n_repeats, m, n, nb_alg,
                        A, TU, TV, TTU, TTV, tu, tv, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Bidiag_UT_u( variant, FLA_ALG_BLOCKED, n_repeats, m, n, nb_alg,
                        A, TU, TV, TTU, TTV, tu, tv, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Bidiag_UT_u( variant, FLA_ALG_BLOCKED_FUSED, n_repeats, m, n, nb_alg,
                        A, TU, TV, TTU, TTV, tu, tv, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, "];\n" );
      fflush( stdout );
    }

    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &TU );
    FLA_Obj_free( &TV );
    FLA_Obj_free( &TTU );
    FLA_Obj_free( &TTV );
    FLA_Obj_free( &tu );
    FLA_Obj_free( &tv );
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
  fprintf( stdout, "title( 'FLAME Bidiag performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc bidiag_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

