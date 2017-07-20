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


void time_Tridiag_UT_l(
               int variant, int type, int n_repeats, int m, int nb_alg,
               FLA_Obj A, FLA_Obj U, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj TT, FLA_Obj t,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    m_input,
    m,
    p_first, p_last, p_inc,
    p,
    nb_alg,
    variant,
    n_repeats,
    i,
    datatype,
    n_variants = 3;
  
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
    A, t, U, Y, Z, T, TT;
  

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


  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;

    if( m < 0 ) m = p / f2c_abs(m_input);

    //datatype = FLA_FLOAT;
    datatype = FLA_DOUBLE;
    //datatype = FLA_COMPLEX;
    //datatype = FLA_DOUBLE_COMPLEX;

    FLA_Obj_create( datatype, m,      m, 0, 0, &A );
    FLA_Obj_create( datatype, m,      m, 0, 0, &U );
    FLA_Obj_create( datatype, m,      m, 0, 0, &Y );
    FLA_Obj_create( datatype, m,      m, 0, 0, &Z );
    FLA_Obj_create( datatype, m,      m, 0, 0, &T );
    FLA_Obj_create( datatype, nb_alg, m, 0, 0, &TT );
    FLA_Obj_create( datatype, m,      1, 0, 0, &t );

    //FLA_Random_herm_matrix( FLA_LOWER_TRIANGULAR, A );
    FLA_Random_spd_matrix( FLA_LOWER_TRIANGULAR, A );
    FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A );


    time_Tridiag_UT_l( 0, FLA_ALG_REFERENCE, n_repeats, m, nb_alg,
                       A, U, Y, Z, T, TT, t, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REFb( %d, 1:3 ) = [ %d %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );


    for ( variant = 1; variant <= n_variants; variant++ ){
      
      fprintf( stdout, "data_var%d( %d, 1:11 ) = [ %d ", variant, i, p );
      fflush( stdout );

      time_Tridiag_UT_l( variant, FLA_ALG_UNBLOCKED, n_repeats, m, nb_alg,
                         A, U, Y, Z, T, TT, t, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Tridiag_UT_l( variant, FLA_ALG_UNB_OPT, n_repeats, m, nb_alg,
                         A, U, Y, Z, T, TT, t, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Tridiag_UT_l( variant, FLA_ALG_UNB_OPT_FUSED, n_repeats, m, nb_alg,
                         A, U, Y, Z, T, TT, t, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Tridiag_UT_l( variant, FLA_ALG_BLOCKED, n_repeats, m, nb_alg,
                         A, U, Y, Z, T, TT, t, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Tridiag_UT_l( variant, FLA_ALG_BLOCKED_FUSED, n_repeats, m, nb_alg,
                         A, U, Y, Z, T, TT, t, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, "];\n" );
      fflush( stdout );
    }

    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &U );
    FLA_Obj_free( &Y );
    FLA_Obj_free( &Z );
    FLA_Obj_free( &T );
    FLA_Obj_free( &TT );
    FLA_Obj_free( &t );
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
  fprintf( stdout, "title( 'FLAME Tridiag performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc tridiag_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

