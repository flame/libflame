/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define N_PARAM_COMBOS    2

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

char* pc_str[N_PARAM_COMBOS] = { "l", "u" };

void time_Tridiag_UT(
                int param_combo, int type, int n_repeats, int m,
                FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj W,
                double *dtime, double *diff, double *gflops );

int main(int argc, char *argv[])
{
  int 
    m_input,
    m,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    //n_param_combos = N_PARAM_COMBOS;
    n_param_combos = 1;

  FLA_Datatype datatype;

  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];

  double max_gflops=6.0;

  dim_t b_alg;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, t, T, W;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c Enter blocking size:", '%' );
  scanf( "%u", &b_alg );
  fprintf( stdout, "%c %u\n", '%', b_alg );

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

    if( m < 0 ) m = p * f2c_abs(m_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){
      
      FLA_Obj_create( datatype, m,     m, 0, 0, &A );
      FLA_Obj_create( datatype, m-1,   1, 0, 0, &t );

      if ( b_alg == 0 ) FLA_Tridiag_UT_create_T( A, &T );
      else              FLA_Obj_create( datatype, b_alg, m, 0, 0, &T );

      FLA_Obj_create( datatype, FLA_Obj_length( T ), m, 0, 0, &W );

      if ( pc_str[param_combo][0] == 'l' )
      {
        FLA_Random_spd_matrix( FLA_LOWER_TRIANGULAR, A );
        FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A );
      }
      else
      {
        FLA_Random_spd_matrix( FLA_UPPER_TRIANGULAR, A );
        FLA_Hermitianize( FLA_UPPER_TRIANGULAR, A );
      }


      fprintf( stdout, "data_tridiag_ut( %d, 1:3 ) = [ %d  ", i, p );
      fflush( stdout );

      time_Tridiag_UT( param_combo, FLA_ALG_REFERENCE, n_repeats, m,
                       A, t, T, W, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      time_Tridiag_UT( param_combo, FLA_ALG_FRONT, n_repeats, m,
                       A, t, T, W, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &A );
      FLA_Obj_free( &t );
      FLA_Obj_free( &T );
      FLA_Obj_free( &W );
    }

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
  fprintf( stdout, "title( 'FLAME Tridiag_UT front-end performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc tridiag_ut_front_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

