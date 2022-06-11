/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define N_PARAM_COMBOS    4

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

char* pc_str[N_PARAM_COMBOS] = { "nn", "nt", "tn", "tt" };

void time_Sylv(
                int param_combo, int type, int nrepeats, int m, int n,
                FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref, FLA_Obj scale,
                double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    n_threads,
    m_input, n_input,
    m, n,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i, j,
    n_param_combos = N_PARAM_COMBOS;

  dim_t nb_alg;

  int sign;
  
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
    A, B, C, C_ref, scale, isgn, norm;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter FLASH blocksize: ", '%' );
  scanf( "%lu", &nb_alg );
  fprintf( stdout, "%c %lu\n", '%', nb_alg );

  fprintf( stdout, "%c enter problem size first, last, inc: ", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c Enter sign (-1 or 1):", '%' );
  scanf( "%d", &sign );
  fprintf( stdout, "%c %d\n", '%', sign );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d %d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );

  fprintf( stdout, "%c enter the number of SuperMatrix threads: ", '%' );
  scanf( "%d", &n_threads );
  fprintf( stdout, "%c %d\n", '%', n_threads );


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

  if ( 0 < sign )
    isgn = FLA_ONE;
  else
    isgn = FLA_MINUS_ONE;

  //datatype = FLA_FLOAT;
  datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  //datatype = FLA_DOUBLE_COMPLEX;

  if ( datatype == FLA_DOUBLE || datatype == FLA_DOUBLE_COMPLEX )
  {
    FLA_Obj_create( FLA_DOUBLE, 1, 1, &scale );
    FLA_Obj_create( FLA_DOUBLE, 1, 1, &norm );
  }
  else if ( datatype == FLA_FLOAT || datatype == FLA_COMPLEX )
  {
    FLA_Obj_create( FLA_FLOAT, 1, 1, &scale );
    FLA_Obj_create( FLA_FLOAT, 1, 1, &norm );
  }

  FLASH_Queue_set_num_threads( n_threads );

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;
    n = n_input;

    if( m < 0 ) m = p / f2c_abs(m_input);
    if( n < 0 ) n = p / f2c_abs(n_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){

      FLASH_Obj_create( datatype, m, m, 1, &nb_alg, &A );
      FLASH_Obj_create( datatype, n, n, 1, &nb_alg, &B );
      FLASH_Obj_create( datatype, m, n, 1, &nb_alg, &C );
      FLASH_Obj_create( datatype, m, n, 1, &nb_alg, &C_ref );

      FLASH_Random_matrix( A );
      FLASH_Random_matrix( B );
      FLASH_Random_matrix( C );

      FLASH_Norm1( A, norm );
      FLASH_Obj_shift_diagonal( FLA_NO_CONJUGATE, norm, A );

      FLASH_Norm1( B, norm );
      if ( FLA_Obj_is( isgn, FLA_MINUS_ONE ) )
        FLA_Negate( norm );
      FLASH_Obj_shift_diagonal( FLA_NO_CONJUGATE, norm, B );


      FLASH_Copy( C, C_ref );

      fprintf( stdout, "data_sylv_%s( %d, 1:5 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );


      time_Sylv( param_combo, FLA_ALG_REFERENCE, n_repeats, m, n,
                 isgn, A, B, C, C_ref, scale, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Sylv( param_combo, FLA_ALG_FRONT, n_repeats, m, n,
                 isgn, A, B, C, C_ref, scale, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLASH_Obj_free( &A );
      FLASH_Obj_free( &B );
      FLASH_Obj_free( &C );
      FLASH_Obj_free( &C_ref );
    }

    fprintf( stdout, "\n" );
  }

  FLA_Obj_free( &scale );
  FLA_Obj_free( &norm );

  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_sylv_%s( :,1 ), data_sylv_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_sylv_%s( :,1 ), data_sylv_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_sylv\\_%s', 'fla\\_sylv\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME sylv front-end performance (%s)' );\n", m_dim_desc );
  fprintf( stdout, "print -depsc sylv_front_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );
}

