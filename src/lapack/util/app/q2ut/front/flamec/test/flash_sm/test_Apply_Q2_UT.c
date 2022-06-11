/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


#define N_PARAM_COMBOS    1

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

char* pc_str[N_PARAM_COMBOS] = { "" };

void time_Apply_Q2_UT(
               int param_combo, int type, int nrepeats, int m,
               FLA_Obj A, FLA_Obj A_ref, FLA_Obj B, FLA_Obj C, FLA_Obj D, FLA_Obj E,
               FLA_Obj T, FLA_Obj TB, FLA_Obj TD, FLA_Obj TE, FLA_Obj W,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    n_threads,
    m_input,
    m,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    n_param_combos = N_PARAM_COMBOS;

  dim_t
    nb_flash,
    nb_alg;
  
  char *colors = "brkgmcbrkgmcbrkgmc";
  char *ticks  = "o+*xso+*xso+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj A_flat;
  FLA_Obj A, B, C, D, E;
  FLA_Obj T, TB, TD, TE, W;
  FLA_Obj A_ref;
  
  FLA_Init( );


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter algorithmic blocksize: ", '%' );
  scanf( "%lu", &nb_alg );
  fprintf( stdout, "%c %lu\n", '%', nb_alg );

  fprintf( stdout, "%c enter FLASH blocksize: ", '%' );
  scanf( "%lu", &nb_flash );
  fprintf( stdout, "%c %lu\n", '%', nb_flash );

  fprintf( stdout, "%c enter problem size first, last, inc: ", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m (-1 means bind to problem size): ", '%' );
  scanf( "%d", &m_input );
  fprintf( stdout, "%c %d\n", '%', m_input );

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

  //datatype = FLA_FLOAT;
  //datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  datatype = FLA_DOUBLE_COMPLEX;

  FLASH_Queue_set_num_threads( n_threads );

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;

    if( m < 0 ) m = p / f2c_abs(m_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ )
    {

      FLA_Obj_create( datatype, m, m, &A_flat );
      FLA_Random_matrix( A_flat );

      FLASH_QR_UT_inc_create_hier_matrices( A_flat, 1, &nb_flash, nb_alg, &A, &T );
      FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_ref );

      FLA_Part_2x2( A,    &B, &C,
                          &D, &E,     1, 1, FLA_TL );

      FLA_Part_2x2( T,    &TB, &W,
                          &TD, &TE,     1, 1, FLA_TL );

      FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, *(FLASH_OBJ_PTR_AT(B)) );



      fprintf( stdout, "data_qr2ut_%s( %d, 1:5 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );

      FLASH_Queue_set_num_threads( 1 );
      time_Apply_Q2_UT( param_combo, FLA_ALG_REFERENCE, n_repeats, m,
                          A, A_ref, B, C, D, E, T, TB, TD, TE, W, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      FLASH_Queue_set_num_threads( n_threads );
      time_Apply_Q2_UT( param_combo, FLA_ALG_FRONT, n_repeats, m,
                          A, A_ref, B, C, D, E, T, TB, TD, TE, W, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &A_flat );

      FLASH_Obj_free( &A );
      FLASH_Obj_free( &T );
    }

    fprintf( stdout, "\n" );
  }

  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_qr2ut_%s( :,1 ), data_qr2ut_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_qr2ut_%s( :,1 ), data_qr2ut_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_qr2ut\\_%s', 'fla\\_qr2ut\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );


  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME qr2ut front-end performance (%s)' );\n",
           m_dim_desc );
  fprintf( stdout, "print -depsc qr2ut_front_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );

  return 0;
}

