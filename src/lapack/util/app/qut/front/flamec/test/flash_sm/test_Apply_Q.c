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

char* pc_str[N_PARAM_COMBOS] = { "ltc" };

void time_Trmm(
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A, FLA_Obj B, FLA_Obj B_ref, FLA_Obj t, FLA_Obj T, FLA_Obj W,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    precision,
    bm, bn,
    n_threads,
    m_input, n_input,
    m, n,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    n_param_combos = N_PARAM_COMBOS;
 
  dim_t nb_alg;
 
  char *colors = "brkgmcbrkgmcbrkgmc";
  char *ticks  = "o+*xso+*xso+*xso+*xs";
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
    A, A_save, A_flat, B, B_ref, T, T_flat, W, t;
  
  FLA_Init( );


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter FLASH blocksize: ", '%' );
  scanf( "%lu", &nb_alg );
  fprintf( stdout, "%c %lu\n", '%', nb_alg );

  fprintf( stdout, "%c enter problem size first, last, inc: ", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d%d", &m_input, &n_input );
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

  //precision = FLA_SINGLE_PRECISION;
  precision = FLA_DOUBLE_PRECISION;

  FLASH_Queue_set_num_threads( n_threads );
  //FLASH_Queue_set_verbose_output( TRUE );

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;
    n = n_input;

    if( m < 0 ) m = p / f2c_abs(m_input);
    if( n < 0 ) n = p / f2c_abs(n_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){

      // Determine datatype based on trans argument.
      if ( pc_str[param_combo][1] == 'c' )
      {
        if ( precision == FLA_SINGLE_PRECISION )
          datatype = FLA_COMPLEX;
        else
          datatype = FLA_DOUBLE_COMPLEX;
      }
      else
      {
        if ( precision == FLA_SINGLE_PRECISION )
          datatype = FLA_FLOAT;
        else
          datatype = FLA_DOUBLE;
      }

      bm = nb_alg / 4;
      bn = nb_alg;

      // If multiplying Q on the left, A is m x m; ...on the right, A is n x n.
      if ( pc_str[param_combo][0] == 'l' )
      {
        FLA_Obj_create( datatype, nb_alg, nb_alg, &A_flat );
        FLASH_Obj_create( datatype, nb_alg, nb_alg, 1, &nb_alg, &A );
        FLASH_Obj_create( datatype, nb_alg, nb_alg, 1, &nb_alg, &A_save );

        FLA_Obj_create( datatype, bm, bn, &T_flat );
        FLASH_Obj_create_ext( datatype, bm, bn, 1, &bm, &bn, &T );
        FLASH_Obj_create_ext( datatype, bm, n,  1, &bm, &bn, &W );
      }
      else
      {
        FLASH_Obj_create( datatype, n, n, 1, &nb_alg, &A );
      }

      FLASH_Obj_create( datatype, nb_alg, n, 1, &nb_alg, &B );
      FLASH_Obj_create( datatype, nb_alg, n, 1, &nb_alg, &B_ref );

      FLA_Obj_create( datatype, nb_alg, 1, &t );

      FLASH_Random_matrix( A );
      FLASH_Random_matrix( B );

      fprintf( stdout, "data_applyq_%s( %d, 1:5 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );

      FLASH_Copy( A, A_save );

      FLASH_Obj_flatten( A, A_flat );
      FLA_QR_blk_external( A_flat, t );
      FLASH_Obj_hierarchify( A_flat, A );

      time_Apply_Q( param_combo, FLA_ALG_REFERENCE, n_repeats, m, n,
                    A, B, B_ref, t, T, W, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      FLASH_Copy( A_save, A );

      FLASH_Obj_flatten( A, A_flat );
      FLA_QR_UT( A_flat, t, T_flat );
      FLASH_Obj_hierarchify( A_flat, A );
      FLASH_Obj_hierarchify( T_flat, T );

      time_Apply_Q( param_combo, FLA_ALG_FRONT, n_repeats, m, n,
                    A, B, B_ref, t, T, W, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLASH_Obj_free( &A );
      FLASH_Obj_free( &A_save );
      FLA_Obj_free( &A_flat );
      FLASH_Obj_free( &B );
      FLASH_Obj_free( &B_ref );
      FLA_Obj_free( &t );
      FLASH_Obj_free( &T );
      FLA_Obj_free( &T_flat );
      FLASH_Obj_free( &W );
    }

    fprintf( stdout, "\n" );
  }

  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_applyq_%s( :,1 ), data_applyq_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_applyq_%s( :,1 ), data_applyq_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_applyq\\_%s', 'fla\\_applyq\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );


  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME applyq front-end performance (%s, %s)' );\n",
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc applyq_front_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );

  return 0;
}

