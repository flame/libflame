/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


#define N_PARAM_COMBOS    6

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

char* pc_str[N_PARAM_COMBOS] = { "lc", "ln", "lt", 
                                 "uc", "un", "ut" }; 

void time_Trsv(
               integer param_combo, integer type, integer n_repeats, integer m,
               FLA_Obj A, FLA_Obj x, FLA_Obj x_ref,
               double *dtime, double *diff, double *gflops );


int main(integer argc, char *argv[])
{
  integer 
    datatype,
    precision,
    nb_alg,
    n_threads,
    m_input,
    m,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    n_param_combos = N_PARAM_COMBOS;
  integer one = 1;
  
  char *colors = "brkgmcbrkgmcbrkgmc";
  char *ticks  = "o+*xso+*xso+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, x, x_ref;
  
  FLA_Init( );


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter FLASH blocksize: ", '%' );
  scanf( "%d", &nb_alg );
  fprintf( stdout, "%c %d\n", '%', nb_alg );

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

  //precision = FLA_SINGLE_PRECISION;
  precision = FLA_DOUBLE_PRECISION;

  FLASH_Queue_set_num_threads( n_threads );

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;

    if( m < 0 ) m = p / f2c_abs(m_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){

      // Determine datatype based on trans argument.
      if ( pc_str[param_combo][2] == 'c' )
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

      FLASH_Obj_create( datatype, m, m, 1, &nb_alg, &A );

      FLASH_Obj_create_ext( datatype, m, 1, 1, &nb_alg, &one, &x );
      FLASH_Obj_create_ext( datatype, m, 1, 1, &nb_alg, &one, &x_ref );

      if ( pc_str[param_combo][1] == 'l' )
        FLASH_Random_spd_matrix( FLA_LOWER_TRIANGULAR, A );
      else
        FLASH_Random_spd_matrix( FLA_UPPER_TRIANGULAR, A );

      FLASH_Random_matrix( x );

      FLASH_Copy( x, x_ref );
      
      fprintf( stdout, "data_trsv_%s( %d, 1:5 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );

      time_Trsv( param_combo, FLA_ALG_REFERENCE, n_repeats, m,
                 A, x, x_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Trsv( param_combo, FLA_ALG_FRONT, n_repeats, m,
                 A, x, x_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLASH_Obj_free( &A );
      FLASH_Obj_free( &x );
      FLASH_Obj_free( &x_ref );
    }

    fprintf( stdout, "\n" );
  }

  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_trsv_%s( :,1 ), data_trsv_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_trsv_%s( :,1 ), data_trsv_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_trsv\\_%s', 'fla\\_trsv\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );


  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME trsv front-end performance (%s)' );\n",
           m_dim_desc );
  fprintf( stdout, "print -depsc trsv_front_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );

  return 0;
}

