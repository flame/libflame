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

char* pc_str[N_PARAM_COMBOS] = { "il", "iu", "nl", "nu" };

void time_Eig_gest(
                int param_combo, int type, int n_repeats, int n,
                FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj norm,
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
    n_param_combos = N_PARAM_COMBOS;

  dim_t b_flash;
  dim_t n_threads;

  FLA_Datatype datatype;
  FLA_Uplo     uplo;
  FLA_Inv      inv;
  
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
    A, B, norm;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter FLASH blocksize: ", '%' );
  scanf( "%lu", &b_flash );
  fprintf( stdout, "%c %lu\n", '%', b_flash );

  fprintf( stdout, "%c enter problem size first, last, inc: ", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m (-1 means bind to problem size): ", '%' );
  scanf( "%d", &m_input );
  fprintf( stdout, "%c %d\n", '%', m_input );

  fprintf( stdout, "%c enter the number of SuperMatrix threads: ", '%' );
  scanf( "%lu", &n_threads );
  fprintf( stdout, "%c %lu\n", '%', n_threads );

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


  //datatype = FLA_FLOAT;
  datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  //datatype = FLA_DOUBLE_COMPLEX;

  FLASH_Queue_set_num_threads( n_threads );

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;

    if( m < 0 ) m = p / f2c_abs(m_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){

      if ( pc_str[param_combo][0] == 'i' ) inv = FLA_INVERSE;
      else                                 inv = FLA_NO_INVERSE;

      if ( pc_str[param_combo][1] == 'l' ) uplo = FLA_LOWER_TRIANGULAR;
      else                                 uplo = FLA_UPPER_TRIANGULAR;

      FLASH_Obj_create( datatype, m, m, 1, &b_flash, &A );
      FLASH_Obj_create( datatype, m, m, 1, &b_flash, &B );

      FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

      FLASH_Random_spd_matrix( uplo, A );
      FLASH_Hermitianize( uplo, A );
      
      FLASH_Random_spd_matrix( uplo, B );
      FLASH_Chol( uplo, B );

      
      fprintf( stdout, "data_eig_gest_%s( %d, 1:3 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );

      time_Eig_gest( param_combo, FLA_ALG_FRONT, n_repeats, m,
                     inv, uplo, A, B, norm, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLASH_Obj_free( &A );
      FLASH_Obj_free( &B );
      FLA_Obj_free( &norm );
    }

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_eig_gest_%s( :,1 ), data_eig_gest_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_eig_gest_%s( :,1 ), data_eig_gest_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_eig_gest\\_%s', 'fla\\_eig_gest\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME eig_gest front-end performance (%s)' );\n", m_dim_desc );
  fprintf( stdout, "print -depsc eig_gest_front_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize();

  return 0;
}

