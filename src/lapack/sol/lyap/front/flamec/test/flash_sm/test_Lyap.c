
#include "FLAME.h"

#define N_PARAM_COMBOS    2

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

char* pc_str[N_PARAM_COMBOS] = { "n", "h" };

void time_Lyap(
                int param_combo, int type, int nrepeats, int m,
                FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale,
                double *dtime, double *diff, double *gflops );


int main( int argc, char *argv[] )
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
    i, j,
    n_param_combos = N_PARAM_COMBOS;

  int sign;

  dim_t b_flash;
  
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
    A, A_flat, C, C_flat, scale, isgn, norm;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c Enter sign (-1 or 1):", '%' );
  scanf( "%d", &sign );
  fprintf( stdout, "%c %d\n", '%', sign );

  fprintf( stdout, "%c enter FLASH blocksize: ", '%' );
  scanf( "%u", &b_flash );
  fprintf( stdout, "%c %u\n", '%', b_flash );

  fprintf( stdout, "%c enter problem size first, last, inc: ", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m (-1 means bind to problem size): ", '%' );
  scanf( "%d", &m_input );
  fprintf( stdout, "%c %d\n", '%', m_input );

  fprintf( stdout, "%c enter the number of SuperMatrix threads: ", '%' );
  scanf( "%d", &n_threads );
  fprintf( stdout, "%c %d\n", '%', n_threads );


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

  if ( 0 < sign )
    isgn = FLA_ONE;
  else
    isgn = FLA_MINUS_ONE;

  //datatype = FLA_FLOAT;
  datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  //datatype = FLA_DOUBLE_COMPLEX;

  FLASH_Queue_set_num_threads( n_threads );
  //FLASH_Queue_disable();

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;

    if( m < 0 ) m = p / abs(m_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){

      FLA_Obj_create( datatype, m, m, 0, 0, &A_flat );
      FLA_Obj_create( datatype, m, m, 0, 0, &C_flat );

      FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A_flat ), 1, 1, 0, 0, &scale );
      FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A_flat ), 1, 1, 0, 0, &norm );

      FLA_Random_tri_matrix( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A_flat );
      FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A_flat );
      FLA_Norm1( A_flat, norm );
      FLA_Shift_diag( FLA_NO_CONJUGATE, norm, A_flat );

      FLA_Random_matrix( C_flat );
      FLA_Hermitianize( FLA_UPPER_TRIANGULAR, C_flat );

      FLASH_Obj_create_hier_copy_of_flat( A_flat, 1, &b_flash, &A );
      FLASH_Obj_create_hier_copy_of_flat( C_flat, 1, &b_flash, &C );


      fprintf( stdout, "data_lyap_%s( %d, 1:3 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );

/*
      time_Lyap( param_combo, FLA_ALG_REFERENCE, n_repeats, m, n,
                 isgn, A, B, C, C_ref, scale, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );
*/
      time_Lyap( param_combo, FLA_ALG_FRONT, n_repeats, m,
                 isgn, A, C, scale, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLASH_Obj_free( &A );
      FLASH_Obj_free( &C );
      FLA_Obj_free( &A_flat );
      FLA_Obj_free( &C_flat );
      FLA_Obj_free( &scale );
      FLA_Obj_free( &norm );
    }

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_lyap_%s( :,1 ), data_lyap_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_lyap_%s( :,1 ), data_lyap_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_lyap\\_%s', 'fla\\_lyap\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME lyap front-end performance (%s)' );\n", m_dim_desc );
  fprintf( stdout, "print -depsc lyap_front_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/
  FLA_Finalize( );

  return 0;
}

