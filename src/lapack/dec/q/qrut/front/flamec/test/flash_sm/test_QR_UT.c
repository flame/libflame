
#include "FLAME.h"


#define N_PARAM_COMBOS    1

#define FLA_ALG_FRONT     0

char* pc_str[N_PARAM_COMBOS] = { "" };

void time_QR_UT(
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A, FLA_Obj TW, FLA_Obj b, FLA_Obj x,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    n_threads,
    m_input,
    m,
    n_input,
    n,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    n_param_combos = N_PARAM_COMBOS;

  dim_t
    b_flash,
    b_alg;

  char *colors = "brkgmcbrkgmcbrkgmc";
  char *ticks  = "o+*xso+*xso+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj A, TW, b, x;
  FLA_Obj A_flat, b_flat, x_flat;
  
  FLA_Init( );


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter FLASH blocksize: ", '%' );
  scanf( "%u", &b_flash );
  fprintf( stdout, "%c %u\n", '%', b_flash );

  fprintf( stdout, "%c enter problem size first, last, inc: ", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d %d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );

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

  //datatype = FLA_FLOAT;
  datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  //datatype = FLA_DOUBLE_COMPLEX;

  FLASH_Queue_set_num_threads( n_threads );
  //FLASH_Queue_set_verbose_output( TRUE );
  //FLA_Check_error_level_set( FLA_NO_ERROR_CHECKING );
  //FLASH_Queue_disable();

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;
    n = n_input;

    if ( m < 0 ) m = p * abs(m_input);
    if ( n < 0 ) n = p * abs(n_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ )
    {
      FLA_Obj_create( datatype, m, n, 0, 0, &A_flat );
      FLA_Obj_create( datatype, n, 1, 0, 0, &x_flat );
      FLA_Obj_create( datatype, m, 1, 0, 0, &b_flat );

      FLA_Random_matrix( A_flat );
      FLA_Random_matrix( b_flat );

      FLASH_QR_UT_create_hier_matrices( A_flat, 1, &b_flash, &A, &TW );
      FLASH_Obj_create_hier_copy_of_flat( b_flat, 1, &b_flash, &b );
      FLASH_Obj_create_hier_copy_of_flat( x_flat, 1, &b_flash, &x );


      fprintf( stdout, "data_qrut_%s( %d, 1:5 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );

      time_QR_UT( param_combo, FLA_ALG_FRONT, n_repeats, m, n,
                  A, TW, b, x, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &A_flat );
      FLA_Obj_free( &b_flat );
      FLA_Obj_free( &x_flat );

      FLASH_Obj_free( &A );
      FLASH_Obj_free( &TW );
      FLASH_Obj_free( &b );
      FLASH_Obj_free( &x );
    }

  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_qrut_%s( :,1 ), data_qrut_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_qrut_%s( :,1 ), data_qrut_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_qrut\\_%s', 'fla\\_qrut\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );


  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME qrut front-end performance (%s)' );\n",
           m_dim_desc );
  fprintf( stdout, "print -depsc qrut_front_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

