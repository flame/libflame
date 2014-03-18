
#include "FLAME.h"

#define N_PARAM_COMBOS    4

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

char* pc_str[N_PARAM_COMBOS] = { "ln", "lu",
                                 "un", "uu" };

void time_Trinv(
                 int param_combo, int type, int n_repeats, int m, FLA_Uplo uplo, FLA_Diag diag,
                 FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj norm,
                 double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    m_input,
    m,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    n_param_combos = N_PARAM_COMBOS;

  FLA_Uplo uplo;
  FLA_Diag diag;
  
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
    A, b, b_orig, norm;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

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

    if( m < 0 ) m = p / abs(m_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){
      
      //FLA_Obj_create( datatype, m, m, 0, 0, &A );
      FLA_Obj_create( datatype, m, m, m, 1, &A );
      FLA_Obj_create( datatype, m, 1, 0, 0, &b );
      FLA_Obj_create( datatype, m, 1, 0, 0, &b_orig );

      if ( FLA_Obj_is_single_precision( A ) )
        FLA_Obj_create( FLA_FLOAT, 1, 1, 0, 0, &norm );
      else
        FLA_Obj_create( FLA_DOUBLE, 1, 1, 0, 0, &norm );

      FLA_Param_map_netlib_to_flame_uplo( &pc_str[param_combo][0], &uplo );
      FLA_Param_map_netlib_to_flame_diag( &pc_str[param_combo][1], &diag );

      FLA_Random_tri_matrix( uplo, diag, A );
      FLA_Random_matrix( b );
      FLA_Copy_external( b, b_orig );

      fprintf( stdout, "data_trinv_%s( %d, 1:5 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );

/*
      time_Trinv( param_combo, FLA_ALG_REFERENCE, n_repeats, m, uplo, diag,
                  A, b, b_orig, norm, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );
*/

      time_Trinv( param_combo, FLA_ALG_FRONT, n_repeats, m, uplo, diag,
                  A, b, b_orig, norm, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &A );
      FLA_Obj_free( &b );
      FLA_Obj_free( &b_orig );
      FLA_Obj_free( &norm );
    }

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_trinv_%s( :,1 ), data_trinv_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_trinv_%s( :,1 ), data_trinv_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_trinv\\_%s', 'fla\\_trinv\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthWest' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME trinv front-end performance (%s)' );\n", m_dim_desc );
  fprintf( stdout, "print -depsc trinv_front_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

