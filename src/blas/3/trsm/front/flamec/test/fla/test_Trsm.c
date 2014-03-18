
#include "FLAME.h"


#define N_PARAM_COMBOS    12

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

char* pc_str[N_PARAM_COMBOS] = { "llh", "lln", "llt", 
                                 "luh", "lun", "lut", 
                                 "rlh", "rln", "rlt", 
                                 "ruh", "run", "rut" };

void time_Trsm(
               int param_combo, int type, int n_repeats, int m, int n,
               FLA_Obj A, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    precision,
    m_input, n_input,
    m, n,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    n_param_combos = N_PARAM_COMBOS;
  
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
    A, C, C_ref;
  
  FLA_Init( );


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d%d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );


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

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;
    n = n_input;

    if( m < 0 ) m = p / abs(m_input);
    if( n < 0 ) n = p / abs(n_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){

      // Determine datatype based on trans argument.
      if ( pc_str[param_combo][2] == 'h' )
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

      // If multiplying A on the left, A is m x m; ...on the right, A is n x n.
      if ( pc_str[param_combo][0] == 'l' )
        FLA_Obj_create( datatype, m, m, 0, 0, &A );
      else
        FLA_Obj_create( datatype, n, n, 0, 0, &A );

      FLA_Obj_create( datatype, m, n, 0, 0, &C );
      FLA_Obj_create( datatype, m, n, 0, 0, &C_ref );

      if ( pc_str[param_combo][1] == 'l' )
      {
        FLA_Random_tri_matrix( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
        FLA_Random_matrix( C );
      }
      else
      {
        FLA_Random_tri_matrix( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
        FLA_Random_matrix( C );
      }

      fprintf( stdout, "data_trsm_%s( %d, 1:3 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );

      time_Trsm( param_combo, FLA_ALG_REFERENCE, n_repeats, m, n,
                 A, C, C_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );
/*
      time_Trsm( param_combo, FLA_ALG_FRONT, n_repeats, m, n,
                 A, C, C_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );
*/

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &A );
      FLA_Obj_free( &C );
      FLA_Obj_free( &C_ref );
    }

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_trsm_%s( :,1 ), data_trsm_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_trsm_%s( :,1 ), data_trsm_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_trsm\\_%s', 'fla\\_trsm\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );


  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME trsm front-end performance (%s, %s)' );\n",
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc trsm_front_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

