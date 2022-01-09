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

char* pc_str[N_PARAM_COMBOS] = { "lhfc", "lnfr" };

void time_Chol(
                integer param_combo, integer type, integer nrepeats, integer m, integer n,
                FLA_Obj A, FLA_Obj A_ref, FLA_Obj T, FLA_Obj t_ref, FLA_Obj B, FLA_Obj B_ref, FLA_Obj X, FLA_Obj X_ref, FLA_Obj W,
                double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  integer 
    m_input, n_input,
    m, n,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    n_param_combos = N_PARAM_COMBOS;

  FLA_Datatype datatype;
  
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
    A, A_ref, t_ref, t, T, T_ref, B, B_ref, X, X_ref, W;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter problem size first, last, inc: ", '%' );
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


  //datatype = FLA_FLOAT;
  //datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  datatype = FLA_DOUBLE_COMPLEX;

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;
    n = n_input;

    if( m < 0 ) m = p / f2c_abs(m_input);
    if( n < 0 ) n = p / f2c_abs(n_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){

      if ( pc_str[param_combo][0] == 'l' )
      {
        FLA_Obj_create( datatype, m, m, 0, 0, &A );
        FLA_Obj_create( datatype, m, m, 0, 0, &A_ref );

        FLA_Obj_create( datatype, m, 1, 0, 0, &t );
        FLA_Obj_create( datatype, m, 1, 0, 0, &t_ref );

        FLA_Obj_create( datatype, m, n, 0, 0, &B );
        FLA_Obj_create( datatype, m, n, 0, 0, &B_ref );

        FLA_Obj_create( datatype, m, n, 0, 0, &X );
        FLA_Obj_create( datatype, m, n, 0, 0, &X_ref );
      }
      else
      {
        FLA_Obj_create( datatype, n, n, 0, 0, &A );
        FLA_Obj_create( datatype, n, n, 0, 0, &A_ref );

        FLA_Obj_create( datatype, n, 1, 0, 0, &t );
        FLA_Obj_create( datatype, n, 1, 0, 0, &t_ref );

        FLA_Obj_create( datatype, m, n, 0, 0, &B );
        FLA_Obj_create( datatype, m, n, 0, 0, &B_ref );

        FLA_Obj_create( datatype, m, n, 0, 0, &X );
        FLA_Obj_create( datatype, m, n, 0, 0, &X_ref );
      }


      FLA_Obj_create( datatype, 32, m, 0, 0, &T );
      FLA_Obj_create( datatype, 32, n, 0, 0, &W );


/*
      FLA_Obj_create( datatype, 4, m, &T );
      FLA_Obj_create( datatype, 4, m, &T_ref );
      FLA_Obj_create( datatype, 4, n, &W );
*/

/*
      FLA_Obj_create( datatype, 2, m, &T );
      FLA_Obj_create( datatype, 2, m, &T_ref );
      FLA_Obj_create( datatype, 2, n, &W );
*/

      FLA_Random_matrix( A );
      FLA_Copy_external( A, A_ref );

      //FLA_Obj_show( "A_orig:", A, "%12.4e", "" );

      FLA_Random_matrix( B );
      FLA_Copy_external( B, B_ref );

      if ( pc_str[param_combo][2] == 'c' )
      {
        //FLA_QR_blk_external( A_ref, t_ref );
        FLA_QR_UT( A, T );
      }
      else if ( pc_str[param_combo][2] == 'r' )
      {
        //FLA_LQ_blk_external( A_ref, t_ref );
        FLA_LQ_UT( A, T );
      }

      fprintf( stdout, "data_applyq_%s( %d, 1:5 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );


      time_Chol( param_combo, FLA_ALG_REFERENCE, n_repeats, m, n,
                 A, A_ref, T, t_ref, B, B_ref, X, X_ref, W, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Chol( param_combo, FLA_ALG_FRONT, n_repeats, m, n,
                 A, A_ref, T, t_ref, B, B_ref, X, X_ref, W, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &A );
      FLA_Obj_free( &A_ref );
      FLA_Obj_free( &t );
      FLA_Obj_free( &t_ref );
      FLA_Obj_free( &B );
      FLA_Obj_free( &B_ref );
      FLA_Obj_free( &X );
      FLA_Obj_free( &X_ref );
      FLA_Obj_free( &T );
      FLA_Obj_free( &W );
    }

    fprintf( stdout, "\n" );
  }

/*
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
  fprintf( stdout, "title( 'FLAME apply_q front-end performance (%s, %s)' );\n", m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc apply_q_front_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize();

  return 0;
}

