/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


#define N_PARAM_COMBOS    9

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

char* pc_str[N_PARAM_COMBOS] = { "cc", "cn", "ct", 
                                 "nc", "nn", "nt", 
                                 "tc", "tn", "tt" };

integer sgemm_p;
integer cgemm_p;
integer dgemm_p;
integer zgemm_p;
integer sgemm_q;
integer cgemm_q;
integer dgemm_q;
integer zgemm_q;
integer sgemm_r;
integer cgemm_r;
integer dgemm_r;
integer zgemm_r;

void time_Gemm(
               integer param_combo, integer type, integer n_repeats, integer m, integer k, integer n,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


int main(integer argc, char *argv[])
{
  integer 
    datatype,
    precision,
    m_input, k_input, n_input,
    m, k, n,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    n_param_combos = N_PARAM_COMBOS;
  
  char *colors = "brkgmcbrkgmcbrkgmc";
  char *ticks  = "o+*xso+*xso+*xso+*xs";
  char m_dim_desc[14];
  char k_dim_desc[14];
  char n_dim_desc[14];
  char m_dim_tag[10];
  char k_dim_tag[10];
  char n_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, Ad, Az, B, Bd, Bz, C, Cd, Cz, C_ref, indexd, indexz;
  FLA_Obj alpha0d, alpha0z, alpha1d, alpha1z, normd, normz;
  FLA_Obj alphad, alphaz, betad, betaz, rhod, rhoz;
  FLA_Obj xd, xz, yd, yz;
  
  FLA_Init( );



  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m k n (-1 means bind to problem size): ", '%' );
  scanf( "%d%d%d", &m_input, &k_input, &n_input );
  fprintf( stdout, "%c %d %d %d\n", '%', m_input, k_input, n_input );


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
  if     ( k_input >  0 ) {
    sprintf( k_dim_desc, "k = %d", k_input );
    sprintf( k_dim_tag,  "k%dc", k_input);
  }
  else if( k_input <  -1 ) {
    sprintf( k_dim_desc, "k = p/%d", -k_input );
    sprintf( k_dim_tag,  "k%dp", -k_input );
  }
  else if( k_input == -1 ) {
    sprintf( k_dim_desc, "k = p" );
    sprintf( k_dim_tag,  "k%dp", 1 );
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
    k = k_input;
    n = n_input;

    if( m < 0 ) m = p / f2c_abs(m_input);
    if( k < 0 ) k = p / f2c_abs(k_input);
    if( n < 0 ) n = p / f2c_abs(n_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){

      // Determine datatype based on trans argument.
      if ( pc_str[param_combo][0] == 'c' ||
           pc_str[param_combo][1] == 'c' )
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

      // If transposing A, switch dimensions.
      if ( pc_str[param_combo][0] == 'n' )
        FLA_Obj_create( datatype, m, k, 0, 0, &A );
      else
        FLA_Obj_create( datatype, k, m, 0, 0, &A );
      
      // If transposing B, switch dimensions.
      if ( pc_str[param_combo][1] == 'n' )
        FLA_Obj_create( datatype, k, n, 0, 0, &B );
      else
        FLA_Obj_create( datatype, n, k, 0, 0, &B );

      FLA_Obj_create( datatype, m, n, 0, 0, &C );
      FLA_Obj_create( datatype, m, n, 0, 0, &C_ref );

      FLA_Random_matrix( A );
      FLA_Random_matrix( B );
      FLA_Random_matrix( C );

      FLA_Copy_external( C, C_ref );

      
      fprintf( stdout, "data_gemm_%s( %d, 1:5 ) = [ %4d %4d %4d  ", pc_str[param_combo], i, m, k, n );
      fflush( stdout );

      time_Gemm( param_combo, FLA_ALG_REFERENCE, n_repeats, m, k, n,
                 A, B, C, C_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );
/*
      time_Gemm( param_combo, FLA_ALG_FRONT, n_repeats, m, k, n,
                 A, B, C, C_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );
*/

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &A );
      FLA_Obj_free( &B );
      FLA_Obj_free( &C );
      FLA_Obj_free( &C_ref );
    }

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_gemm_%s( :,1 ), data_gemm_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_gemm_%s( :,1 ), data_gemm_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_gemm\\_%s', 'fla\\_gemm\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );


  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME gemm front-end performance (%s, %s, %s)' );\n",
           m_dim_desc, k_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc gemm_front_%s_%s_%s.eps\n", m_dim_tag, k_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

