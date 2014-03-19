/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE     0
#define FLA_ALG_UNBLOCKED     1
#define FLA_ALG_UNB_OPT       2
#define FLA_ALG_UNB_ASM       3
#define FLA_ALG_BLOCKED       4

void fill_cs( FLA_Obj G );
void fill_a( int ij, FLA_Obj A );

void time_Apply_G_rf(
               int variant, int type, int n_repeats, int m, int k, int n, int b_alg,
               FLA_Obj A, FLA_Obj A_ref, FLA_Obj G, FLA_Obj P,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    m_input, k_input, n_input,
    m, k, n,
    p_first, p_last, p_inc,
    p,
    b_alg,
    variant,
    n_repeats,
    i,
    datatype, dt_real, dt_comp,
    n_variants = 9;
  
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];
  char k_dim_desc[14];
  char k_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, A_ref, G, P;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c Enter blocking size:", '%' );
  scanf( "%d", &b_alg );
  fprintf( stdout, "%c %d\n", '%', b_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n k (-1 means bind to problem size): ", '%' );
  scanf( "%d %d %d", &m_input, &n_input, &k_input );
  fprintf( stdout, "%c %d %d %d\n", '%', m_input, n_input, k_input );


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


  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;
    k = k_input;
    n = n_input;

    if( m < 0 ) m = p / abs(m_input);
    if( k < 0 ) k = p / abs(k_input);
    if( n < 0 ) n = p / abs(n_input);

    //datatype = FLA_FLOAT;
    //datatype = FLA_DOUBLE;
    //datatype = FLA_COMPLEX;
    datatype = FLA_DOUBLE_COMPLEX;


    FLA_Obj_create( datatype, m,   n, 0, 0, &A );
    FLA_Obj_create( datatype, m,   n, 0, 0, &A_ref );

	if ( FLA_Obj_is_double_precision( A ) ) dt_comp = FLA_DOUBLE_COMPLEX;
	else                                    dt_comp = FLA_COMPLEX;

    FLA_Obj_create( dt_comp, n-1, k, 0, 0, &G );
    FLA_Obj_create( dt_comp, n-1, k, 0, 0, &P );

    FLA_Random_matrix( A );
    //FLA_Set_to_identity( A );
    //FLA_Set( FLA_ZERO, A );
    //FLA_Set_diag( FLA_TWO, A );
    //FLA_Random_tri_matrix( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
    //FLA_Random_tri_matrix( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
    FLA_Random_matrix( G );
    //fill_cs( G );

/*
{
  FLA_Obj GTL, GTR;
  FLA_Obj GBL, GBR;

  FLA_Part_2x2( G,   &GTL, &GTR,
                     &GBL, &GBR,     6, 1, FLA_TL );
  FLA_Obj_show( "GTL", GTL, "%9.2e %9.2e ", "" );
}
*/

/*
    time_Apply_G_rf( 0, FLA_ALG_REFERENCE, n_repeats, m, nb_alg,
               A, t, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REF( %d, 1:2 ) = [ %d  %6.3lf %6.2le ]; \n", i, p, gflops, diff );
    fflush( stdout );
*/

    for ( variant = 1; variant <= n_variants; variant++ ){
      
      fprintf( stdout, "data_var%d( %d, 1:7 ) = [ %d ", variant, i, p );
      fflush( stdout );

      time_Apply_G_rf( variant, FLA_ALG_UNB_OPT, n_repeats, m, k, n, b_alg,
                       A, A_ref, G, P, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Apply_G_rf( variant, FLA_ALG_UNB_ASM, n_repeats, m, k, n, b_alg,
                       A, A_ref, G, P, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Apply_G_rf( variant, FLA_ALG_BLOCKED, n_repeats, m, k, n, b_alg,
                       A, A_ref, G, P, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, "];\n" );
      fflush( stdout );
    }

    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &A_ref );
    FLA_Obj_free( &G );
    FLA_Obj_free( &P );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), '-' ); \n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 1; i <= n_variants; i++ ) {
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 2 ), '%c:%c' ); \n",
            i, i, colors[ i-1 ], ticks[ i-1 ] );
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 4 ), '%c-.%c' ); \n",
            i, i, colors[ i-1 ], ticks[ i-1 ] );
  }

  fprintf( stdout, "legend( ... \n" );
  fprintf( stdout, "'Reference', ... \n" );

  for ( i = 1; i < n_variants; i++ )
    fprintf( stdout, "'unb\\_var%d', 'blk\\_var%d', ... \n", i, i );
  fprintf( stdout, "'unb\\_var%d', 'blk\\_var%d' ); \n", i, i );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME Apply_G performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc tridiag_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

void fill_cs( FLA_Obj G )
{
  FLA_Obj GL,    GR,       G0,  g1,  G2;

  FLA_Obj g1T,
          g1B;

  FLA_Part_1x2( G,    &GL,  &GR,      0, FLA_LEFT );

  while ( FLA_Obj_width( GL ) < FLA_Obj_width( G ) ){

    FLA_Repart_1x2_to_1x3( GL,  /**/ GR,        &G0, /**/ &g1, &G2,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( g1,  &g1T,
                       &g1B,    FLA_Obj_width( G0 ), FLA_TOP );
    FLA_Set( FLA_ONE,  g1T );

    FLA_Part_2x1( g1,  &g1T,
                       &g1B,    FLA_Obj_width( G0 ), FLA_BOTTOM );
//printf( "n(G0) = %d\n", FLA_Obj_width( G0 ) );
//printf( "m(g1B) = %d\n", FLA_Obj_length( g1B ) );
    FLA_Set( FLA_ONE,  g1B );
//if ( FLA_Obj_length( g1B ) == 8 ) FLA_Obj_show( "g1", g1, "%9.2e + %9.2e ", "" );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &GL,  /**/ &GR,        G0, g1, /**/ G2,
                              FLA_LEFT );
  }

}

