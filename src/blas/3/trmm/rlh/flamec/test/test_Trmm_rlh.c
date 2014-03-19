/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"



#define N_VARIANTS 4

#define FLA_ALG_REFERENCE  0
#define FLA_ALG_UNBLOCKED  1
#define FLA_ALG_BLOCKED    2
#define FLA_ALG_OPTIMIZED  3


void time_Trmm_rlh(
               int variant, int type, int n_repeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    m_input, n_input,
    m, n,
    p_first, p_last, p_inc,
    p,
    nb_alg,
    n_repeats,
    variant,
    i, j,
    datatype,
    n_variants = N_VARIANTS;
  
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
    A, B, C, C_ref;
  
  /* Initialize FLAME */
  FLA_Init( );


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c Enter blocking size: ", '%' );
  scanf( "%d", &nb_alg );
  fprintf( stdout, "%c %d\n", '%', nb_alg );

  fprintf( stdout, "%c enter problem size first, last, inc: ", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d%d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );


  /* Delete all existing data structures */
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



  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;
    n = n_input;

    if( m < 0 ) m = p / abs(m_input);
    if( n < 0 ) n = p / abs(n_input);

    //datatype = FLA_COMPLEX;
    datatype = FLA_DOUBLE_COMPLEX;

    /* Allocate space for the matrices */
    FLA_Obj_create( datatype, n, n, &A );
    FLA_Obj_create( datatype, m, n, &C );
    FLA_Obj_create( datatype, m, n, &C_ref );

    /* Generate random matrices A, C */
    FLA_Random_tri_matrix( FLA_UPPER_TRIANGULAR, FLA_UNIT_DIAG, A );
    FLA_Random_matrix( C );

    FLA_Copy_external( C, C_ref );


    /* Time the reference implementation */
    time_Trmm_rlh( 0, FLA_ALG_REFERENCE, n_repeats, p, nb_alg,
                  A, B, C, C_ref, &dtime, &diff, &gflops );

    fprintf( stdout, "data_REF( %d, 1:2 ) = [ %d  %6.3lf ]; \n", i, p, gflops );
    fflush( stdout );

    for ( variant = 1; variant <= n_variants; variant++ ){
      
      //fprintf( stdout, "data_var%d( %d, 1:7 ) = [ %d  ", variant, i, p );
      fprintf( stdout, "data_var%d( %d, 1:5 ) = [ %d  ", variant, i, p );
      fflush( stdout );


      time_Trmm_rlh( variant, FLA_ALG_UNBLOCKED, n_repeats, p, nb_alg,
                    A, B, C, C_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Trmm_rlh( variant, FLA_ALG_BLOCKED, n_repeats, p, nb_alg,
                    A, B, C, C_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      //time_Trmm_rlh( variant, FLA_ALG_OPTIMIZED, n_repeats, p, nb_alg,
      //              A, B, C, C_ref, &dtime, &diff, &gflops );

      //fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      //fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );
    }
    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &C );
    FLA_Obj_free( &C_ref );
  }

  /* Print the MATLAB commands to plot the data */

  /* Delete all existing figures */
  fprintf( stdout, "figure;\n" );

  /* Plot the performance of the reference implementation */
  fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), '-' ); \n" );

  /* Indicate that you want to add to the existing plot */
  fprintf( stdout, "hold on;\n" );

  /* Plot the data for the other numbers of threads */
  for ( i = 1; i <= n_variants; i++ ) {
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 2 ), '%c:%c' ); \n",
            i, i, colors[ i-1 ], ticks[ i-1 ] );
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 4 ), '%c-.%c' ); \n",
            i, i, colors[ i-1 ], ticks[ i-1 ] );
    //fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 6 ), '%c--%c' ); \n",
    //        i, i, colors[ i-1 ], ticks[ i-1 ] );
  }

  fprintf( stdout, "legend( ... \n" );
  fprintf( stdout, "'Reference', ... \n" );

  for ( i = 1; i < n_variants; i++ )
    //fprintf( stdout, "'unb\\_var%d', 'blk\\_var%d', 'opt\\_var%d', ... \n", i, i, i );
    fprintf( stdout, "'unb\\_var%d', 'blk\\_var%d', ... \n", i, i );

  i = n_variants;
  fprintf( stdout, "'unb\\_var%d', 'blk\\_var%d' ); \n", i, i );


  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME trmm\\_rlc performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc trmm_rlc_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );
}

