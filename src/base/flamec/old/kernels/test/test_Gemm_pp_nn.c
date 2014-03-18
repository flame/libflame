
#include "FLAME.h"


#define N_VARIANTS 1

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_UNBLOCKED 1
#define FLA_ALG_BLOCKED   2
#define FLA_ALG_RECURSIVE 3
#define FLA_ALG_OPTIMIZED 4


void time_Gemm_pp_nn(
		     int variant, int type, int nrepeats, int n, int nb_alg,
		     FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj Cref,
		     double *dtime, double *diff, double *mflops );


int main(int argc, char *argv[])
{
  int 
    m_input, k, m, n,
    nfirst, nlast, ninc,
    nb_alg,
    nrepeats,
    ntimings,
    variant,
    uplo,
    i,
    nvariants = N_VARIANTS;

  char *colors = "brkgykkk";
  char *ticks = "o+*x-+*x";

  int max_mflops=3600;

  double
    dtime,
    mflops,
    diff,
    d_n;

  FLA_Obj
    A, B, C, Cref;
  
  /* Initialize FLAME */
  FLA_Init( );

  /* Every time trial is repeated "repeat" times */
  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &nrepeats );
  fprintf( stdout, "%c %d\n", '%', nrepeats );


  /* Blocking size used for blocked algorithms */
  fprintf( stdout, "%c Enter blocking size:", '%' );
  scanf( "%d", &nb_alg );
  fprintf( stdout, "%c %d\n", '%', nb_alg );

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed */
  fprintf( stdout, "%c enter nfirst, nlast, ninc:", '%' );
  scanf( "%d%d%d", &nfirst, &nlast, &ninc );
  fprintf( stdout, "%c %d %d %d\n", '%', nfirst, nlast, ninc );

  /* Matrix dimensions */
  fprintf( stdout, "%c enter m: (-1 means m=n)", '%' );
  scanf( "%d", &m_input );
  fprintf( stdout, "%c %d\n", '%', m_input );

  fprintf( stdout, "%c enter k: ", '%' );
  scanf( "%d", &k );
  fprintf( stdout, "%c %d\n", '%', k );


  for ( n=nfirst; n<= nlast; n+=ninc ){
    
    if ( m_input == -1 ) 
      m = n;
    else
      m = m_input;

    /* Allocate space for the matrices */

    /* Note: for some operations, there is no matrix B.  We are creating
       it anyway, for uniformity sake.
       Note: If your operation does something like B := A B than think of
       it as C := A C instead:  C is always the matrix that is the output. */
    FLA_Obj_create( FLA_DOUBLE, m, k, 0, 0, &A );
    FLA_Obj_create( FLA_DOUBLE, k, n, 0, 0, &B );
    FLA_Obj_create( FLA_DOUBLE, m, n, 0, 0, &C );
    FLA_Obj_create( FLA_DOUBLE, m, n, 0, 0, &Cref );

    /* Generate random matrices A, C */
    FLA_Random_matrix( A );
    FLA_Random_matrix( B );
    FLA_Random_matrix( C );

    FLA_Copy_external( C, Cref );

    /* Time the reference implementation */
    time_Gemm_pp_nn( 0, FLA_ALG_REFERENCE, nrepeats, n, nb_alg,
		     A, B, C, Cref, &dtime, &diff, &mflops );

    fprintf( stdout, "data_REF( %d, 1:2 ) = [ %d  %6.2lf ]; \n",
	    i, n, mflops );
    fflush( stdout );


    for ( variant=1; variant<=nvariants; variant++ ){

      fprintf( stdout, "data_var%d( %d, 1:5 ) = [ %d  ", variant, i, n );
      fflush( stdout );

      time_Gemm_pp_nn( variant, FLA_ALG_UNBLOCKED, nrepeats, n, nb_alg,
		       A, B, C, Cref, &dtime, &diff, &mflops );

      fprintf( stdout, "%6.2lf %6.2le ", mflops, diff );
      fflush( stdout );

      /*       time_Gemm_pp_nn( variant, FLA_ALG_BLOCKED, nrepeats, n, nb_alg,
		    A, B, C, Cref, &dtime, &diff, &mflops );

      fprintf( stdout, "%6.2lf %6.2le ", mflops, diff );
      fflush( stdout ); */

      fprintf( stdout, " ]; \n" );
      fflush( stdout );
    }

    FLA_Obj_free( &A );
    FLA_Obj_free( &B );
    FLA_Obj_free( &C );
    FLA_Obj_free( &Cref );
    fprintf( stdout, "\n" );
  }

  /* Print the MATLAB commands to plot the data */

  /* Delete all existing figures */
  fprintf( stdout, "close all\n" );

  /* Plot the performance of the reference implementation */
  fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), '-' ); \n" );

  /* Indicate that you want to add to the existing plot */
  fprintf( stdout, "hold on\n" );

  /* Plot the data for the other variants */
  for ( i=1; i<=nvariants; i++ ){
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 2 ), '%c:%c' ); \n", 
	    i, i, colors[ i-1 ], ticks[ i-1 ] );
    //    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 4 ), '%c-.%c' ",
    //	    i, i, colors[ i-1 ], ticks[ i-1 ] );
    fprintf( stdout, " ); \n" );
  }

  fprintf( stdout, "legend( 'Reference'           ); \n" );

  for ( i=1; i<nvariants; i++ ){
    fprintf( stdout, "'unb\\_var%d', ", i, i );
    //    fprintf( stdout, "'blk\\_var%d', ", i, i );
    fprintf( stdout, "... \n", i, i );
  }

  i = nvariants;
  fprintf( stdout, "'unb\\_var%d', ", i, i );
  //    fprintf( stdout, "'blk\\_var%d', ", i, i );
  fprintf( stdout, "... \n", i, i );

  fprintf( stdout, "xlabel( 'matrix dimension n' );\n");
  fprintf( stdout, "ylabel( 'MFLOPS/sec.' );\n");
  fprintf( stdout, "axis( [ 0 %d 0 %d ] ); \n", nlast, max_mflops );
  fprintf( stdout, "print -depsc graphs.eps\n");
  fprintf( stdout, "hold off\n");

  FLA_Finalize( );
}
