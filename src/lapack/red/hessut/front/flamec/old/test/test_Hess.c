
#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

void time_Hess(
              int variant, int type, int n_repeats, int m,
              int nfc, int nlc,
              FLA_Obj C, FLA_Obj C_ref, FLA_Obj t,
              double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    m_input,
    m,
    p_first, p_last, p_inc,
    p,
    variant,
    n_repeats,
    i, j,
    nb_alg,
    nfc, nlc,
    n_variants = 1;
  
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
    C, C_ref, t;
  

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

  fprintf( stdout, "%c enter nfc, nlc (number of columns, initial and trailing, not processed): ", '%' );
  scanf( "%d %d", &nfc, &nlc );
  fprintf( stdout, "%c %d %d\n", '%', nfc, nlc );


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
  datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  //datatype = FLA_DOUBLE_COMPLEX;

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;

    if( m < 0 ) m = p / abs(m_input);

    for ( variant = 0; variant < n_variants; variant++ ){

      FLA_Obj_create( datatype, m, m, &C );
      FLA_Obj_create( datatype, m, m, &C_ref );
      FLA_Obj_create( datatype, m, 1, &t );

      FLA_Random_matrix( C );

      FLA_Copy_external( C, C_ref );

      fprintf( stdout, "data_hess( %d, 1:5 ) = [ %d  ", i, p );
      fflush( stdout );

      time_Hess( variant, FLA_ALG_REFERENCE, n_repeats, m, nfc, nlc,
                 C, C_ref, t, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Hess( variant, FLA_ALG_FRONT, n_repeats, m, nfc, nlc,
                 C, C_ref, t, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &C );
      FLA_Obj_free( &C_ref );
      FLA_Obj_free( &t );
    }

    fprintf( stdout, "\n" );
  }

  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  fprintf( stdout, "plot( data_hess( :,1 ), data_hess( :, 2 ), '%c:%c' ); \n",
          colors[ 0 ], ticks[ 0 ] );
  fprintf( stdout, "plot( data_hess( :,1 ), data_hess( :, 4 ), '%c:%c' ); \n",
          colors[ 1 ], ticks[ 1 ] );

  fprintf( stdout, "legend( ... \n" );

  fprintf( stdout, "'ref\\_hess', 'fla\\_hess', ... \n" );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME Hessenberg reduction front-end performance (%s)' );\n", 
           m_dim_desc );
  fprintf( stdout, "print -depsc hess_front_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );
}

