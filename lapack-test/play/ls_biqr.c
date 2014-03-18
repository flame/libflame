#include "FLAME.h"
#define EPS 1.e-10

#define FORT_FUNC(prefix, name) prefix ## name

typedef float testtype;
#define TESTTYPE FLA_FLOAT

//typedef double testtype;
//#define TESTTYPE FLA_DOUBLE

//typedef dcomplex testtype;
//#define TESTTYPE FLA_DOUBLE_COMPLEX


int main( int argc, char** argv ) {
  FLA_Datatype datatype = TESTTYPE;
  FLA_Obj      A, A_flame, A_lapack;
  int          m;
  FLA_Error    init_result; 

  FLA_Obj TU, TV, U_flame, V_flame, d_flame, e_flame, B_flame;
  FLA_Obj tauq, taup, d_lapack, e_lapack, U_lapack, V_lapack, W, B_lapack;
  testtype *buff_tauq, *buff_taup, *buff_d_lapack, *buff_e_lapack, 
    *buff_W, *buff_A_lapack, *buff_U_lapack, *buff_V_lapack;
  int lwork, info, is_flame;
  
  if ( argc == 3 ) {
    m = atoi(argv[1]);
    is_flame = atoi(argv[2]);
  } else {
    fprintf(stderr, "       \n");
    fprintf(stderr, "Usage: %s m is_flame\n", argv[0]);
    fprintf(stderr, "       m : matrix length\n");
    fprintf(stderr, "       is_flame : 1 yes, 0 no\n");
    fprintf(stderr, "       \n");
    return -1;
  }
  if ( m == 0 )
    return 0;

  FLA_Init_safe( &init_result );          

  fprintf( stdout, "lapack2flame: %d x %d: \n", m, m);

  FLA_Obj_create( datatype, m, m, 0, 0, &A );
  FLA_Random_matrix( A ); 
  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_flame  );
  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_lapack );

  if ( is_flame ) {
    fprintf( stdout, " flame executed\n");
    FLA_Bidiag_UT_create_T( A_flame, &TU, &TV );

    FLA_Bidiag_UT( A_flame, TU, TV );
    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A_flame, &U_flame );
    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A_flame, &V_flame );

    FLA_Bidiag_UT_form_U( U_flame, TU, U_flame );
    FLA_Bidiag_UT_form_V( V_flame, TV, V_flame );
    
    FLA_Obj_create( datatype, m,      1, 0, 0, &d_flame );
    FLA_Obj_create( datatype, m - 1,  1, 0, 0, &e_flame );
    FLA_Bidiag_UT_extract_diagonals( A_flame, d_flame, e_flame );

    FLA_Obj_create( datatype, m, m, 0, 0, &B_flame ); FLA_Set( FLA_ZERO, B_flame );
    {
      FLA_Obj BTL, BTR, BBL, BBR;
      FLA_Part_2x2( B_flame, &BTL, &BTR, &BBL, &BBR, 1,1, FLA_BL );
      FLA_Set_diagonal_matrix( d_flame, B_flame );
      FLA_Set_diagonal_matrix( e_flame, BTR );
    }

    if (1) {
      fprintf( stdout, " - FLAME ----------\n");
      FLA_Obj_fshow( stdout, " - Given A - ", A, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - A - ", A_flame, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - U - ", U_flame, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - V - ", V_flame, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - d - ", d_flame, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - e - ", e_flame, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - B - ", B_flame, "% 6.4e", "------");
    }
  } else {
    fprintf( stdout, " lapack executed\n");

    FLA_Obj_create( datatype, m, 1, 0, 0, &tauq );
    FLA_Obj_create( datatype, m, 1, 0, 0, &taup );
    FLA_Obj_create( datatype, m,      1, 0, 0, &d_lapack );
    FLA_Obj_create( datatype, m - 1,  1, 0, 0, &e_lapack );

    buff_A_lapack = (testtype*)FLA_Obj_buffer_at_view( A_lapack );
    buff_tauq     = (testtype*)FLA_Obj_buffer_at_view( tauq );
    buff_taup     = (testtype*)FLA_Obj_buffer_at_view( taup );
    buff_d_lapack = (testtype*)FLA_Obj_buffer_at_view( d_lapack );
    buff_e_lapack = (testtype*)FLA_Obj_buffer_at_view( e_lapack );

    lwork = 32*m;
    
    FLA_Obj_create( datatype, lwork, 1, 0, 0, &W );
    buff_W = (testtype*)FLA_Obj_buffer_at_view( W );
    FORT_FUNC(s, gebrd_)( &m, &m, 
                          buff_A_lapack, &m,
                          buff_d_lapack,
                          buff_e_lapack,
                          buff_tauq,
                          buff_taup,
                          buff_W,
                          &lwork,
                          &info );
    
    FLA_Obj_create( datatype, m, m, 0, 0, &U_lapack );
    FLA_Obj_create( datatype, m, m, 0, 0, &V_lapack );
    
    FLA_Copy( A_lapack, U_lapack );
    FLA_Copy( A_lapack, V_lapack );

    buff_U_lapack = (testtype*)FLA_Obj_buffer_at_view( U_lapack );
    buff_V_lapack = (testtype*)FLA_Obj_buffer_at_view( V_lapack );

    FORT_FUNC(s, orgbr_)( "Q", &m, &m, &m,
                          buff_U_lapack, &m,
                          buff_tauq, 
                          buff_W,
                          &lwork,
                          &info );
    
    FORT_FUNC(s, orgbr_)( "P", &m, &m, &m,
                          buff_V_lapack, &m,
                          buff_taup,
                          buff_W,
                          &lwork,
                          &info );
    
    FLA_Obj_create( datatype, m, m, 0, 0, &B_lapack ); FLA_Set( FLA_ZERO, B_lapack );
    {
      FLA_Obj BTL, BTR, BBL, BBR;
      FLA_Part_2x2( B_lapack, &BTL, &BTR, &BBL, &BBR, 1,1, FLA_BL );
      FLA_Set_diagonal_matrix( d_lapack, B_lapack );
      FLA_Set_diagonal_matrix( e_lapack, BTR );
    }
    
    FLA_Obj_free( &W );    


    if (1) {
      fprintf( stdout, " - LAPACK ----------\n");
      FLA_Obj_fshow( stdout, " - Given A - ", A, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - A - ", A_lapack, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - U - ", U_lapack, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - V - ", V_lapack, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - d - ", d_lapack, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - e - ", e_lapack, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - B - ", B_lapack, "% 6.4e", "------");
    }
  }

  {
    testtype     dummy;
    int          zero = 0, one = 1;
    FLA_Obj      D_lapack;

    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &D_lapack ); FLA_Set( FLA_ZERO, D_lapack );

    if ( is_flame ) {
      buff_d_lapack = (testtype*)FLA_Obj_buffer_at_view( d_flame );
      buff_e_lapack = (testtype*)FLA_Obj_buffer_at_view( e_flame );
      buff_U_lapack = (testtype*)FLA_Obj_buffer_at_view( U_flame );
      buff_V_lapack = (testtype*)FLA_Obj_buffer_at_view( V_flame );
    }

    FLA_Obj_create( datatype, 4*m, 1, 0, 0, &W );
    buff_W = (testtype*)FLA_Obj_buffer_at_view( W );
    FORT_FUNC(s,bdsqr_)( "U", &m, &m, &m, &zero, 
                         buff_d_lapack, buff_e_lapack, 
                         buff_V_lapack, &m, 
                         buff_U_lapack, &m, 
                         &dummy, &one, 
                         buff_W, &info );
    FLA_Obj_free( &W );
    if (info != 0)
      printf( " Error info = %d\n", info );

    if ( is_flame )
      FLA_Set_diagonal_matrix( d_flame, D_lapack );
    else
      FLA_Set_diagonal_matrix( d_lapack, D_lapack );

    if ( is_flame ) {
      fprintf( stdout, " - FLAME ----------\n");
      FLA_Obj_fshow( stdout, " - U - ", U_flame, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - V - ", V_flame, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - d - ", d_flame, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - e - ", e_flame, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - D - ", D_lapack, "% 6.4e", "------");
    } else {
      fprintf( stdout, " - LAPACK ----------\n");
      FLA_Obj_fshow( stdout, " - U - ", U_lapack, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - V - ", V_lapack, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - d - ", d_lapack, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - e - ", e_lapack, "% 6.4e", "------");
      FLA_Obj_fshow( stdout, " - D - ", D_lapack, "% 6.4e", "------");
    }

    FLA_Obj_free( &D_lapack );
  }

  if ( is_flame ) {
    FLA_Obj_free( &TU );
    FLA_Obj_free( &TV );
    FLA_Obj_free( &U_flame );
    FLA_Obj_free( &V_flame );
    FLA_Obj_free( &d_flame );
    FLA_Obj_free( &e_flame );
    FLA_Obj_free( &B_flame );
  } else {
    FLA_Obj_free( &tauq );
    FLA_Obj_free( &taup );
    FLA_Obj_free( &d_lapack );
    FLA_Obj_free( &e_lapack );
    FLA_Obj_free( &U_lapack );
    FLA_Obj_free( &V_lapack );
    FLA_Obj_free( &B_lapack );
  }
  FLA_Obj_free( &A );
  FLA_Obj_free( &A_flame );
  FLA_Obj_free( &A_lapack );

  FLA_Finalize_safe( init_result );     
}
