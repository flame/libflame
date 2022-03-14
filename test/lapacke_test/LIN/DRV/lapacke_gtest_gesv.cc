#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"
#include <math.h>
#include "../lapacke_gtest_lin_common.hh"

TEST(gesv, dgesv1) {
   /* LAPACKE DGESV prototype */
typedef int (*Fptr_NL_LAPACKE_dgesv)( int matrix_layout, int n, int nrhs,
              double *a, int lda, int *ipiv, double *b, int ldb  );
                           
   void *hModule, *dModule;
   Fptr_NL_LAPACKE_dgesv DGESV;     
   double diff;
   int ipiv_diff;
   
   lin_solver_double_parameters dgesv_obj(LAPACK_COL_MAJOR, 'L',150, 40, 150, 150);

   //cgesv_double_parameters dgesv_obj(LAPACK_COL_MAJOR, 500, 40, 500, 500);
                                      
   dModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
  
   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
  
   DGESV = (Fptr_NL_LAPACKE_dgesv)dlsym(hModule, "LAPACKE_dgesv");

   if (NULL == DGESV)
   {
      printf("Could not get the DGESV symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   
   /* Compute the reference o/p by invoking Netlib-Lapack's API */
#if 0   
   dgesv_obj.inforef = DGESV( dgesv_obj.matrix_layout, 5,
                               1, dgesv_obj.aref, 5,
                 dgesv_obj.ipivref, dgesv_obj.bref, 5);
#endif
#if 1                            
   /* Compute the reference o/p by invoking Netlib-Lapack's API */  
   dgesv_obj.inforef = DGESV( dgesv_obj.matrix_layout, dgesv_obj.n,
                             dgesv_obj.nrhs, dgesv_obj.aref, dgesv_obj.lda,
                 dgesv_obj.ipivref, dgesv_obj.bref, dgesv_obj.ldb
                             );
   /* Compute libflame's Lapacke o/p  */
   dgesv_obj.info = LAPACKE_dgesv( dgesv_obj.matrix_layout ,dgesv_obj.n,
                            dgesv_obj.nrhs, dgesv_obj.a, dgesv_obj.lda,
                dgesv_obj.ipiv, dgesv_obj.b, dgesv_obj.ldb
                             );
#endif   

    if( dgesv_obj.info < 0 ) {
        printf( "The i:%d th argument with LibFlame is wrong\n", dgesv_obj.info );
        //exit( -1 );
    }
    if( dgesv_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LapackE is wrong\n", dgesv_obj.inforef );
        //exit( -1 );
    }

  /* info check validation is not needed here  */
    ipiv_diff = computeDiff_i( dgesv_obj.n, dgesv_obj.ipiv, dgesv_obj.ipivref );
    if (ipiv_diff >0){
        printf("\n dgesv Test failed wrt ipiv info");
    }
        
    diff =  computeDiff_d( (dgesv_obj.ldb*dgesv_obj.nrhs), dgesv_obj.b, dgesv_obj.bref );
    EXPECT_NEAR(0.0, diff, DOUBLE_DIFF_THRESHOLD);
}

/** test case for float version **/
TEST(gesv, sgesv1) {
   /* LAPACKE SGESV prototype */
typedef int (*Fptr_NL_LAPACKE_sgesv)( int matrix_layout, int n, int nrhs,
              float *a, int lda, int *ipiv, float *b, int ldb  );
                           
   void *hModule, *dModule;
   Fptr_NL_LAPACKE_sgesv SGESV;     
   float diff;
   int ipiv_diff;
   
   lin_solver_float_parameters sgesv_obj(LAPACK_COL_MAJOR, 'L',500, 40, 500, 500);

   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
  
   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
  
   SGESV = (Fptr_NL_LAPACKE_sgesv)dlsym(hModule, "LAPACKE_sgesv");

   if (NULL == SGESV)
   {
      printf("Could not get the SGESV symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   
   /* Compute the reference o/p by invoking Netlib-Lapack's API */
   sgesv_obj.inforef = SGESV( sgesv_obj.matrix_layout, sgesv_obj.n,
                             sgesv_obj.nrhs, sgesv_obj.aref, sgesv_obj.lda,
                 sgesv_obj.ipivref, sgesv_obj.bref, sgesv_obj.ldb
                             );
   /* Compute libflame's Lapacke o/p  */
   sgesv_obj.info = LAPACKE_sgesv( sgesv_obj.matrix_layout ,sgesv_obj.n,
                            sgesv_obj.nrhs, sgesv_obj.a, sgesv_obj.lda,
                sgesv_obj.ipiv, sgesv_obj.b, sgesv_obj.ldb
                             );

    if( sgesv_obj.info < 0 ) {
        printf( "The i:%d th argument with LibFlame is wrong\n", sgesv_obj.info );
        //exit( -1 );
    }
    if( sgesv_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LapackE is wrong\n", sgesv_obj.inforef );
        //exit( -1 );
    }

  /* info check validation is not needed here  */
    ipiv_diff = computeDiff_i( sgesv_obj.n, sgesv_obj.ipiv, sgesv_obj.ipivref );
    if (ipiv_diff >0){
        printf("\n sgesv Test failed wrt ipiv info");
    }
        
    diff =  computeDiff_s( (sgesv_obj.ldb*sgesv_obj.nrhs), sgesv_obj.b, sgesv_obj.bref );
    EXPECT_NEAR(0.0, diff, FLOAT_DIFF_THRESHOLD);
}


TEST(gesv, cgesv1) {
   /* LAPACKE CGESV prototype */
typedef int (*Fptr_NL_LAPACKE_cgesv)( int matrix_layout, int n, int nrhs,
              lapack_complex_float *a, int lda, int *ipiv, 
              lapack_complex_float *b, int ldb  );
                           
   void *hModule, *dModule;
   Fptr_NL_LAPACKE_cgesv CGESV;     
   float diff;
   int ipiv_diff;
   

   lin_solver_scomplex_parameters cgesv_obj(LAPACK_COL_MAJOR, 'L',150, 40, 150, 150);
                                      
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
  
   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
  
   CGESV = (Fptr_NL_LAPACKE_cgesv)dlsym(hModule, "LAPACKE_cgesv");

   if (NULL == CGESV)
   {
      printf("Could not get the CGESV symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   
   /* Compute the reference o/p by invoking Netlib-Lapack's API */  
   cgesv_obj.inforef = CGESV( cgesv_obj.matrix_layout, cgesv_obj.n,
                             cgesv_obj.nrhs, cgesv_obj.aref, cgesv_obj.lda,
                 cgesv_obj.ipivref, cgesv_obj.bref, cgesv_obj.ldb
                             );
   /* Compute libflame's Lapacke o/p  */
   cgesv_obj.info = LAPACKE_cgesv( cgesv_obj.matrix_layout ,cgesv_obj.n,
                            cgesv_obj.nrhs, cgesv_obj.a, cgesv_obj.lda,
                cgesv_obj.ipiv, cgesv_obj.b, cgesv_obj.ldb
                             );

    if( cgesv_obj.info < 0 ) {
        printf( "The i:%d th argument with LibFlame is wrong\n", cgesv_obj.info );
        //exit( -1 );
    }
    if( cgesv_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LapackE is wrong\n", cgesv_obj.inforef );
        //exit( -1 );
    }

  /* info check validation is not needed here  */
    ipiv_diff = computeDiff_i( cgesv_obj.n, cgesv_obj.ipiv, cgesv_obj.ipivref );
    if (ipiv_diff >0){
        printf("\n sgesv Test failed wrt ipiv info");
    }
        
    diff =  computeDiff_c( (cgesv_obj.ldb*cgesv_obj.nrhs), cgesv_obj.b, cgesv_obj.bref );
    EXPECT_NEAR(0.0, diff, FLOAT_DIFF_THRESHOLD);
}

TEST(gesv, zgesv1) {
   /* LAPACKE ZGESV prototype */
typedef int (*Fptr_NL_LAPACKE_zgesv)( int matrix_layout, int n, int nrhs,
              lapack_complex_double *a, int lda, int *ipiv, 
              lapack_complex_double *b, int ldb  );
                           
   void *hModule, *dModule;
   Fptr_NL_LAPACKE_zgesv ZGESV;     
   double diff;
   int ipiv_diff;
   

   lin_solver_dcomplex_parameters zgesv_obj(LAPACK_COL_MAJOR, 'L',250, 100, 250, 250);
                                      
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
  
   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
  
   ZGESV = (Fptr_NL_LAPACKE_zgesv)dlsym(hModule, "LAPACKE_zgesv");

   if (NULL == ZGESV)
   {
      printf("Could not get the ZGESV symbol. Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   
   /* Compute the reference o/p by invoking Netlib-Lapack's API */  
   zgesv_obj.inforef = ZGESV( zgesv_obj.matrix_layout, zgesv_obj.n,
                             zgesv_obj.nrhs, zgesv_obj.aref, zgesv_obj.lda,
                             zgesv_obj.ipivref, zgesv_obj.bref, zgesv_obj.ldb
                             );
   /* Compute libflame's Lapacke o/p  */
   zgesv_obj.info = LAPACKE_zgesv( zgesv_obj.matrix_layout ,zgesv_obj.n,
                            zgesv_obj.nrhs, zgesv_obj.a, zgesv_obj.lda,
                            zgesv_obj.ipiv, zgesv_obj.b, zgesv_obj.ldb
                             );

    if( zgesv_obj.info < 0 ) {
        printf( "The i:%d th argument with LibFlame is wrong\n", zgesv_obj.info );
        //exit( -1 );
    }
    if( zgesv_obj.inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LapackE is wrong\n", zgesv_obj.inforef );
        //exit( -1 );
    }

  /* info check validation is not needed here  */
    ipiv_diff = computeDiff_i( zgesv_obj.n, zgesv_obj.ipiv, zgesv_obj.ipivref );
    if (ipiv_diff >0){
        printf("\n sgesv Test failed wrt ipiv info");
    }
        
    diff =  computeDiff_z( (zgesv_obj.ldb*zgesv_obj.nrhs), zgesv_obj.b, zgesv_obj.bref );
    EXPECT_NEAR(0.0, diff, FLOAT_DIFF_THRESHOLD);
}
