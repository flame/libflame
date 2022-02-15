#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

static int idx = 0; // parameter index of the parameters structure array

//  Test fixture class definition
class dgetrf_test  : public  ::testing::Test {
public:
   double_common_parameters  *dgetrf_obj;
   void SetUp();   
   void TearDown () { delete dgetrf_obj; }
};

void dgetrf_test::SetUp(){

    dgetrf_obj = new  double_common_parameters(lin_solver_paramslist[idx].m,
                                              lin_solver_paramslist[idx].n);
    idx = Circular_Increment_Index(idx);

    /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgetrf) ( int matrix_layout,lapack_int m,
                 lapack_int n, double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_dgetrf DGETRF;

    dgetrf_obj->dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    dgetrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == dgetrf_obj->hModule) || (NULL == dgetrf_obj->dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit (-1);
    }

    DGETRF = (Fptr_NL_LAPACKE_dgetrf)dlsym(dgetrf_obj->hModule,"LAPACKE_dgetrf");
    if (NULL == DGETRF)
    {
      printf("Could not get the symbol. Exiting...\n");
      dlclose(dgetrf_obj->hModule);
      dlclose(dgetrf_obj->dModule);
      exit (-1);
    }
    printf(" dgetrf_obj->m: %d ,dgetrf_obj->n : %d \n", dgetrf_obj->m,dgetrf_obj->n);
       
    /* Compute libflame's Lapacke o/p  */
    dgetrf_obj->info  = LAPACKE_dgetrf( LAPACK_COL_MAJOR,dgetrf_obj->m,dgetrf_obj->n,
                                    dgetrf_obj->A, dgetrf_obj->lda,dgetrf_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    dgetrf_obj->inforef = DGETRF( LAPACK_COL_MAJOR,dgetrf_obj->m,dgetrf_obj->n,
                         dgetrf_obj->Aref, dgetrf_obj->lda,dgetrf_obj->ipivref);

    if( dgetrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgetrf \
                                                is wrong\n", dgetrf_obj->info );
    }
    if( dgetrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgetrf is wrong\n",
        dgetrf_obj->inforef );
    }
    dlclose(dgetrf_obj->hModule);
    dlclose(dgetrf_obj->dModule);
}

TEST_F(dgetrf_test, dgetrf1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgetrf_obj->m*dgetrf_obj->n,dgetrf_obj->A,
                                                  dgetrf_obj->Aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgetrf_test, dgetrf2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgetrf_obj->m*dgetrf_obj->n,dgetrf_obj->A,
                                                  dgetrf_obj->Aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgetrf_test, dgetrf3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgetrf_obj->m*dgetrf_obj->n,dgetrf_obj->A,
                                                  dgetrf_obj->Aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgetrf_test, dgetrf4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgetrf_obj->m*dgetrf_obj->n,dgetrf_obj->A,
                                                  dgetrf_obj->Aref );

    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

class sgetrf_test  : public  ::testing::Test {
public:
   float_common_parameters  *sgetrf_obj;
   void SetUp();   
   void TearDown () { delete sgetrf_obj; }
};


void sgetrf_test::SetUp(){

    sgetrf_obj = new  float_common_parameters(lin_solver_paramslist[idx].m,
                                              lin_solver_paramslist[idx].n);
    idx = Circular_Increment_Index(idx);

    /* LAPACKE sgetrf prototype */
    typedef int (*Fptr_NL_LAPACKE_sgetrf) ( int matrix_layout,lapack_int m,
                  lapack_int n, float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_sgetrf SGETRF;


    sgetrf_obj->dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    sgetrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == sgetrf_obj->hModule) || (NULL == sgetrf_obj->dModule) )
    {
            printf("Load Library failed. Exiting ....\n");
            exit( -1 );
    }

    SGETRF = (Fptr_NL_LAPACKE_sgetrf)dlsym(sgetrf_obj->hModule,"LAPACKE_sgetrf");
    if (NULL == SGETRF)
    {
          printf("Could not get the symbol. Exiting...\n");
          dlclose(sgetrf_obj->hModule);
          dlclose(sgetrf_obj->dModule);
          exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    sgetrf_obj->info     = LAPACKE_sgetrf( LAPACK_COL_MAJOR,sgetrf_obj->m,
           sgetrf_obj->n,sgetrf_obj->A, sgetrf_obj->lda,sgetrf_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    sgetrf_obj->inforef = SGETRF( LAPACK_COL_MAJOR,sgetrf_obj->m,
      sgetrf_obj->n,sgetrf_obj->Aref, sgetrf_obj->lda,sgetrf_obj->ipivref);

    if( sgetrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgetrf \
                                               is wrong\n", sgetrf_obj->info );
    }
    if( sgetrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgetrf is wrong\n",
        sgetrf_obj->inforef );
    }
    dlclose(sgetrf_obj->hModule);
    dlclose(sgetrf_obj->dModule);
}

TEST_F(sgetrf_test, sgetrf1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sgetrf_obj->m*sgetrf_obj->n,sgetrf_obj->A,
                                                 sgetrf_obj->Aref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgetrf_test, sgetrf2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sgetrf_obj->m*sgetrf_obj->n,sgetrf_obj->A,
                                                  sgetrf_obj->Aref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

class cgetrf_test  : public  ::testing::Test {
public:
   scomplex_common_parameters  *cgetrf_obj;
   void SetUp();   
   void TearDown () { delete cgetrf_obj; }
};

void cgetrf_test::SetUp(){

    cgetrf_obj = new  scomplex_common_parameters(lin_solver_paramslist[idx].m,
                                              lin_solver_paramslist[idx].n);
    idx = Circular_Increment_Index(idx);

    /* LAPACKE cgetrf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgetrf) ( int matrix_layout,lapack_int m,
     lapack_int n, lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_cgetrf cgetrf;

    cgetrf_obj->dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    cgetrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == cgetrf_obj->hModule) || (NULL == cgetrf_obj->dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    cgetrf = (Fptr_NL_LAPACKE_cgetrf)dlsym(cgetrf_obj->hModule,"LAPACKE_cgetrf");
    if (NULL == cgetrf)
    {
        printf("Could not get the symbol. Exiting...\n");
        dlclose(cgetrf_obj->hModule);
        dlclose(cgetrf_obj->dModule);
        exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    cgetrf_obj->info     = LAPACKE_cgetrf( LAPACK_COL_MAJOR,cgetrf_obj->m,
           cgetrf_obj->n,cgetrf_obj->A, cgetrf_obj->lda,cgetrf_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cgetrf_obj->inforef = cgetrf( LAPACK_COL_MAJOR,cgetrf_obj->m,cgetrf_obj->n,cgetrf_obj->Aref,
                                        cgetrf_obj->lda,cgetrf_obj->ipivref);

    if( cgetrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgetrf is wrong\n",
                    cgetrf_obj->info );
    }
    if( cgetrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgetrf is wrong\n",
        cgetrf_obj->inforef );
    }
    dlclose(cgetrf_obj->hModule);
    dlclose(cgetrf_obj->dModule);
}

TEST_F(cgetrf_test, cgetrf1) {
    /* Compute Difference between libflame and Netlib o/ps  */
    float diff =  computeDiff_c( cgetrf_obj->m*cgetrf_obj->n,cgetrf_obj->A,
                                                        cgetrf_obj->Aref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgetrf_test, cgetrf2) {
    /* Compute Difference between libflame and Netlib o/ps  */
    float diff =  computeDiff_c( cgetrf_obj->m*cgetrf_obj->n,
                             cgetrf_obj->A,cgetrf_obj->Aref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

class zgetrf_test  : public  ::testing::Test {
public:
   dcomplex_common_parameters  *zgetrf_obj;
   void SetUp();   
   void TearDown () { delete zgetrf_obj; }
};

void zgetrf_test::SetUp(){

    zgetrf_obj = new dcomplex_common_parameters(lin_solver_paramslist[idx].m,
                                              lin_solver_paramslist[idx].n);
    idx = Circular_Increment_Index(idx);
    /* LAPACKE zgetrf prototype */
    typedef int (*Fptr_NL_LAPACKE_zgetrf) (int matrix_layout,lapack_int m,
                                    lapack_int n, lapack_complex_double *a,
                                         lapack_int lda,lapack_int* ipiv );
    Fptr_NL_LAPACKE_zgetrf ZGETRF2;


    zgetrf_obj->dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
    zgetrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
    if ((NULL == zgetrf_obj->hModule) || (NULL == zgetrf_obj->dModule) )
    {
        printf("Load Library failed. Exiting ....\n");
        exit( -1 );
    }

    ZGETRF2 = (Fptr_NL_LAPACKE_zgetrf)dlsym(zgetrf_obj->hModule,"LAPACKE_zgetrf");
    if (NULL == ZGETRF2)
    {
        printf("Could not get the symbol. Exiting...\n");
        dlclose(zgetrf_obj->hModule);
        dlclose(zgetrf_obj->dModule);
        exit( -1 );
    }

    /* Compute libflame's Lapacke o/p  */
    zgetrf_obj->info     = LAPACKE_zgetrf( LAPACK_COL_MAJOR,zgetrf_obj->m,
            zgetrf_obj->n, zgetrf_obj->A,zgetrf_obj->lda,zgetrf_obj->ipiv);

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zgetrf_obj->inforef = ZGETRF2( LAPACK_COL_MAJOR,zgetrf_obj->m,
       zgetrf_obj->n, zgetrf_obj->Aref,zgetrf_obj->lda,zgetrf_obj->ipivref);

    if( zgetrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
                   LAPACKE_zgetrf is wrong\n", zgetrf_obj->info );
    }
    if( zgetrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgetrf is wrong\n",
        zgetrf_obj->inforef );
    }
    dlclose(zgetrf_obj->hModule);
    dlclose(zgetrf_obj->dModule);
}

TEST_F(zgetrf_test, zgetrf1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgetrf_obj->m*zgetrf_obj->n,
                      zgetrf_obj->A,zgetrf_obj->Aref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgetrf_test, zgetrf2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgetrf_obj->m*zgetrf_obj->n,
                      zgetrf_obj->A,zgetrf_obj->Aref );
    EXPECT_NEAR(0.0,diff,LAPACKE_GTEST_THRESHOLD);
}
