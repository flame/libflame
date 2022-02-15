#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define LAPACKE_TEST_VERBOSE (1)
#define getrs_free() \
       if (ipiv != NULL) free (ipiv); \
       if (bref != NULL) free (bref); \
       if (b != NULL)    free (b   ); \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (ipivref != NULL)free (ipivref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin getrs_double_parameters  class definition */
class getrs_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      getrs_double_parameters ( int matrix_layout_i, char trans_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~getrs_double_parameters ();
};  /* end of getrs_double_parameters  class definition */


/* Constructor getrs_double_parameters definition */
getrs_double_parameters:: getrs_double_parameters ( int matrix_layout_i,
                         char trans_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    trans = trans_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n getrs Double:  n: %d, trans: %c lda: %d ldb: %d nrhs: %d \n",
             n, trans, lda, ldb, nrhs);
#endif


    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "getrs_double_parameters object: malloc error.";
       getrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

getrs_double_parameters:: ~getrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" getrs_double_parameters object: destructor invoked. \n");
#endif
   getrs_free();
}


//  Test fixture class definition
class dgetrs_test  : public  ::testing::Test {
public:
   getrs_double_parameters  *dgetrs_obj;
   void SetUp();  
   void TearDown () { delete dgetrs_obj; }
};


void dgetrs_test::SetUp(){

    /* LAPACKE DGETRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dgetrs) ( int matrix_layout, char trans,
                          lapack_int n, lapack_int nrhs, const double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                                            double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dgetrs DGETRS;

     /* LAPACKE DGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_dgetrf DGETRF;

    dgetrs_obj = new  getrs_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dgetrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgetrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgetrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgetrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGETRS = (Fptr_NL_LAPACKE_dgetrs)dlsym(dgetrs_obj->hModule, "LAPACKE_dgetrs");
    ASSERT_TRUE(DGETRS != NULL) << "failed to get the Netlib LAPACKE_dgetrs symbol";

    DGETRF = (Fptr_NL_LAPACKE_dgetrf)dlsym(dgetrs_obj->hModule,"LAPACKE_dgetrf");
    ASSERT_TRUE(DGETRF != NULL) << "failed to get the Netlib LAPACKE_dgetrf symbol";

    /* Pre condition: need to call getrf - before calling getrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    dgetrs_obj->inforef = DGETRF( dgetrs_obj->matrix_layout,
                                    dgetrs_obj->n, dgetrs_obj->n,
                                     dgetrs_obj->aref,
                               dgetrs_obj->lda, dgetrs_obj->ipivref);

    dgetrs_obj->inforef = DGETRS( dgetrs_obj->matrix_layout,
                                  dgetrs_obj->trans, dgetrs_obj->n,
                                  dgetrs_obj->nrhs,
                                  (const double *)dgetrs_obj->aref,
                                  dgetrs_obj->lda, dgetrs_obj->ipivref,
                                  dgetrs_obj->bref, dgetrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dgetrs_obj->info = LAPACKE_dgetrf( dgetrs_obj->matrix_layout,
                                    dgetrs_obj->n, dgetrs_obj->n,
                                     dgetrs_obj->a,
                               dgetrs_obj->lda, dgetrs_obj->ipiv);

    dgetrs_obj->info = LAPACKE_dgetrs( dgetrs_obj->matrix_layout,
                dgetrs_obj->trans, dgetrs_obj->n, dgetrs_obj->nrhs,
                                  (const double *)dgetrs_obj->a,
                               dgetrs_obj->lda, dgetrs_obj->ipiv,
                                 dgetrs_obj->b, dgetrs_obj->ldb );


    if( dgetrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgetrs is wrong\n",
                    dgetrs_obj->info );
    }
    if( dgetrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgetrs is wrong\n",
        dgetrs_obj->inforef );
    }
}

TEST_F(dgetrs_test, dgetrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgetrs_obj->b_bufsize,
                           dgetrs_obj->b, dgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dgetrs_test, dgetrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgetrs_obj->b_bufsize,
                           dgetrs_obj->b, dgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dgetrs_test, dgetrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgetrs_obj->b_bufsize,
                           dgetrs_obj->b, dgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dgetrs_test, dgetrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dgetrs_obj->b_bufsize,
                           dgetrs_obj->b, dgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin getrs_float_parameters  class definition */
class getrs_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      getrs_float_parameters (int matrix_layout_i, char trans_i,
                               lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~getrs_float_parameters ();
};  /* end of getrs_float_parameters  class definition */


/* Constructor getrs_float_parameters definition */
getrs_float_parameters:: getrs_float_parameters ( int matrix_layout_i,
                       char trans_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    trans = trans_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n getrs float:  n: %d, trans: %c lda: %d ldb: %d nrhs: %d \n",
             n, trans, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (lda*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "getrs_double_parameters object: malloc error.";
       getrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

getrs_float_parameters:: ~getrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" getrs_float_parameters object: destructor invoked. \n");
#endif
   getrs_free();
}


//  Test fixture class definition
class sgetrs_test  : public  ::testing::Test {
public:
   getrs_float_parameters  *sgetrs_obj;
   void SetUp();  
   void TearDown () { delete sgetrs_obj; }
};


void sgetrs_test::SetUp(){

    /* LAPACKE SGETRS prototype */
    typedef int (*Fptr_NL_LAPACKE_sgetrs) ( int matrix_layout, char trans,
                          lapack_int n, lapack_int nrhs, const float * a,
                                  lapack_int lda, const lapack_int * ipiv,
                                            float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_sgetrs SGETRS;

     /* LAPACKE SGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_sgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_sgetrf SGETRF;

    sgetrs_obj = new  getrs_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    sgetrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgetrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgetrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgetrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGETRS = (Fptr_NL_LAPACKE_sgetrs)dlsym(sgetrs_obj->hModule, "LAPACKE_sgetrs");
    ASSERT_TRUE(SGETRS != NULL) << "failed to get the Netlib LAPACKE_sgetrs symbol";

    SGETRF = (Fptr_NL_LAPACKE_sgetrf)dlsym(sgetrs_obj->hModule,"LAPACKE_sgetrf");
    ASSERT_TRUE(SGETRF != NULL) << "failed to get the Netlib LAPACKE_sgetrf symbol";

    /* Pre condition: need to call getrf - before calling getrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    sgetrs_obj->inforef = SGETRF( sgetrs_obj->matrix_layout,
                                    sgetrs_obj->n, sgetrs_obj->n,
                                     sgetrs_obj->aref,
                               sgetrs_obj->lda, sgetrs_obj->ipivref);

    sgetrs_obj->inforef = SGETRS( sgetrs_obj->matrix_layout,
                                  sgetrs_obj->trans, sgetrs_obj->n,
                                  sgetrs_obj->nrhs,
                                  (const float *)sgetrs_obj->aref,
                                  sgetrs_obj->lda, sgetrs_obj->ipivref,
                                  sgetrs_obj->bref, sgetrs_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    sgetrs_obj->info = LAPACKE_sgetrf( sgetrs_obj->matrix_layout,
                                    sgetrs_obj->n, sgetrs_obj->n,
                                     sgetrs_obj->a,
                               sgetrs_obj->lda, sgetrs_obj->ipiv);

    sgetrs_obj->info = LAPACKE_sgetrs( sgetrs_obj->matrix_layout,
                sgetrs_obj->trans, sgetrs_obj->n, sgetrs_obj->nrhs,
                                  (const float *)sgetrs_obj->a,
                               sgetrs_obj->lda, sgetrs_obj->ipiv,
                                 sgetrs_obj->b, sgetrs_obj->ldb );
    if( sgetrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgetrs is wrong\n",
                    sgetrs_obj->info );
    }
    if( sgetrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgetrs is wrong\n",
        sgetrs_obj->inforef );
    }
}

TEST_F(sgetrs_test, sgetrs1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sgetrs_obj->b_bufsize,
                           sgetrs_obj->b, sgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin getrs_scomplex_parameters  class definition */
class getrs_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      getrs_scomplex_parameters ( int matrix_layout_i, char trans_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~getrs_scomplex_parameters ();
};  /* end of getrs_scomplex_parameters  class definition */


/* Constructor getrs_scomplex_parameters definition */
getrs_scomplex_parameters:: getrs_scomplex_parameters ( int matrix_layout_i,
                         char trans_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    trans = trans_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n getrs scomplex:  n: %d, trans: %c lda: %d ldb: %d nrhs: %d \n",
             n, trans, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "getrs_scomplex_parameters object: malloc error.";
       getrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

getrs_scomplex_parameters:: ~getrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" getrs_scomplex_parameters object: destructor invoked. \n");
#endif
   getrs_free();
}


//  Test fixture class definition
class cgetrs_test  : public  ::testing::Test {
public:
   getrs_scomplex_parameters  *cgetrs_obj;
   void SetUp();  
   void TearDown () { delete cgetrs_obj; }
};


void cgetrs_test::SetUp(){

    /* LAPACKE CGETRS prototype */
    typedef int (*Fptr_NL_LAPACKE_cgetrs) ( int matrix_layout, char trans,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_float * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_cgetrs CGETRS;

     /* LAPACKE CGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_cgetrf CGETRF;


    cgetrs_obj = new  getrs_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    cgetrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgetrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgetrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgetrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGETRS = (Fptr_NL_LAPACKE_cgetrs)dlsym(cgetrs_obj->hModule, "LAPACKE_cgetrs");
    ASSERT_TRUE(CGETRS != NULL) << "failed to get the Netlib LAPACKE_cgetrs symbol";

    CGETRF = (Fptr_NL_LAPACKE_cgetrf)dlsym(cgetrs_obj->hModule,"LAPACKE_cgetrf");
    ASSERT_TRUE(CGETRF != NULL) << "failed to get the Netlib LAPACKE_cgetrf symbol";

    /* Pre condition: need to call getrf - before calling getrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    cgetrs_obj->inforef = CGETRF( cgetrs_obj->matrix_layout,
                                    cgetrs_obj->n, cgetrs_obj->n,
                                     cgetrs_obj->aref,
                               cgetrs_obj->lda, cgetrs_obj->ipivref);

    cgetrs_obj->inforef = CGETRS( cgetrs_obj->matrix_layout,
                                  cgetrs_obj->trans, cgetrs_obj->n,
                                  cgetrs_obj->nrhs,
                                  (const lapack_complex_float *)cgetrs_obj->aref,
                                  cgetrs_obj->lda, cgetrs_obj->ipivref,
                                  cgetrs_obj->bref, cgetrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    cgetrs_obj->info = LAPACKE_cgetrf( cgetrs_obj->matrix_layout,
                                    cgetrs_obj->n, cgetrs_obj->n,
                                     cgetrs_obj->a,
                               cgetrs_obj->lda, cgetrs_obj->ipiv);

    cgetrs_obj->info = LAPACKE_cgetrs( cgetrs_obj->matrix_layout,
                cgetrs_obj->trans, cgetrs_obj->n, cgetrs_obj->nrhs,
                                  (const lapack_complex_float *)cgetrs_obj->a,
                               cgetrs_obj->lda, cgetrs_obj->ipiv,
                                 cgetrs_obj->b, cgetrs_obj->ldb );


    if( cgetrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgetrs is wrong\n",
                    cgetrs_obj->info );
    }
    if( cgetrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgetrs is wrong\n",
        cgetrs_obj->inforef );
    }
}

TEST_F(cgetrs_test, cgetrs1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgetrs_obj->b_bufsize,
                           cgetrs_obj->b, cgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cgetrs_test, cgetrs2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgetrs_obj->b_bufsize,
                           cgetrs_obj->b, cgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cgetrs_test, cgetrs3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgetrs_obj->b_bufsize,
                           cgetrs_obj->b, cgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cgetrs_test, cgetrs4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cgetrs_obj->b_bufsize,
                           cgetrs_obj->b, cgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin getrs_dcomplex_parameters  class definition */
class getrs_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      getrs_dcomplex_parameters ( int matrix_layout_i, char trans_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~getrs_dcomplex_parameters ();
};  /* end of getrs_dcomplex_parameters  class definition */


/* Constructor getrs_dcomplex_parameters definition */
getrs_dcomplex_parameters:: getrs_dcomplex_parameters ( int matrix_layout_i,
                         char trans_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    trans = trans_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n getrs DComplex:  n: %d, trans: %c lda: %d ldb: %d nrhs: %d \n",
             n, trans, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "getrs_dcomplex_parameters object: malloc error.";
       getrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

getrs_dcomplex_parameters:: ~getrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" getrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   getrs_free();
}


//  Test fixture class definition
class zgetrs_test  : public  ::testing::Test {
public:
   getrs_dcomplex_parameters  *zgetrs_obj;
   void SetUp();  
   void TearDown () { delete zgetrs_obj; }
};


void zgetrs_test::SetUp(){

    /* LAPACKE ZGETRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zgetrs) ( int matrix_layout, char trans,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zgetrs ZGETRS;

     /* LAPACKE ZGETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zgetrf) ( int matrix_layout,lapack_int m,lapack_int n,
                                    lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zgetrf ZGETRF;


    zgetrs_obj = new  getrs_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zgetrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgetrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgetrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgetrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGETRS = (Fptr_NL_LAPACKE_zgetrs)dlsym(zgetrs_obj->hModule, "LAPACKE_zgetrs");
    ASSERT_TRUE(ZGETRS != NULL) << "failed to get the Netlib LAPACKE_zgetrs symbol";

    ZGETRF = (Fptr_NL_LAPACKE_zgetrf)dlsym(zgetrs_obj->hModule,"LAPACKE_zgetrf");
    ASSERT_TRUE(ZGETRF != NULL) << "failed to get the Netlib LAPACKE_zgetrf symbol";

    /* Pre condition: need to call getrf - before calling getrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zgetrs_obj->inforef = ZGETRF( zgetrs_obj->matrix_layout,
                                    zgetrs_obj->n, zgetrs_obj->n,
                                     zgetrs_obj->aref,
                               zgetrs_obj->lda, zgetrs_obj->ipivref);

    zgetrs_obj->inforef = ZGETRS( zgetrs_obj->matrix_layout,
                                  zgetrs_obj->trans, zgetrs_obj->n,
                                  zgetrs_obj->nrhs,
                                  (const lapack_complex_double *)zgetrs_obj->aref,
                                  zgetrs_obj->lda, zgetrs_obj->ipivref,
                                  zgetrs_obj->bref, zgetrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zgetrs_obj->info = LAPACKE_zgetrf( zgetrs_obj->matrix_layout,
                                    zgetrs_obj->n, zgetrs_obj->n,
                                     zgetrs_obj->a,
                               zgetrs_obj->lda, zgetrs_obj->ipiv);

    zgetrs_obj->info = LAPACKE_zgetrs( zgetrs_obj->matrix_layout,
                zgetrs_obj->trans, zgetrs_obj->n, zgetrs_obj->nrhs,
                                  (const lapack_complex_double *)zgetrs_obj->a,
                               zgetrs_obj->lda, zgetrs_obj->ipiv,
                                 zgetrs_obj->b, zgetrs_obj->ldb );


    if( zgetrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgetrs is wrong\n",
                    zgetrs_obj->info );
    }
    if( zgetrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgetrs is wrong\n",
        zgetrs_obj->inforef );
    }
}

TEST_F(zgetrs_test, zgetrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgetrs_obj->b_bufsize,
                           zgetrs_obj->b, zgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zgetrs_test, zgetrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgetrs_obj->b_bufsize,
                           zgetrs_obj->b, zgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zgetrs_test, zgetrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgetrs_obj->b_bufsize,
                           zgetrs_obj->b, zgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zgetrs_test, zgetrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zgetrs_obj->b_bufsize,
                           zgetrs_obj->b, zgetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
