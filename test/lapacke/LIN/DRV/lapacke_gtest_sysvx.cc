#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define sysvx_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (x != NULL)    free (x  ); \
  if (xref != NULL) free (xref); \
  if (af != NULL)    free (af  ); \
  if (afref != NULL) free (afref); \
  if (ferr != NULL)    free (ferr  ); \
  if (ferrref != NULL) free (ferrref); \
  if (berr != NULL)    free (berr  ); \
  if (berrref != NULL) free (berrref); \
  if (ipiv != NULL)    free (ipiv  ); \
  if (ipivref != NULL) free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin sysvx_float_parameters  class definition */
class sysvx_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
      int ipiv_diff;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.
      float *a, *aref; //The array ab contains the matrix A
      float *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sysvx_float_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~sysvx_float_parameters (); 
};  /* end of sysvx_float_parameters  class definition */


/* Constructor sysvx_float_parameters definition */
sysvx_float_parameters:: sysvx_float_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n sysvx float:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);

    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       sysvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(float)));
    memcpy(afref, a, (n*n*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_float_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

sysvx_float_parameters:: ~sysvx_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysvx_float_parameters object: destructor invoked. \n");
#endif
   sysvx_free();
}

//  Test fixture class definition
class ssysvx_test  : public  ::testing::Test {
public:
   sysvx_float_parameters  *ssysvx_obj;
   void SetUp();  
   void TearDown () { delete ssysvx_obj; }
};


void ssysvx_test::SetUp(){

    /* LAPACKE SSYSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_ssysvx) (int matrix_layout, char fact,
            char uplo, lapack_int n, lapack_int nrhs, const float* a,
            lapack_int lda, float* af, lapack_int ldaf, lapack_int* ipiv,
            const float* b, lapack_int ldb, float* x, lapack_int ldx,
            float* rcond, float* ferr, float* berr);

    Fptr_NL_LAPACKE_ssysvx SSYSVX;

     /* LAPACKE SSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf) ( int matrix_layout , char uplo,
            lapack_int n , float * a , lapack_int lda , lapack_int * ipiv );

    Fptr_NL_LAPACKE_ssytrf SSYTRF;

    ssysvx_obj = new sysvx_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    ssysvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssysvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssysvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssysvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SSYSVX = (Fptr_NL_LAPACKE_ssysvx)dlsym(ssysvx_obj->hModule, "LAPACKE_ssysvx");
    ASSERT_TRUE(SSYSVX != NULL) << "failed to syt the Netlib LAPACKE_ssysvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the sytrf API to compute the factorized A.  */
    if(ssysvx_obj->fact == 'F') {
        SSYTRF = (Fptr_NL_LAPACKE_ssytrf)dlsym(ssysvx_obj->hModule,"LAPACKE_ssytrf");
        ASSERT_TRUE(SSYTRF != NULL) << "failed to syt the Netlib LAPACKE_ssytrf symbol";
            
        ssysvx_obj->inforef = SSYTRF( ssysvx_obj->matrix_layout,
                                      ssysvx_obj->uplo, ssysvx_obj->n,
                                      ssysvx_obj->afref,
                                      ssysvx_obj->lda,
                                      ssysvx_obj->ipivref);
                               
        ssysvx_obj->info = LAPACKE_ssytrf( ssysvx_obj->matrix_layout,
                                           ssysvx_obj->uplo, ssysvx_obj->n,
                                           ssysvx_obj->af,
                                           ssysvx_obj->lda,
                                           ssysvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    ssysvx_obj->inforef = SSYSVX( ssysvx_obj->matrix_layout, ssysvx_obj->fact,
                                  ssysvx_obj->uplo, ssysvx_obj->n,
                                  ssysvx_obj->nrhs,
                                  ssysvx_obj->aref, ssysvx_obj->lda, 
                                  ssysvx_obj->afref, ssysvx_obj->ldaf,
                                  ssysvx_obj->ipivref,
                                  ssysvx_obj->bref, ssysvx_obj->ldb,
                                  ssysvx_obj->xref, ssysvx_obj->ldx,
                                  &ssysvx_obj->rcondref, 
                                  ssysvx_obj->ferrref,
                                  ssysvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    ssysvx_obj->info = LAPACKE_ssysvx( ssysvx_obj->matrix_layout, ssysvx_obj->fact,
                                  ssysvx_obj->uplo, ssysvx_obj->n,
                                  ssysvx_obj->nrhs,
                                  ssysvx_obj->a, ssysvx_obj->lda, 
                                  ssysvx_obj->af, ssysvx_obj->ldaf,
                                  ssysvx_obj->ipiv,
                                  ssysvx_obj->b, ssysvx_obj->ldb,
                                  ssysvx_obj->x, ssysvx_obj->ldx,
                                  &ssysvx_obj->rcond, 
                                  ssysvx_obj->ferr,
                                  ssysvx_obj->berr);

    if( ssysvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssysvx is wrong\n", ssysvx_obj->info );
    }
    if( ssysvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssysvx is wrong\n", 
        ssysvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ssysvx_obj->ipiv_diff = computeDiff_i( ssysvx_obj->n, ssysvx_obj->ipiv, ssysvx_obj->ipivref );
    
    ssysvx_obj->diff =  computeDiff_s( ssysvx_obj->n * ssysvx_obj->n, 
                ssysvx_obj->af, ssysvx_obj->afref );

    ssysvx_obj->diff_xerr =  computeDiff_s( ssysvx_obj->x_bufsize, 
                ssysvx_obj->x, ssysvx_obj->xref );

    ssysvx_obj->diff_berr =  computeDiff_s( ssysvx_obj->nrhs, 
                ssysvx_obj->berr, ssysvx_obj->berrref );
                
    ssysvx_obj->diff_ferr =  computeDiff_s( ssysvx_obj->nrhs, 
                ssysvx_obj->ferr, ssysvx_obj->ferrref );
}

TEST_F(ssysvx_test, ssysvx1) {
    EXPECT_NEAR(0.0, ssysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (ssysvx_obj->rcond - ssysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, ssysvx_obj->ipiv_diff );
}

TEST_F(ssysvx_test, ssysvx2) {
    EXPECT_NEAR(0.0, ssysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (ssysvx_obj->rcond - ssysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysvx_test, ssysvx3) {
    EXPECT_NEAR(0.0, ssysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (ssysvx_obj->rcond - ssysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysvx_test, ssysvx4) {
    EXPECT_NEAR(0.0, ssysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, ssysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (ssysvx_obj->rcond - ssysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin sysvx_double_parameters  class definition */
class sysvx_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
      int ipiv_diff;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.
      double *a, *aref; //The array ab contains the matrix A
      double *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sysvx_double_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~sysvx_double_parameters (); 
};  /* end of sysvx_double_parameters  class definition */


/* Constructor sysvx_double_parameters definition */
sysvx_double_parameters:: sysvx_double_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n sysvx double:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);

    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       sysvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(double)));
    memcpy(afref, a, (n*n*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_double_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

sysvx_double_parameters:: ~sysvx_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysvx_double_parameters object: destructor invoked. \n");
#endif
   sysvx_free();
}

//  Test fixture class definition
class dsysvx_test  : public  ::testing::Test {
public:
   sysvx_double_parameters  *dsysvx_obj;
   void SetUp();  
   void TearDown () { delete dsysvx_obj; }
};


void dsysvx_test::SetUp(){

    /* LAPACKE DSYSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_dsysvx) (int matrix_layout, char fact,
            char uplo, lapack_int n, lapack_int nrhs, const double* a,
            lapack_int lda, double* af, lapack_int ldaf, lapack_int* ipiv,
            const double* b, lapack_int ldb, double* x, lapack_int ldx,
            double* rcond, double* ferr, double* berr);

    Fptr_NL_LAPACKE_dsysvx DSYSVX;

     /* LAPACKE DSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf) ( int matrix_layout , char uplo,
            lapack_int n , double * a , lapack_int lda , lapack_int * ipiv );

    Fptr_NL_LAPACKE_dsytrf DSYTRF;

    dsysvx_obj = new sysvx_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dsysvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsysvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsysvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsysvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DSYSVX = (Fptr_NL_LAPACKE_dsysvx)dlsym(dsysvx_obj->hModule, "LAPACKE_dsysvx");
    ASSERT_TRUE(DSYSVX != NULL) << "failed to syt the Netlib LAPACKE_dsysvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the sytrf API to compute the factorized A.  */
    if(dsysvx_obj->fact == 'F') {
        DSYTRF = (Fptr_NL_LAPACKE_dsytrf)dlsym(dsysvx_obj->hModule,"LAPACKE_dsytrf");
        ASSERT_TRUE(DSYTRF != NULL) << "failed to syt the Netlib LAPACKE_dsytrf symbol";
            
        dsysvx_obj->inforef = DSYTRF( dsysvx_obj->matrix_layout,
                                      dsysvx_obj->uplo, dsysvx_obj->n,
                                      dsysvx_obj->afref,
                                      dsysvx_obj->lda,
                                      dsysvx_obj->ipivref);
                               
        dsysvx_obj->info = LAPACKE_dsytrf( dsysvx_obj->matrix_layout,
                                           dsysvx_obj->uplo, dsysvx_obj->n,
                                           dsysvx_obj->af,
                                           dsysvx_obj->lda,
                                           dsysvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dsysvx_obj->inforef = DSYSVX( dsysvx_obj->matrix_layout, dsysvx_obj->fact,
                                  dsysvx_obj->uplo, dsysvx_obj->n,
                                  dsysvx_obj->nrhs,
                                  dsysvx_obj->aref, dsysvx_obj->lda, 
                                  dsysvx_obj->afref, dsysvx_obj->ldaf,
                                  dsysvx_obj->ipivref,
                                  dsysvx_obj->bref, dsysvx_obj->ldb,
                                  dsysvx_obj->xref, dsysvx_obj->ldx,
                                  &dsysvx_obj->rcondref, 
                                  dsysvx_obj->ferrref,
                                  dsysvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dsysvx_obj->info = LAPACKE_dsysvx( dsysvx_obj->matrix_layout, dsysvx_obj->fact,
                                  dsysvx_obj->uplo, dsysvx_obj->n,
                                  dsysvx_obj->nrhs,
                                  dsysvx_obj->a, dsysvx_obj->lda, 
                                  dsysvx_obj->af, dsysvx_obj->ldaf,
                                  dsysvx_obj->ipiv,
                                  dsysvx_obj->b, dsysvx_obj->ldb,
                                  dsysvx_obj->x, dsysvx_obj->ldx,
                                  &dsysvx_obj->rcond, 
                                  dsysvx_obj->ferr,
                                  dsysvx_obj->berr);

    if( dsysvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsysvx is wrong\n", dsysvx_obj->info );
    }
    if( dsysvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsysvx is wrong\n", 
        dsysvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dsysvx_obj->ipiv_diff = computeDiff_i( dsysvx_obj->n, dsysvx_obj->ipiv, dsysvx_obj->ipivref );
    
    dsysvx_obj->diff =  computeDiff_d( dsysvx_obj->n * dsysvx_obj->n, 
                dsysvx_obj->af, dsysvx_obj->afref );

    dsysvx_obj->diff_xerr =  computeDiff_d( dsysvx_obj->x_bufsize, 
                dsysvx_obj->x, dsysvx_obj->xref );

    dsysvx_obj->diff_berr =  computeDiff_d( dsysvx_obj->nrhs, 
                dsysvx_obj->berr, dsysvx_obj->berrref );
                
    dsysvx_obj->diff_ferr =  computeDiff_d( dsysvx_obj->nrhs, 
                dsysvx_obj->ferr, dsysvx_obj->ferrref );
}

TEST_F(dsysvx_test, dsysvx1) {
    EXPECT_NEAR(0.0, dsysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dsysvx_obj->rcond - dsysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, dsysvx_obj->ipiv_diff );
}

TEST_F(dsysvx_test, dsysvx2) {
    EXPECT_NEAR(0.0, dsysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dsysvx_obj->rcond - dsysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsysvx_test, dsysvx3) {
    EXPECT_NEAR(0.0, dsysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dsysvx_obj->rcond - dsysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsysvx_test, dsysvx4) {
    EXPECT_NEAR(0.0, dsysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dsysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dsysvx_obj->rcond - dsysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin sysvx_scomplex_parameters  class definition */
class sysvx_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
      int ipiv_diff;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_float *a, *aref; //The array ab contains the matrix A
      lapack_complex_float *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sysvx_scomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~sysvx_scomplex_parameters (); 
};  /* end of sysvx_scomplex_parameters  class definition */


/* Constructor sysvx_scomplex_parameters definition */
sysvx_scomplex_parameters:: sysvx_scomplex_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n sysvx lapack_complex_float:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);

    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       sysvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(lapack_complex_float)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

sysvx_scomplex_parameters:: ~sysvx_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysvx_scomplex_parameters object: destructor invoked. \n");
#endif
   sysvx_free();
}

//  Test fixture class definition
class csysvx_test  : public  ::testing::Test {
public:
   sysvx_scomplex_parameters  *csysvx_obj;
   void SetUp();  
   void TearDown () { delete csysvx_obj; }
};


void csysvx_test::SetUp(){

    /* LAPACKE CSYSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_csysvx) (int matrix_layout, char fact,
                char uplo, lapack_int n, lapack_int nrhs, 
                const lapack_complex_float* a, lapack_int lda, 
                lapack_complex_float* af, lapack_int ldaf, lapack_int* ipiv,
                const lapack_complex_float* b, lapack_int ldb,
                lapack_complex_float* x, lapack_int ldx,
                float* rcond, float* ferr, float* berr);

    Fptr_NL_LAPACKE_csysvx CSYSVX;

     /* LAPACKE CSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf) ( int matrix_layout , char uplo,
            lapack_int n , lapack_complex_float * a , lapack_int lda , lapack_int * ipiv );

    Fptr_NL_LAPACKE_csytrf CSYTRF;

    csysvx_obj = new sysvx_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    csysvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csysvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csysvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csysvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CSYSVX = (Fptr_NL_LAPACKE_csysvx)dlsym(csysvx_obj->hModule, "LAPACKE_csysvx");
    ASSERT_TRUE(CSYSVX != NULL) << "failed to syt the Netlib LAPACKE_csysvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the sytrf API to compute the factorized A.  */
    if(csysvx_obj->fact == 'F') {
        CSYTRF = (Fptr_NL_LAPACKE_csytrf)dlsym(csysvx_obj->hModule,"LAPACKE_csytrf");
        ASSERT_TRUE(CSYTRF != NULL) << "failed to syt the Netlib LAPACKE_csytrf symbol";
            
        csysvx_obj->inforef = CSYTRF( csysvx_obj->matrix_layout,
                                      csysvx_obj->uplo, csysvx_obj->n,
                                      csysvx_obj->afref,
                                      csysvx_obj->lda,
                                      csysvx_obj->ipivref);
                               
        csysvx_obj->info = LAPACKE_csytrf( csysvx_obj->matrix_layout,
                                           csysvx_obj->uplo, csysvx_obj->n,
                                           csysvx_obj->af,
                                           csysvx_obj->lda,
                                           csysvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    csysvx_obj->inforef = CSYSVX( csysvx_obj->matrix_layout, csysvx_obj->fact,
                                  csysvx_obj->uplo, csysvx_obj->n,
                                  csysvx_obj->nrhs,
                                  csysvx_obj->aref, csysvx_obj->lda, 
                                  csysvx_obj->afref, csysvx_obj->ldaf,
                                  csysvx_obj->ipivref,
                                  csysvx_obj->bref, csysvx_obj->ldb,
                                  csysvx_obj->xref, csysvx_obj->ldx,
                                  &csysvx_obj->rcondref, 
                                  csysvx_obj->ferrref,
                                  csysvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    csysvx_obj->info = LAPACKE_csysvx( csysvx_obj->matrix_layout, csysvx_obj->fact,
                                  csysvx_obj->uplo, csysvx_obj->n,
                                  csysvx_obj->nrhs,
                                  csysvx_obj->a, csysvx_obj->lda, 
                                  csysvx_obj->af, csysvx_obj->ldaf,
                                  csysvx_obj->ipiv,
                                  csysvx_obj->b, csysvx_obj->ldb,
                                  csysvx_obj->x, csysvx_obj->ldx,
                                  &csysvx_obj->rcond, 
                                  csysvx_obj->ferr,
                                  csysvx_obj->berr);

    if( csysvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csysvx is wrong\n", csysvx_obj->info );
    }
    if( csysvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csysvx is wrong\n", 
        csysvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    csysvx_obj->ipiv_diff = computeDiff_i( csysvx_obj->n, csysvx_obj->ipiv, csysvx_obj->ipivref );
    
    csysvx_obj->diff =  computeDiff_c( csysvx_obj->n * csysvx_obj->n, 
                csysvx_obj->af, csysvx_obj->afref );

    csysvx_obj->diff_xerr =  computeDiff_c( csysvx_obj->x_bufsize, 
                csysvx_obj->x, csysvx_obj->xref );

    csysvx_obj->diff_berr =  computeDiff_s( csysvx_obj->nrhs, 
                csysvx_obj->berr, csysvx_obj->berrref );
                
    csysvx_obj->diff_ferr =  computeDiff_s( csysvx_obj->nrhs, 
                csysvx_obj->ferr, csysvx_obj->ferrref );
}

TEST_F(csysvx_test, csysvx1) {
    EXPECT_NEAR(0.0, csysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (csysvx_obj->rcond - csysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, csysvx_obj->ipiv_diff );
}

TEST_F(csysvx_test, csysvx2) {
    EXPECT_NEAR(0.0, csysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (csysvx_obj->rcond - csysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csysvx_test, csysvx3) {
    EXPECT_NEAR(0.0, csysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (csysvx_obj->rcond - csysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csysvx_test, csysvx4) {
    EXPECT_NEAR(0.0, csysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, csysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (csysvx_obj->rcond - csysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin sysvx_dcomplex_parameters  class definition */
class sysvx_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
      int ipiv_diff;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_double *a, *aref; //The array ab contains the matrix A
      lapack_complex_double *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.
      lapack_int *ipiv, *ipivref; // pivot buffer

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sysvx_dcomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~sysvx_dcomplex_parameters (); 
};  /* end of sysvx_dcomplex_parameters  class definition */


/* Constructor sysvx_dcomplex_parameters definition */
sysvx_dcomplex_parameters:: sysvx_dcomplex_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    n = n_i;
    nrhs = nrhs_i;
    lda = n;
    ldaf = n;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n sysvx lapack_complex_double:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);

    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       sysvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(lapack_complex_double)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

sysvx_dcomplex_parameters:: ~sysvx_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysvx_dcomplex_parameters object: destructor invoked. \n");
#endif
   sysvx_free();
}

//  Test fixture class definition
class zsysvx_test  : public  ::testing::Test {
public:
   sysvx_dcomplex_parameters  *zsysvx_obj;
   void SetUp();  
   void TearDown () { delete zsysvx_obj; }
};


void zsysvx_test::SetUp(){

    /* LAPACKE ZSYSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zsysvx) (int matrix_layout, char fact,
                char uplo, lapack_int n, lapack_int nrhs, 
                const lapack_complex_double* a, lapack_int lda, 
                lapack_complex_double* af, lapack_int ldaf, lapack_int* ipiv,
                const lapack_complex_double* b, lapack_int ldb,
                lapack_complex_double* x, lapack_int ldx,
                double* rcond, double* ferr, double* berr);

    Fptr_NL_LAPACKE_zsysvx ZSYSVX;

     /* LAPACKE ZSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf) ( int matrix_layout , char uplo,
            lapack_int n , lapack_complex_double * a , lapack_int lda , lapack_int * ipiv );

    Fptr_NL_LAPACKE_zsytrf ZSYTRF;

    zsysvx_obj = new sysvx_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zsysvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsysvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsysvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsysvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZSYSVX = (Fptr_NL_LAPACKE_zsysvx)dlsym(zsysvx_obj->hModule, "LAPACKE_zsysvx");
    ASSERT_TRUE(ZSYSVX != NULL) << "failed to syt the Netlib LAPACKE_zsysvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the sytrf API to compute the factorized A.  */
    if(zsysvx_obj->fact == 'F') {
        ZSYTRF = (Fptr_NL_LAPACKE_zsytrf)dlsym(zsysvx_obj->hModule,"LAPACKE_zsytrf");
        ASSERT_TRUE(ZSYTRF != NULL) << "failed to syt the Netlib LAPACKE_zsytrf symbol";
            
        zsysvx_obj->inforef = ZSYTRF( zsysvx_obj->matrix_layout,
                                      zsysvx_obj->uplo, zsysvx_obj->n,
                                      zsysvx_obj->afref,
                                      zsysvx_obj->lda,
                                      zsysvx_obj->ipivref);
                               
        zsysvx_obj->info = LAPACKE_zsytrf( zsysvx_obj->matrix_layout,
                                           zsysvx_obj->uplo, zsysvx_obj->n,
                                           zsysvx_obj->af,
                                           zsysvx_obj->lda,
                                           zsysvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zsysvx_obj->inforef = ZSYSVX( zsysvx_obj->matrix_layout, zsysvx_obj->fact,
                                  zsysvx_obj->uplo, zsysvx_obj->n,
                                  zsysvx_obj->nrhs,
                                  zsysvx_obj->aref, zsysvx_obj->lda, 
                                  zsysvx_obj->afref, zsysvx_obj->ldaf,
                                  zsysvx_obj->ipivref,
                                  zsysvx_obj->bref, zsysvx_obj->ldb,
                                  zsysvx_obj->xref, zsysvx_obj->ldx,
                                  &zsysvx_obj->rcondref, 
                                  zsysvx_obj->ferrref,
                                  zsysvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zsysvx_obj->info = LAPACKE_zsysvx( zsysvx_obj->matrix_layout, zsysvx_obj->fact,
                                  zsysvx_obj->uplo, zsysvx_obj->n,
                                  zsysvx_obj->nrhs,
                                  zsysvx_obj->a, zsysvx_obj->lda, 
                                  zsysvx_obj->af, zsysvx_obj->ldaf,
                                  zsysvx_obj->ipiv,
                                  zsysvx_obj->b, zsysvx_obj->ldb,
                                  zsysvx_obj->x, zsysvx_obj->ldx,
                                  &zsysvx_obj->rcond, 
                                  zsysvx_obj->ferr,
                                  zsysvx_obj->berr);

    if( zsysvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsysvx is wrong\n", zsysvx_obj->info );
    }
    if( zsysvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsysvx is wrong\n", 
        zsysvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zsysvx_obj->ipiv_diff = computeDiff_i( zsysvx_obj->n, zsysvx_obj->ipiv, zsysvx_obj->ipivref );
    
    zsysvx_obj->diff =  computeDiff_z( zsysvx_obj->n * zsysvx_obj->n, 
                zsysvx_obj->af, zsysvx_obj->afref );

    zsysvx_obj->diff_xerr =  computeDiff_z( zsysvx_obj->x_bufsize, 
                zsysvx_obj->x, zsysvx_obj->xref );

    zsysvx_obj->diff_berr =  computeDiff_d( zsysvx_obj->nrhs, 
                zsysvx_obj->berr, zsysvx_obj->berrref );
                
    zsysvx_obj->diff_ferr =  computeDiff_d( zsysvx_obj->nrhs, 
                zsysvx_obj->ferr, zsysvx_obj->ferrref );
}

TEST_F(zsysvx_test, zsysvx1) {
    EXPECT_NEAR(0.0, zsysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zsysvx_obj->rcond - zsysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, zsysvx_obj->ipiv_diff );
}

TEST_F(zsysvx_test, zsysvx2) {
    EXPECT_NEAR(0.0, zsysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zsysvx_obj->rcond - zsysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsysvx_test, zsysvx3) {
    EXPECT_NEAR(0.0, zsysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zsysvx_obj->rcond - zsysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsysvx_test, zsysvx4) {
    EXPECT_NEAR(0.0, zsysvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zsysvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zsysvx_obj->rcond - zsysvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}
