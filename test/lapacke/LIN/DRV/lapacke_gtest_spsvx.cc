#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define spsvx_free() \
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

/* Begin spsvx_float_parameters  class definition */
class spsvx_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      int a_bufsize;
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
      spsvx_float_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~spsvx_float_parameters (); 
};  /* end of spsvx_float_parameters  class definition */


/* Constructor spsvx_float_parameters definition */
spsvx_float_parameters:: spsvx_float_parameters ( int matrix_layout_i, 
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
    
    a_bufsize = n*(n+1)/2;
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
   printf(" \n spsvx float:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &af, &afref, a_bufsize);
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
       spsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand ( a, aref, a_bufsize);
    memcpy(af, a, (a_bufsize*sizeof(float)));
    memcpy(afref, a, (a_bufsize*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_float_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

spsvx_float_parameters:: ~spsvx_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" spsvx_float_parameters object: destructor invoked. \n");
#endif
   spsvx_free();
}

//  Test fixture class definition
class sspsvx_test  : public  ::testing::Test {
public:
   spsvx_float_parameters  *sspsvx_obj;
   void SetUp();  
   void TearDown () { delete sspsvx_obj; }
};


void sspsvx_test::SetUp(){

    /* LAPACKE SSPSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_sspsvx) (int matrix_layout, char fact,
            char uplo, lapack_int n, lapack_int nrhs, const float* ap,
            float* afp, lapack_int* ipiv, const float* b, lapack_int ldb,
            float* x, lapack_int ldx, float* rcond, float* ferr, float* berr);

    Fptr_NL_LAPACKE_sspsvx SSPSVX;

     /* LAPACKE SSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_ssptrf) ( int matrix_layout , char uplo,
                            lapack_int n , float * ap , lapack_int * ipiv );

    Fptr_NL_LAPACKE_ssptrf SSPTRF;

    sspsvx_obj = new spsvx_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    sspsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sspsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sspsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sspsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SSPSVX = (Fptr_NL_LAPACKE_sspsvx)dlsym(sspsvx_obj->hModule, "LAPACKE_sspsvx");
    ASSERT_TRUE(SSPSVX != NULL) << "failed to syt the Netlib LAPACKE_sspsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the sptrf API to compute the factorized A.  */
    if(sspsvx_obj->fact == 'F') {
        SSPTRF = (Fptr_NL_LAPACKE_ssptrf)dlsym(sspsvx_obj->hModule,"LAPACKE_ssptrf");
        ASSERT_TRUE(SSPTRF != NULL) << "failed to syt the Netlib LAPACKE_ssptrf symbol";
        sspsvx_obj->inforef = SSPTRF( sspsvx_obj->matrix_layout,
                                      sspsvx_obj->uplo, sspsvx_obj->n,
                                      sspsvx_obj->afref,
                                      sspsvx_obj->ipivref);
                               
        sspsvx_obj->info = LAPACKE_ssptrf( sspsvx_obj->matrix_layout,
                                           sspsvx_obj->uplo, sspsvx_obj->n,
                                           sspsvx_obj->af,
                                           sspsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    sspsvx_obj->inforef = SSPSVX( sspsvx_obj->matrix_layout, sspsvx_obj->fact,
                                  sspsvx_obj->uplo, sspsvx_obj->n,
                                  sspsvx_obj->nrhs,
                                  sspsvx_obj->aref,
                                  sspsvx_obj->afref,
                                  sspsvx_obj->ipivref,
                                  sspsvx_obj->bref, sspsvx_obj->ldb,
                                  sspsvx_obj->xref, sspsvx_obj->ldx,
                                  &sspsvx_obj->rcondref, 
                                  sspsvx_obj->ferrref,
                                  sspsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    sspsvx_obj->info = LAPACKE_sspsvx( sspsvx_obj->matrix_layout, sspsvx_obj->fact,
                                  sspsvx_obj->uplo, sspsvx_obj->n,
                                  sspsvx_obj->nrhs,
                                  sspsvx_obj->a,
                                  sspsvx_obj->af,
                                  sspsvx_obj->ipiv,
                                  sspsvx_obj->b, sspsvx_obj->ldb,
                                  sspsvx_obj->x, sspsvx_obj->ldx,
                                  &sspsvx_obj->rcond, 
                                  sspsvx_obj->ferr,
                                  sspsvx_obj->berr);

    if( sspsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sspsvx is wrong\n", sspsvx_obj->info );
    }
    if( sspsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sspsvx is wrong\n", 
        sspsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sspsvx_obj->ipiv_diff = computeDiff_i( sspsvx_obj->n, sspsvx_obj->ipiv, sspsvx_obj->ipivref );
    
    sspsvx_obj->diff =  computeDiff_s( sspsvx_obj->a_bufsize, 
                sspsvx_obj->af, sspsvx_obj->afref );

    sspsvx_obj->diff_xerr =  computeDiff_s( sspsvx_obj->x_bufsize, 
                sspsvx_obj->x, sspsvx_obj->xref );

    sspsvx_obj->diff_berr =  computeDiff_s( sspsvx_obj->nrhs, 
                sspsvx_obj->berr, sspsvx_obj->berrref );
                
    sspsvx_obj->diff_ferr =  computeDiff_s( sspsvx_obj->nrhs, 
                sspsvx_obj->ferr, sspsvx_obj->ferrref );
}

TEST_F(sspsvx_test, sspsvx1) {
    EXPECT_NEAR(0.0, sspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sspsvx_obj->rcond - sspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, sspsvx_obj->ipiv_diff );
}

TEST_F(sspsvx_test, sspsvx2) {
    EXPECT_NEAR(0.0, sspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sspsvx_obj->rcond - sspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sspsvx_test, sspsvx3) {
    EXPECT_NEAR(0.0, sspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sspsvx_obj->rcond - sspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sspsvx_test, sspsvx4) {
    EXPECT_NEAR(0.0, sspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sspsvx_obj->rcond - sspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin spsvx_double_parameters  class definition */
class spsvx_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      int a_bufsize;
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
      spsvx_double_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~spsvx_double_parameters (); 
};  /* end of spsvx_double_parameters  class definition */


/* Constructor spsvx_double_parameters definition */
spsvx_double_parameters:: spsvx_double_parameters ( int matrix_layout_i, 
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
    
    a_bufsize = n*(n+1)/2;
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
   printf(" \n spsvx double:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &af, &afref, a_bufsize);
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
       spsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand ( a, aref, a_bufsize);
    memcpy(af, a, (a_bufsize*sizeof(double)));
    memcpy(afref, a, (a_bufsize*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_double_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

spsvx_double_parameters:: ~spsvx_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" spsvx_double_parameters object: destructor invoked. \n");
#endif
   spsvx_free();
}

//  Test fixture class definition
class dspsvx_test  : public  ::testing::Test {
public:
   spsvx_double_parameters  *dspsvx_obj;
   void SetUp();  
   void TearDown () { delete dspsvx_obj; }
};


void dspsvx_test::SetUp(){

    /* LAPACKE DSPSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_dspsvx) (int matrix_layout, char fact,
            char uplo, lapack_int n, lapack_int nrhs, const double* ap,
            double* afp, lapack_int* ipiv, const double* b, lapack_int ldb,
            double* x, lapack_int ldx, double* rcond, double* ferr, double* berr);

    Fptr_NL_LAPACKE_dspsvx DSPSVX;

     /* LAPACKE DSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dsptrf) ( int matrix_layout , char uplo,
                            lapack_int n , double * ap , lapack_int * ipiv );

    Fptr_NL_LAPACKE_dsptrf DSPTRF;

    dspsvx_obj = new spsvx_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dspsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dspsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dspsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dspsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DSPSVX = (Fptr_NL_LAPACKE_dspsvx)dlsym(dspsvx_obj->hModule, "LAPACKE_dspsvx");
    ASSERT_TRUE(DSPSVX != NULL) << "failed to syt the Netlib LAPACKE_dspsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the sptrf API to compute the factorized A.  */
    if(dspsvx_obj->fact == 'F') {
        DSPTRF = (Fptr_NL_LAPACKE_dsptrf)dlsym(dspsvx_obj->hModule,"LAPACKE_dsptrf");
        ASSERT_TRUE(DSPTRF != NULL) << "failed to syt the Netlib LAPACKE_dsptrf symbol";
        dspsvx_obj->inforef = DSPTRF( dspsvx_obj->matrix_layout,
                                      dspsvx_obj->uplo, dspsvx_obj->n,
                                      dspsvx_obj->afref,
                                      dspsvx_obj->ipivref);
                               
        dspsvx_obj->info = LAPACKE_dsptrf( dspsvx_obj->matrix_layout,
                                           dspsvx_obj->uplo, dspsvx_obj->n,
                                           dspsvx_obj->af,
                                           dspsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dspsvx_obj->inforef = DSPSVX( dspsvx_obj->matrix_layout, dspsvx_obj->fact,
                                  dspsvx_obj->uplo, dspsvx_obj->n,
                                  dspsvx_obj->nrhs,
                                  dspsvx_obj->aref,
                                  dspsvx_obj->afref,
                                  dspsvx_obj->ipivref,
                                  dspsvx_obj->bref, dspsvx_obj->ldb,
                                  dspsvx_obj->xref, dspsvx_obj->ldx,
                                  &dspsvx_obj->rcondref, 
                                  dspsvx_obj->ferrref,
                                  dspsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dspsvx_obj->info = LAPACKE_dspsvx( dspsvx_obj->matrix_layout, dspsvx_obj->fact,
                                  dspsvx_obj->uplo, dspsvx_obj->n,
                                  dspsvx_obj->nrhs,
                                  dspsvx_obj->a,
                                  dspsvx_obj->af,
                                  dspsvx_obj->ipiv,
                                  dspsvx_obj->b, dspsvx_obj->ldb,
                                  dspsvx_obj->x, dspsvx_obj->ldx,
                                  &dspsvx_obj->rcond, 
                                  dspsvx_obj->ferr,
                                  dspsvx_obj->berr);

    if( dspsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dspsvx is wrong\n", dspsvx_obj->info );
    }
    if( dspsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dspsvx is wrong\n", 
        dspsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dspsvx_obj->ipiv_diff = computeDiff_i( dspsvx_obj->n, dspsvx_obj->ipiv, dspsvx_obj->ipivref );
    
    dspsvx_obj->diff =  computeDiff_d( dspsvx_obj->a_bufsize, 
                dspsvx_obj->af, dspsvx_obj->afref );

    dspsvx_obj->diff_xerr =  computeDiff_d( dspsvx_obj->x_bufsize, 
                dspsvx_obj->x, dspsvx_obj->xref );

    dspsvx_obj->diff_berr =  computeDiff_d( dspsvx_obj->nrhs, 
                dspsvx_obj->berr, dspsvx_obj->berrref );
                
    dspsvx_obj->diff_ferr =  computeDiff_d( dspsvx_obj->nrhs, 
                dspsvx_obj->ferr, dspsvx_obj->ferrref );
}

TEST_F(dspsvx_test, dspsvx1) {
    EXPECT_NEAR(0.0, dspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dspsvx_obj->rcond - dspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, dspsvx_obj->ipiv_diff );
}

TEST_F(dspsvx_test, dspsvx2) {
    EXPECT_NEAR(0.0, dspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dspsvx_obj->rcond - dspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dspsvx_test, dspsvx3) {
    EXPECT_NEAR(0.0, dspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dspsvx_obj->rcond - dspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dspsvx_test, dspsvx4) {
    EXPECT_NEAR(0.0, dspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dspsvx_obj->rcond - dspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin spsvx_scomplex_parameters  class definition */
class spsvx_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      int a_bufsize;
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
      spsvx_scomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~spsvx_scomplex_parameters (); 
};  /* end of spsvx_scomplex_parameters  class definition */


/* Constructor spsvx_scomplex_parameters definition */
spsvx_scomplex_parameters:: spsvx_scomplex_parameters ( int matrix_layout_i, 
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
    
    a_bufsize = n*(n+1)/2;
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
   printf(" \n spsvx lapack_complex_float:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &af, &afref, a_bufsize);
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
       spsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand ( a, aref, a_bufsize);
    memcpy(af, a, (a_bufsize*sizeof(lapack_complex_float)));
    memcpy(afref, a, (a_bufsize*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

spsvx_scomplex_parameters:: ~spsvx_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" spsvx_scomplex_parameters object: destructor invoked. \n");
#endif
   spsvx_free();
}

//  Test fixture class definition
class cspsvx_test  : public  ::testing::Test {
public:
   spsvx_scomplex_parameters  *cspsvx_obj;
   void SetUp();  
   void TearDown () { delete cspsvx_obj; }
};


void cspsvx_test::SetUp(){

    /* LAPACKE CSPSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_cspsvx) (int matrix_layout, char fact,
                               char uplo, lapack_int n, lapack_int nrhs,
              const lapack_complex_float* ap, lapack_complex_float* afp,
        lapack_int* ipiv, const lapack_complex_float* b, lapack_int ldb,
                  lapack_complex_float* x, lapack_int ldx, float* rcond,
                                            float* ferr, float* berr);

    Fptr_NL_LAPACKE_cspsvx CSPSVX;

     /* LAPACKE CSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_csptrf) ( int matrix_layout , char uplo,
            lapack_int n , lapack_complex_float * a , lapack_int * ipiv );

    Fptr_NL_LAPACKE_csptrf CSPTRF;

    cspsvx_obj = new spsvx_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    cspsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cspsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cspsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cspsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CSPSVX = (Fptr_NL_LAPACKE_cspsvx)dlsym(cspsvx_obj->hModule, "LAPACKE_cspsvx");
    ASSERT_TRUE(CSPSVX != NULL) << "failed to syt the Netlib LAPACKE_cspsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the sptrf API to compute the factorized A.  */
    if(cspsvx_obj->fact == 'F') {
        CSPTRF = (Fptr_NL_LAPACKE_csptrf)dlsym(cspsvx_obj->hModule,"LAPACKE_csptrf");
        ASSERT_TRUE(CSPTRF != NULL) << "failed to syt the Netlib LAPACKE_csptrf symbol";
            
        cspsvx_obj->inforef = CSPTRF( cspsvx_obj->matrix_layout,
                                      cspsvx_obj->uplo, cspsvx_obj->n,
                                      cspsvx_obj->afref,
                                      cspsvx_obj->ipivref);
                               
        cspsvx_obj->info = LAPACKE_csptrf( cspsvx_obj->matrix_layout,
                                           cspsvx_obj->uplo, cspsvx_obj->n,
                                           cspsvx_obj->af,
                                           cspsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cspsvx_obj->inforef = CSPSVX( cspsvx_obj->matrix_layout, cspsvx_obj->fact,
                                  cspsvx_obj->uplo, cspsvx_obj->n,
                                  cspsvx_obj->nrhs,
                                  cspsvx_obj->aref,
                                  cspsvx_obj->afref,
                                  cspsvx_obj->ipivref,
                                  cspsvx_obj->bref, cspsvx_obj->ldb,
                                  cspsvx_obj->xref, cspsvx_obj->ldx,
                                  &cspsvx_obj->rcondref, 
                                  cspsvx_obj->ferrref,
                                  cspsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cspsvx_obj->info = LAPACKE_cspsvx( cspsvx_obj->matrix_layout, cspsvx_obj->fact,
                                  cspsvx_obj->uplo, cspsvx_obj->n,
                                  cspsvx_obj->nrhs,
                                  cspsvx_obj->a,
                                  cspsvx_obj->af,
                                  cspsvx_obj->ipiv,
                                  cspsvx_obj->b, cspsvx_obj->ldb,
                                  cspsvx_obj->x, cspsvx_obj->ldx,
                                  &cspsvx_obj->rcond, 
                                  cspsvx_obj->ferr,
                                  cspsvx_obj->berr);

    if( cspsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cspsvx is wrong\n", cspsvx_obj->info );
    }
    if( cspsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cspsvx is wrong\n", 
        cspsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cspsvx_obj->ipiv_diff = computeDiff_i( cspsvx_obj->n, cspsvx_obj->ipiv, cspsvx_obj->ipivref );
    
    cspsvx_obj->diff =  computeDiff_c( cspsvx_obj->a_bufsize, 
                cspsvx_obj->af, cspsvx_obj->afref );

    cspsvx_obj->diff_xerr =  computeDiff_c( cspsvx_obj->x_bufsize, 
                cspsvx_obj->x, cspsvx_obj->xref );

    cspsvx_obj->diff_berr =  computeDiff_s( cspsvx_obj->nrhs, 
                cspsvx_obj->berr, cspsvx_obj->berrref );
                
    cspsvx_obj->diff_ferr =  computeDiff_s( cspsvx_obj->nrhs, 
                cspsvx_obj->ferr, cspsvx_obj->ferrref );
}

TEST_F(cspsvx_test, cspsvx1) {
    EXPECT_NEAR(0.0, cspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cspsvx_obj->rcond - cspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, cspsvx_obj->ipiv_diff );
}

TEST_F(cspsvx_test, cspsvx2) {
    EXPECT_NEAR(0.0, cspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cspsvx_obj->rcond - cspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cspsvx_test, cspsvx3) {
    EXPECT_NEAR(0.0, cspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cspsvx_obj->rcond - cspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cspsvx_test, cspsvx4) {
    EXPECT_NEAR(0.0, cspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cspsvx_obj->rcond - cspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin spsvx_dcomplex_parameters  class definition */
class spsvx_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      int a_bufsize;
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
      spsvx_dcomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~spsvx_dcomplex_parameters (); 
};  /* end of spsvx_dcomplex_parameters  class definition */


/* Constructor spsvx_dcomplex_parameters definition */
spsvx_dcomplex_parameters:: spsvx_dcomplex_parameters ( int matrix_layout_i, 
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
    
    a_bufsize = n*(n+1)/2;
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
   printf(" \n spsvx lapack_complex_double:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &af, &afref, a_bufsize);
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
       spsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand ( a, aref, a_bufsize);
    memcpy(af, a, (a_bufsize*sizeof(lapack_complex_double)));
    memcpy(afref, a, (a_bufsize*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

spsvx_dcomplex_parameters:: ~spsvx_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" spsvx_dcomplex_parameters object: destructor invoked. \n");
#endif
   spsvx_free();
}

//  Test fixture class definition
class zspsvx_test  : public  ::testing::Test {
public:
   spsvx_dcomplex_parameters  *zspsvx_obj;
   void SetUp();  
   void TearDown () { delete zspsvx_obj; }
};


void zspsvx_test::SetUp(){

    /* LAPACKE ZSPSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zspsvx) (int matrix_layout, char fact,
                               char uplo, lapack_int n, lapack_int nrhs,
              const lapack_complex_double* ap, lapack_complex_double* afp,
        lapack_int* ipiv, const lapack_complex_double* b, lapack_int ldb,
                  lapack_complex_double* x, lapack_int ldx, double* rcond,
                                            double* ferr, double* berr);

    Fptr_NL_LAPACKE_zspsvx ZSPSVX;

     /* LAPACKE ZSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zsptrf) ( int matrix_layout , char uplo,
            lapack_int n , lapack_complex_double * a , lapack_int * ipiv );

    Fptr_NL_LAPACKE_zsptrf ZSPTRF;

    zspsvx_obj = new spsvx_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zspsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zspsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zspsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zspsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZSPSVX = (Fptr_NL_LAPACKE_zspsvx)dlsym(zspsvx_obj->hModule, "LAPACKE_zspsvx");
    ASSERT_TRUE(ZSPSVX != NULL) << "failed to syt the Netlib LAPACKE_zspsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the sptrf API to compute the factorized A.  */
    if(zspsvx_obj->fact == 'F') {
        ZSPTRF = (Fptr_NL_LAPACKE_zsptrf)dlsym(zspsvx_obj->hModule,"LAPACKE_zsptrf");
        ASSERT_TRUE(ZSPTRF != NULL) << "failed to syt the Netlib LAPACKE_zsptrf symbol";
            
        zspsvx_obj->inforef = ZSPTRF( zspsvx_obj->matrix_layout,
                                      zspsvx_obj->uplo, zspsvx_obj->n,
                                      zspsvx_obj->afref,
                                      zspsvx_obj->ipivref);
                               
        zspsvx_obj->info = LAPACKE_zsptrf( zspsvx_obj->matrix_layout,
                                           zspsvx_obj->uplo, zspsvx_obj->n,
                                           zspsvx_obj->af,
                                           zspsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zspsvx_obj->inforef = ZSPSVX( zspsvx_obj->matrix_layout, zspsvx_obj->fact,
                                  zspsvx_obj->uplo, zspsvx_obj->n,
                                  zspsvx_obj->nrhs,
                                  zspsvx_obj->aref,
                                  zspsvx_obj->afref,
                                  zspsvx_obj->ipivref,
                                  zspsvx_obj->bref, zspsvx_obj->ldb,
                                  zspsvx_obj->xref, zspsvx_obj->ldx,
                                  &zspsvx_obj->rcondref, 
                                  zspsvx_obj->ferrref,
                                  zspsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zspsvx_obj->info = LAPACKE_zspsvx( zspsvx_obj->matrix_layout, zspsvx_obj->fact,
                                  zspsvx_obj->uplo, zspsvx_obj->n,
                                  zspsvx_obj->nrhs,
                                  zspsvx_obj->a,
                                  zspsvx_obj->af,
                                  zspsvx_obj->ipiv,
                                  zspsvx_obj->b, zspsvx_obj->ldb,
                                  zspsvx_obj->x, zspsvx_obj->ldx,
                                  &zspsvx_obj->rcond, 
                                  zspsvx_obj->ferr,
                                  zspsvx_obj->berr);

    if( zspsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zspsvx is wrong\n", zspsvx_obj->info );
    }
    if( zspsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zspsvx is wrong\n", 
        zspsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zspsvx_obj->ipiv_diff = computeDiff_i( zspsvx_obj->n, zspsvx_obj->ipiv, zspsvx_obj->ipivref );
    
    zspsvx_obj->diff =  computeDiff_z( zspsvx_obj->a_bufsize, 
                zspsvx_obj->af, zspsvx_obj->afref );

    zspsvx_obj->diff_xerr =  computeDiff_z( zspsvx_obj->x_bufsize, 
                zspsvx_obj->x, zspsvx_obj->xref );

    zspsvx_obj->diff_berr =  computeDiff_d( zspsvx_obj->nrhs, 
                zspsvx_obj->berr, zspsvx_obj->berrref );
                
    zspsvx_obj->diff_ferr =  computeDiff_d( zspsvx_obj->nrhs, 
                zspsvx_obj->ferr, zspsvx_obj->ferrref );
}

TEST_F(zspsvx_test, zspsvx1) {
    EXPECT_NEAR(0.0, zspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zspsvx_obj->rcond - zspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, zspsvx_obj->ipiv_diff );
}

TEST_F(zspsvx_test, zspsvx2) {
    EXPECT_NEAR(0.0, zspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zspsvx_obj->rcond - zspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zspsvx_test, zspsvx3) {
    EXPECT_NEAR(0.0, zspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zspsvx_obj->rcond - zspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zspsvx_test, zspsvx4) {
    EXPECT_NEAR(0.0, zspsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zspsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zspsvx_obj->rcond - zspsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}