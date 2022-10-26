#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define hpsvx_free() \
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

/* Begin hpsvx_scomplex_parameters  class definition */
class hpsvx_scomplex_parameters{
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
      hpsvx_scomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~hpsvx_scomplex_parameters (); 
};  /* end of hpsvx_scomplex_parameters  class definition */


/* Constructor hpsvx_scomplex_parameters definition */
hpsvx_scomplex_parameters:: hpsvx_scomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n hpsvx lapack_complex_float:  n: %d, fact: %c uplo: %c  lda: %d  \
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
       hpsvx_free();
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

hpsvx_scomplex_parameters:: ~hpsvx_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hpsvx_scomplex_parameters object: destructor invoked. \n");
#endif
   hpsvx_free();
}

//  Test fixture class definition
class chpsvx_test  : public  ::testing::Test {
public:
   hpsvx_scomplex_parameters  *chpsvx_obj;
   void SetUp();  
   void TearDown () { delete chpsvx_obj; }
};


void chpsvx_test::SetUp(){

    /* LAPACKE CHPSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_chpsvx) (int matrix_layout, char fact,
                               char uplo, lapack_int n, lapack_int nrhs,
              const lapack_complex_float* ap, lapack_complex_float* afp,
        lapack_int* ipiv, const lapack_complex_float* b, lapack_int ldb,
                  lapack_complex_float* x, lapack_int ldx, float* rcond,
                                            float* ferr, float* berr);

    Fptr_NL_LAPACKE_chpsvx CHPSVX;

     /* LAPACKE CHPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chptrf) ( int matrix_layout , char uplo,
            lapack_int n , lapack_complex_float * a , lapack_int * ipiv );

    Fptr_NL_LAPACKE_chptrf CHPTRF;

    chpsvx_obj = new hpsvx_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    chpsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chpsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chpsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chpsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CHPSVX = (Fptr_NL_LAPACKE_chpsvx)dlsym(chpsvx_obj->hModule, "LAPACKE_chpsvx");
    ASSERT_TRUE(CHPSVX != NULL) << "failed to syt the Netlib LAPACKE_chpsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the hptrf API to compute the factorized A.  */
    if(chpsvx_obj->fact == 'F') {
        CHPTRF = (Fptr_NL_LAPACKE_chptrf)dlsym(chpsvx_obj->hModule,"LAPACKE_chptrf");
        ASSERT_TRUE(CHPTRF != NULL) << "failed to syt the Netlib LAPACKE_chptrf symbol";
            
        chpsvx_obj->inforef = CHPTRF( chpsvx_obj->matrix_layout,
                                      chpsvx_obj->uplo, chpsvx_obj->n,
                                      chpsvx_obj->afref,
                                      chpsvx_obj->ipivref);
                               
        chpsvx_obj->info = LAPACKE_chptrf( chpsvx_obj->matrix_layout,
                                           chpsvx_obj->uplo, chpsvx_obj->n,
                                           chpsvx_obj->af,
                                           chpsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    chpsvx_obj->inforef = CHPSVX( chpsvx_obj->matrix_layout, chpsvx_obj->fact,
                                  chpsvx_obj->uplo, chpsvx_obj->n,
                                  chpsvx_obj->nrhs,
                                  chpsvx_obj->aref,
                                  chpsvx_obj->afref,
                                  chpsvx_obj->ipivref,
                                  chpsvx_obj->bref, chpsvx_obj->ldb,
                                  chpsvx_obj->xref, chpsvx_obj->ldx,
                                  &chpsvx_obj->rcondref, 
                                  chpsvx_obj->ferrref,
                                  chpsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    chpsvx_obj->info = LAPACKE_chpsvx( chpsvx_obj->matrix_layout, chpsvx_obj->fact,
                                  chpsvx_obj->uplo, chpsvx_obj->n,
                                  chpsvx_obj->nrhs,
                                  chpsvx_obj->a,
                                  chpsvx_obj->af,
                                  chpsvx_obj->ipiv,
                                  chpsvx_obj->b, chpsvx_obj->ldb,
                                  chpsvx_obj->x, chpsvx_obj->ldx,
                                  &chpsvx_obj->rcond, 
                                  chpsvx_obj->ferr,
                                  chpsvx_obj->berr);

    if( chpsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chpsvx is wrong\n", chpsvx_obj->info );
    }
    if( chpsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chpsvx is wrong\n", 
        chpsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chpsvx_obj->ipiv_diff = computeDiff_i( chpsvx_obj->n, chpsvx_obj->ipiv, chpsvx_obj->ipivref );
    
    chpsvx_obj->diff =  computeDiff_c( chpsvx_obj->a_bufsize, 
                chpsvx_obj->af, chpsvx_obj->afref );

    chpsvx_obj->diff_xerr =  computeDiff_c( chpsvx_obj->x_bufsize, 
                chpsvx_obj->x, chpsvx_obj->xref );

    chpsvx_obj->diff_berr =  computeDiff_s( chpsvx_obj->nrhs, 
                chpsvx_obj->berr, chpsvx_obj->berrref );
                
    chpsvx_obj->diff_ferr =  computeDiff_s( chpsvx_obj->nrhs, 
                chpsvx_obj->ferr, chpsvx_obj->ferrref );
}

TEST_F(chpsvx_test, chpsvx1) {
    EXPECT_NEAR(0.0, chpsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (chpsvx_obj->rcond - chpsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, chpsvx_obj->ipiv_diff );
}

TEST_F(chpsvx_test, chpsvx2) {
    EXPECT_NEAR(0.0, chpsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (chpsvx_obj->rcond - chpsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpsvx_test, chpsvx3) {
    EXPECT_NEAR(0.0, chpsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (chpsvx_obj->rcond - chpsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpsvx_test, chpsvx4) {
    EXPECT_NEAR(0.0, chpsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chpsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (chpsvx_obj->rcond - chpsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin hpsvx_dcomplex_parameters  class definition */
class hpsvx_dcomplex_parameters{
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
      hpsvx_dcomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~hpsvx_dcomplex_parameters (); 
};  /* end of hpsvx_dcomplex_parameters  class definition */


/* Constructor hpsvx_dcomplex_parameters definition */
hpsvx_dcomplex_parameters:: hpsvx_dcomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n hpsvx lapack_complex_double:  n: %d, fact: %c uplo: %c  lda: %d  \
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
       hpsvx_free();
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

hpsvx_dcomplex_parameters:: ~hpsvx_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hpsvx_dcomplex_parameters object: destructor invoked. \n");
#endif
   hpsvx_free();
}

//  Test fixture class definition
class zhpsvx_test  : public  ::testing::Test {
public:
   hpsvx_dcomplex_parameters  *zhpsvx_obj;
   void SetUp();  
   void TearDown () { delete zhpsvx_obj; }
};


void zhpsvx_test::SetUp(){

    /* LAPACKE ZHPSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zhpsvx) (int matrix_layout, char fact,
                               char uplo, lapack_int n, lapack_int nrhs,
              const lapack_complex_double* ap, lapack_complex_double* afp,
        lapack_int* ipiv, const lapack_complex_double* b, lapack_int ldb,
                  lapack_complex_double* x, lapack_int ldx, double* rcond,
                                            double* ferr, double* berr);

    Fptr_NL_LAPACKE_zhpsvx ZHPSVX;

     /* LAPACKE ZHPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhptrf) ( int matrix_layout , char uplo,
            lapack_int n , lapack_complex_double * a , lapack_int * ipiv );

    Fptr_NL_LAPACKE_zhptrf ZHPTRF;

    zhpsvx_obj = new hpsvx_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zhpsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhpsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhpsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhpsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZHPSVX = (Fptr_NL_LAPACKE_zhpsvx)dlsym(zhpsvx_obj->hModule, "LAPACKE_zhpsvx");
    ASSERT_TRUE(ZHPSVX != NULL) << "failed to syt the Netlib LAPACKE_zhpsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the hptrf API to compute the factorized A.  */
    if(zhpsvx_obj->fact == 'F') {
        ZHPTRF = (Fptr_NL_LAPACKE_zhptrf)dlsym(zhpsvx_obj->hModule,"LAPACKE_zhptrf");
        ASSERT_TRUE(ZHPTRF != NULL) << "failed to syt the Netlib LAPACKE_zhptrf symbol";
            
        zhpsvx_obj->inforef = ZHPTRF( zhpsvx_obj->matrix_layout,
                                      zhpsvx_obj->uplo, zhpsvx_obj->n,
                                      zhpsvx_obj->afref,
                                      zhpsvx_obj->ipivref);
                               
        zhpsvx_obj->info = LAPACKE_zhptrf( zhpsvx_obj->matrix_layout,
                                           zhpsvx_obj->uplo, zhpsvx_obj->n,
                                           zhpsvx_obj->af,
                                           zhpsvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zhpsvx_obj->inforef = ZHPSVX( zhpsvx_obj->matrix_layout, zhpsvx_obj->fact,
                                  zhpsvx_obj->uplo, zhpsvx_obj->n,
                                  zhpsvx_obj->nrhs,
                                  zhpsvx_obj->aref,
                                  zhpsvx_obj->afref,
                                  zhpsvx_obj->ipivref,
                                  zhpsvx_obj->bref, zhpsvx_obj->ldb,
                                  zhpsvx_obj->xref, zhpsvx_obj->ldx,
                                  &zhpsvx_obj->rcondref, 
                                  zhpsvx_obj->ferrref,
                                  zhpsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zhpsvx_obj->info = LAPACKE_zhpsvx( zhpsvx_obj->matrix_layout, zhpsvx_obj->fact,
                                  zhpsvx_obj->uplo, zhpsvx_obj->n,
                                  zhpsvx_obj->nrhs,
                                  zhpsvx_obj->a,
                                  zhpsvx_obj->af,
                                  zhpsvx_obj->ipiv,
                                  zhpsvx_obj->b, zhpsvx_obj->ldb,
                                  zhpsvx_obj->x, zhpsvx_obj->ldx,
                                  &zhpsvx_obj->rcond, 
                                  zhpsvx_obj->ferr,
                                  zhpsvx_obj->berr);

    if( zhpsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhpsvx is wrong\n", zhpsvx_obj->info );
    }
    if( zhpsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhpsvx is wrong\n", 
        zhpsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhpsvx_obj->ipiv_diff = computeDiff_i( zhpsvx_obj->n, zhpsvx_obj->ipiv, zhpsvx_obj->ipivref );
    
    zhpsvx_obj->diff =  computeDiff_z( zhpsvx_obj->a_bufsize, 
                zhpsvx_obj->af, zhpsvx_obj->afref );

    zhpsvx_obj->diff_xerr =  computeDiff_z( zhpsvx_obj->x_bufsize, 
                zhpsvx_obj->x, zhpsvx_obj->xref );

    zhpsvx_obj->diff_berr =  computeDiff_d( zhpsvx_obj->nrhs, 
                zhpsvx_obj->berr, zhpsvx_obj->berrref );
                
    zhpsvx_obj->diff_ferr =  computeDiff_d( zhpsvx_obj->nrhs, 
                zhpsvx_obj->ferr, zhpsvx_obj->ferrref );
}

TEST_F(zhpsvx_test, zhpsvx1) {
    EXPECT_NEAR(0.0, zhpsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zhpsvx_obj->rcond - zhpsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, zhpsvx_obj->ipiv_diff );
}

TEST_F(zhpsvx_test, zhpsvx2) {
    EXPECT_NEAR(0.0, zhpsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zhpsvx_obj->rcond - zhpsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpsvx_test, zhpsvx3) {
    EXPECT_NEAR(0.0, zhpsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zhpsvx_obj->rcond - zhpsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpsvx_test, zhpsvx4) {
    EXPECT_NEAR(0.0, zhpsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhpsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zhpsvx_obj->rcond - zhpsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}