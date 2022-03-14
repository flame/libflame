#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define hesvx_free() \
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

/* Begin hesvx_scomplex_parameters  class definition */
class hesvx_scomplex_parameters{
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
      hesvx_scomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~hesvx_scomplex_parameters (); 
};  /* end of hesvx_scomplex_parameters  class definition */


/* Constructor hesvx_scomplex_parameters definition */
hesvx_scomplex_parameters:: hesvx_scomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n hesvx lapack_complex_float:  n: %d, fact: %c uplo: %c  lda: %d  \
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
       hesvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, 'H');
    memcpy(af, a, (n*n*sizeof(lapack_complex_float)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

hesvx_scomplex_parameters:: ~hesvx_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hesvx_scomplex_parameters object: destructor invoked. \n");
#endif
   hesvx_free();
}

//  Test fixture class definition
class chesvx_test  : public  ::testing::Test {
public:
   hesvx_scomplex_parameters  *chesvx_obj;
   void SetUp();  
   void TearDown () { delete chesvx_obj; }
};


void chesvx_test::SetUp(){

    /* LAPACKE CHESVX prototype */
    typedef int (*Fptr_NL_LAPACKE_chesvx) (int matrix_layout, char fact,
                char uplo, lapack_int n, lapack_int nrhs, 
                const lapack_complex_float* a, lapack_int lda, 
                lapack_complex_float* af, lapack_int ldaf, lapack_int* ipiv,
                const lapack_complex_float* b, lapack_int ldb,
                lapack_complex_float* x, lapack_int ldx,
                float* rcond, float* ferr, float* berr);

    Fptr_NL_LAPACKE_chesvx CHESVX;

     /* LAPACKE CHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf) ( int matrix_layout , char uplo,
            lapack_int n , lapack_complex_float * a , lapack_int lda , lapack_int * ipiv );

    Fptr_NL_LAPACKE_chetrf CHETRF;

    chesvx_obj = new hesvx_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    chesvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chesvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chesvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chesvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CHESVX = (Fptr_NL_LAPACKE_chesvx)dlsym(chesvx_obj->hModule, "LAPACKE_chesvx");
    ASSERT_TRUE(CHESVX != NULL) << "failed to syt the Netlib LAPACKE_chesvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the hetrf API to compute the factorized A.  */
    if(chesvx_obj->fact == 'F') {
        CHETRF = (Fptr_NL_LAPACKE_chetrf)dlsym(chesvx_obj->hModule,"LAPACKE_chetrf");
        ASSERT_TRUE(CHETRF != NULL) << "failed to syt the Netlib LAPACKE_chetrf symbol";
            
        chesvx_obj->inforef = CHETRF( chesvx_obj->matrix_layout,
                                      chesvx_obj->uplo, chesvx_obj->n,
                                      chesvx_obj->afref,
                                      chesvx_obj->lda,
                                      chesvx_obj->ipivref);
                               
        chesvx_obj->info = LAPACKE_chetrf( chesvx_obj->matrix_layout,
                                           chesvx_obj->uplo, chesvx_obj->n,
                                           chesvx_obj->af,
                                           chesvx_obj->lda,
                                           chesvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    chesvx_obj->inforef = CHESVX( chesvx_obj->matrix_layout, chesvx_obj->fact,
                                  chesvx_obj->uplo, chesvx_obj->n,
                                  chesvx_obj->nrhs,
                                  chesvx_obj->aref, chesvx_obj->lda, 
                                  chesvx_obj->afref, chesvx_obj->ldaf,
                                  chesvx_obj->ipivref,
                                  chesvx_obj->bref, chesvx_obj->ldb,
                                  chesvx_obj->xref, chesvx_obj->ldx,
                                  &chesvx_obj->rcondref, 
                                  chesvx_obj->ferrref,
                                  chesvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    chesvx_obj->info = LAPACKE_chesvx( chesvx_obj->matrix_layout, chesvx_obj->fact,
                                  chesvx_obj->uplo, chesvx_obj->n,
                                  chesvx_obj->nrhs,
                                  chesvx_obj->a, chesvx_obj->lda, 
                                  chesvx_obj->af, chesvx_obj->ldaf,
                                  chesvx_obj->ipiv,
                                  chesvx_obj->b, chesvx_obj->ldb,
                                  chesvx_obj->x, chesvx_obj->ldx,
                                  &chesvx_obj->rcond, 
                                  chesvx_obj->ferr,
                                  chesvx_obj->berr);

    if( chesvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chesvx is wrong\n", chesvx_obj->info );
    }
    if( chesvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chesvx is wrong\n", 
        chesvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chesvx_obj->ipiv_diff = computeDiff_i( chesvx_obj->n, chesvx_obj->ipiv, chesvx_obj->ipivref );
    
    chesvx_obj->diff =  computeDiff_c( chesvx_obj->n * chesvx_obj->n, 
                chesvx_obj->af, chesvx_obj->afref );

    chesvx_obj->diff_xerr =  computeDiff_c( chesvx_obj->x_bufsize, 
                chesvx_obj->x, chesvx_obj->xref );

    chesvx_obj->diff_berr =  computeDiff_s( chesvx_obj->nrhs, 
                chesvx_obj->berr, chesvx_obj->berrref );
                
    chesvx_obj->diff_ferr =  computeDiff_s( chesvx_obj->nrhs, 
                chesvx_obj->ferr, chesvx_obj->ferrref );
}

TEST_F(chesvx_test, chesvx1) {
    EXPECT_NEAR(0.0, chesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (chesvx_obj->rcond - chesvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, chesvx_obj->ipiv_diff );
}

TEST_F(chesvx_test, chesvx2) {
    EXPECT_NEAR(0.0, chesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (chesvx_obj->rcond - chesvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chesvx_test, chesvx3) {
    EXPECT_NEAR(0.0, chesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (chesvx_obj->rcond - chesvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chesvx_test, chesvx4) {
    EXPECT_NEAR(0.0, chesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (chesvx_obj->rcond - chesvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin hesvx_dcomplex_parameters  class definition */
class hesvx_dcomplex_parameters{
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
      hesvx_dcomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~hesvx_dcomplex_parameters (); 
};  /* end of hesvx_dcomplex_parameters  class definition */


/* Constructor hesvx_dcomplex_parameters definition */
hesvx_dcomplex_parameters:: hesvx_dcomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n hesvx lapack_complex_double:  n: %d, fact: %c uplo: %c  lda: %d  \
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
       hesvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, 'H');
    memcpy(af, a, (n*n*sizeof(lapack_complex_double)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

hesvx_dcomplex_parameters:: ~hesvx_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hesvx_dcomplex_parameters object: destructor invoked. \n");
#endif
   hesvx_free();
}

//  Test fixture class definition
class zhesvx_test  : public  ::testing::Test {
public:
   hesvx_dcomplex_parameters  *zhesvx_obj;
   void SetUp();  
   void TearDown () { delete zhesvx_obj; }
};


void zhesvx_test::SetUp(){

    /* LAPACKE ZHESVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zhesvx) (int matrix_layout, char fact,
                char uplo, lapack_int n, lapack_int nrhs, 
                const lapack_complex_double* a, lapack_int lda, 
                lapack_complex_double* af, lapack_int ldaf, lapack_int* ipiv,
                const lapack_complex_double* b, lapack_int ldb,
                lapack_complex_double* x, lapack_int ldx,
                double* rcond, double* ferr, double* berr);

    Fptr_NL_LAPACKE_zhesvx ZHESVX;

     /* LAPACKE ZHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf) ( int matrix_layout , char uplo,
            lapack_int n , lapack_complex_double * a , lapack_int lda , lapack_int * ipiv );

    Fptr_NL_LAPACKE_zhetrf ZHETRF;

    zhesvx_obj = new hesvx_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zhesvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhesvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhesvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhesvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZHESVX = (Fptr_NL_LAPACKE_zhesvx)dlsym(zhesvx_obj->hModule, "LAPACKE_zhesvx");
    ASSERT_TRUE(ZHESVX != NULL) << "failed to syt the Netlib LAPACKE_zhesvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the hetrf API to compute the factorized A.  */
    if(zhesvx_obj->fact == 'F') {
        ZHETRF = (Fptr_NL_LAPACKE_zhetrf)dlsym(zhesvx_obj->hModule,"LAPACKE_zhetrf");
        ASSERT_TRUE(ZHETRF != NULL) << "failed to syt the Netlib LAPACKE_zhetrf symbol";
            
        zhesvx_obj->inforef = ZHETRF( zhesvx_obj->matrix_layout,
                                      zhesvx_obj->uplo, zhesvx_obj->n,
                                      zhesvx_obj->afref,
                                      zhesvx_obj->lda,
                                      zhesvx_obj->ipivref);
                               
        zhesvx_obj->info = LAPACKE_zhetrf( zhesvx_obj->matrix_layout,
                                           zhesvx_obj->uplo, zhesvx_obj->n,
                                           zhesvx_obj->af,
                                           zhesvx_obj->lda,
                                           zhesvx_obj->ipiv);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zhesvx_obj->inforef = ZHESVX( zhesvx_obj->matrix_layout, zhesvx_obj->fact,
                                  zhesvx_obj->uplo, zhesvx_obj->n,
                                  zhesvx_obj->nrhs,
                                  zhesvx_obj->aref, zhesvx_obj->lda, 
                                  zhesvx_obj->afref, zhesvx_obj->ldaf,
                                  zhesvx_obj->ipivref,
                                  zhesvx_obj->bref, zhesvx_obj->ldb,
                                  zhesvx_obj->xref, zhesvx_obj->ldx,
                                  &zhesvx_obj->rcondref, 
                                  zhesvx_obj->ferrref,
                                  zhesvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zhesvx_obj->info = LAPACKE_zhesvx( zhesvx_obj->matrix_layout, zhesvx_obj->fact,
                                  zhesvx_obj->uplo, zhesvx_obj->n,
                                  zhesvx_obj->nrhs,
                                  zhesvx_obj->a, zhesvx_obj->lda, 
                                  zhesvx_obj->af, zhesvx_obj->ldaf,
                                  zhesvx_obj->ipiv,
                                  zhesvx_obj->b, zhesvx_obj->ldb,
                                  zhesvx_obj->x, zhesvx_obj->ldx,
                                  &zhesvx_obj->rcond, 
                                  zhesvx_obj->ferr,
                                  zhesvx_obj->berr);

    if( zhesvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhesvx is wrong\n", zhesvx_obj->info );
    }
    if( zhesvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhesvx is wrong\n", 
        zhesvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhesvx_obj->ipiv_diff = computeDiff_i( zhesvx_obj->n, zhesvx_obj->ipiv, zhesvx_obj->ipivref );
    
    zhesvx_obj->diff =  computeDiff_z( zhesvx_obj->n * zhesvx_obj->n, 
                zhesvx_obj->af, zhesvx_obj->afref );

    zhesvx_obj->diff_xerr =  computeDiff_z( zhesvx_obj->x_bufsize, 
                zhesvx_obj->x, zhesvx_obj->xref );

    zhesvx_obj->diff_berr =  computeDiff_d( zhesvx_obj->nrhs, 
                zhesvx_obj->berr, zhesvx_obj->berrref );
                
    zhesvx_obj->diff_ferr =  computeDiff_d( zhesvx_obj->nrhs, 
                zhesvx_obj->ferr, zhesvx_obj->ferrref );
}

TEST_F(zhesvx_test, zhesvx1) {
    EXPECT_NEAR(0.0, zhesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zhesvx_obj->rcond - zhesvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(0, zhesvx_obj->ipiv_diff );
}

TEST_F(zhesvx_test, zhesvx2) {
    EXPECT_NEAR(0.0, zhesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zhesvx_obj->rcond - zhesvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhesvx_test, zhesvx3) {
    EXPECT_NEAR(0.0, zhesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zhesvx_obj->rcond - zhesvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhesvx_test, zhesvx4) {
    EXPECT_NEAR(0.0, zhesvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhesvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zhesvx_obj->rcond - zhesvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}
