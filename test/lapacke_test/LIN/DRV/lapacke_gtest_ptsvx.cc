#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define ptsvx_free() \
  if (d != NULL)    free (d  ); \
  if (dref != NULL) free (dref); \
  if (df != NULL)    free (df   ); \
  if (dfref != NULL) free (dfref); \
  if (e != NULL)     free (e  ); \
  if (eref != NULL)  free (eref); \
  if (ef != NULL)    free (ef   ); \
  if (efref != NULL) free (efref); \
  if (b != NULL)    free (b  ); \
  if (bref != NULL) free (bref); \
  if (x != NULL)    free (x  ); \
  if (xref != NULL) free (xref); \
  if (ferr != NULL)    free (ferr  ); \
  if (ferrref != NULL) free (ferrref); \
  if (berr != NULL)    free (berr  ); \
  if (berrref != NULL) free (berrref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ptsvx_float_parameters  class definition */
class ptsvx_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff_df, diff_ef; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      float *d, *dref; // n diagonal elements of the tridiagonal matrix A
      float *e, *eref; // (n - 1) subdiagonal elements of the tridiagonal matrix A.
      float *df, *dfref; // n diagonal elements of the diagonal matrix D
      float* ef, *efref; // (n - 1) subdiagonal elements of the unit bidiagonal factor
      float *b, *bref; //right-hand sides for the systems of equations.
  
      /* Output parameters */
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      float* x, *xref; //  solution matrix X to the system of equations
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ptsvx_float_parameters ( int matrix_layout_i, char fact_i,
                               lapack_int n_i, lapack_int nrhs_i);
          
      ~ptsvx_float_parameters (); 
};  /* end of ptsvx_float_parameters  class definition */


/* Constructor ptsvx_float_parameters definition */
ptsvx_float_parameters:: ptsvx_float_parameters ( int matrix_layout_i, 
                char fact_i, lapack_int n_i, lapack_int nrhs_i) {
                                
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    n = n_i;
    nrhs = nrhs_i;

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
   printf(" \n ptsvx float:  matrix_layout: %d  n: %d, fact: %c ldb: %d nrhs: %d  \
ldx: %d \n",  matrix_layout, n, fact, ldb, nrhs, ldx);
#endif
    diff_df = 0.0;
    diff_ef = 0.0;
    diff_berr = 0.0;
    diff_ferr = 0.0;
    diff_xerr = 0.0;
    hModule = NULL;
    dModule = NULL;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_float_buffer_pair( &df, &dfref, n);
    lapacke_gtest_alloc_float_buffer_pair( &e,  &eref,  n-1);
    lapacke_gtest_alloc_float_buffer_pair( &ef, &efref, n-1);
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);

    if( (d==NULL) || (dref==NULL) ||  \
        (df==NULL) || (dfref==NULL) || \
        (e==NULL) || (eref==NULL) || \
        (ef==NULL) || (efref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (b==NULL) || (bref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       ptsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand ( d, dref, n);
    memcpy(df, d, (n*sizeof(float)));
    memcpy(dfref, d, (n*sizeof(float)));

    lapacke_gtest_init_float_buffer_pair_rand ( e, eref, n-1);
    memcpy(ef, e, ( (n-1)*sizeof(float)));
    memcpy(efref, e, ( (n-1)*sizeof(float)));

    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize); 
    lapacke_gtest_init_float_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Constructor  */

ptsvx_float_parameters:: ~ptsvx_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptsvx_float_parameters object: destructor invoked. \n");
#endif
   ptsvx_free();
}

//  Test fixture class definition
class sptsvx_test  : public  ::testing::Test {
public:
   ptsvx_float_parameters  *sptsvx_obj;
   void SetUp();  
   void TearDown () { delete sptsvx_obj; }
};


void sptsvx_test::SetUp(){

    /* LAPACKE SPTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_sptsvx) (int matrix_layout, char fact, 
            lapack_int n, lapack_int nrhs, const float *d, const float *e, 
            float *df, float *ef, const float *b, lapack_int ldb, float *x,
            lapack_int ldx, float *rcond, float *ferr, float *berr);

    Fptr_NL_LAPACKE_sptsvx SPTSVX;

     /* LAPACKE SPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spttrf) ( lapack_int n, float* d, float* e );

    Fptr_NL_LAPACKE_spttrf SPTTRF;

    sptsvx_obj = new ptsvx_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    sptsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sptsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sptsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sptsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SPTSVX = (Fptr_NL_LAPACKE_sptsvx)dlsym(sptsvx_obj->hModule, "LAPACKE_sptsvx");
    ASSERT_TRUE(SPTSVX != NULL) << "failed to ptt the Netlib LAPACKE_sptsvx symbol";

    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pttrf API to compute the factorized A.  */
    if(sptsvx_obj->fact == 'F') {
        SPTTRF = (Fptr_NL_LAPACKE_spttrf)dlsym(sptsvx_obj->hModule,"LAPACKE_spttrf");
        ASSERT_TRUE(SPTTRF != NULL) << "failed to ptt the Netlib LAPACKE_spttrf symbol";
        
        sptsvx_obj->inforef = SPTTRF( sptsvx_obj->n,
                                      sptsvx_obj->dfref,
                                      sptsvx_obj->efref);
                           
        sptsvx_obj->info = LAPACKE_spttrf( sptsvx_obj->n,
                                           sptsvx_obj->df,
                                           sptsvx_obj->ef);
    }
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    sptsvx_obj->inforef = SPTSVX( sptsvx_obj->matrix_layout, sptsvx_obj->fact,
                                  sptsvx_obj->n,
                                  sptsvx_obj->nrhs,
                                  (const float *)sptsvx_obj->dref,
                                  (const float *)sptsvx_obj->eref, 
                                  sptsvx_obj->dfref, sptsvx_obj->efref,
                                  (const float *)sptsvx_obj->bref, 
                                  sptsvx_obj->ldb,
                                  sptsvx_obj->xref, sptsvx_obj->ldx,
                                  &sptsvx_obj->rcondref, 
                                  sptsvx_obj->ferrref,
                                  sptsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    sptsvx_obj->info = LAPACKE_sptsvx( sptsvx_obj->matrix_layout, sptsvx_obj->fact,
                                  sptsvx_obj->n,
                                  sptsvx_obj->nrhs,
                                  (const float *)sptsvx_obj->d,
                                  (const float *)sptsvx_obj->e, 
                                  sptsvx_obj->df, sptsvx_obj->ef,
                                  (const float *)sptsvx_obj->b,
                                  sptsvx_obj->ldb,
                                  sptsvx_obj->x, sptsvx_obj->ldx,
                                  &sptsvx_obj->rcond, 
                                  sptsvx_obj->ferr,
                                  sptsvx_obj->berr);

    if( sptsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sptsvx is wrong\n", sptsvx_obj->info );
    }
    if( sptsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sptsvx is wrong\n", 
        sptsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    if( sptsvx_obj->fact == 'N'){
       sptsvx_obj->diff_df =  computeDiff_s( sptsvx_obj->n, 
                   sptsvx_obj->df, sptsvx_obj->dfref );
   
       sptsvx_obj->diff_ef =  computeDiff_s( sptsvx_obj->n-1, 
                   sptsvx_obj->ef, sptsvx_obj->efref );
    }

    sptsvx_obj->diff_xerr =  computeDiff_s( sptsvx_obj->x_bufsize, 
                sptsvx_obj->x, sptsvx_obj->xref );

    sptsvx_obj->diff_berr =  computeDiff_s( sptsvx_obj->nrhs, 
                sptsvx_obj->berr, sptsvx_obj->berrref );
            
    sptsvx_obj->diff_ferr =  computeDiff_s( sptsvx_obj->nrhs, 
                sptsvx_obj->ferr, sptsvx_obj->ferrref );
}

TEST_F(sptsvx_test, sptsvx1) {
    EXPECT_NEAR(0.0, sptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sptsvx_obj->rcond - sptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sptsvx_test, sptsvx2) {
    EXPECT_NEAR(0.0, sptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sptsvx_obj->rcond - sptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sptsvx_test, sptsvx3) {
    EXPECT_NEAR(0.0, sptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sptsvx_obj->rcond - sptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sptsvx_test, sptsvx4) {
    EXPECT_NEAR(0.0, sptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sptsvx_obj->rcond - sptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin ptsvx_double_parameters  class definition */
class ptsvx_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff_df, diff_ef; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      double *d, *dref; // n diagonal elements of the tridiagonal matrix A
      double *e, *eref; // (n - 1) subdiagonal elements of the tridiagonal matrix A.
      double *df, *dfref; // n diagonal elements of the diagonal matrix D
      double* ef, *efref; // (n - 1) subdiagonal elements of the unit bidiagonal factor
      double *b, *bref; //right-hand sides for the systems of equations.
  
      /* Output parameters */
      double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      double* x, *xref; //  solution matrix X to the system of equations
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ptsvx_double_parameters ( int matrix_layout_i, char fact_i,
                               lapack_int n_i, lapack_int nrhs_i);
          
      ~ptsvx_double_parameters (); 
};  /* end of ptsvx_double_parameters  class definition */


/* Constructor ptsvx_double_parameters definition */
ptsvx_double_parameters:: ptsvx_double_parameters ( int matrix_layout_i, 
                char fact_i, lapack_int n_i, lapack_int nrhs_i) {
                                
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    n = n_i;
    nrhs = nrhs_i;

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
   printf(" \n ptsvx double:  matrix_layout: %d  n: %d, fact: %c ldb: %d nrhs: %d  \
ldx: %d \n",  matrix_layout, n, fact, ldb, nrhs, ldx);
#endif
    diff_df = 0.0;
    diff_ef = 0.0;
    diff_berr = 0.0;
    diff_ferr = 0.0;
    diff_xerr = 0.0;
    hModule = NULL;
    dModule = NULL;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_double_buffer_pair( &df, &dfref, n);
    lapacke_gtest_alloc_double_buffer_pair( &e,  &eref,  n-1);
    lapacke_gtest_alloc_double_buffer_pair( &ef, &efref, n-1);
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);

    if( (d==NULL) || (dref==NULL) ||  \
        (df==NULL) || (dfref==NULL) || \
        (e==NULL) || (eref==NULL) || \
        (ef==NULL) || (efref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (b==NULL) || (bref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       ptsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand ( d, dref, n);
    memcpy(df, d, (n*sizeof(double)));
    memcpy(dfref, d, (n*sizeof(double)));

    lapacke_gtest_init_double_buffer_pair_rand ( e, eref, n-1);
    memcpy(ef, e, ( (n-1)*sizeof(double)));
    memcpy(efref, e, ( (n-1)*sizeof(double)));

    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize); 
    lapacke_gtest_init_double_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Constructor  */

ptsvx_double_parameters:: ~ptsvx_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptsvx_double_parameters object: destructor invoked. \n");
#endif
   ptsvx_free();
}

//  Test fixture class definition
class dptsvx_test  : public  ::testing::Test {
public:
   ptsvx_double_parameters  *dptsvx_obj;
   void SetUp();  
   void TearDown () { delete dptsvx_obj; }
};


void dptsvx_test::SetUp(){

    /* LAPACKE DPTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_dptsvx) (int matrix_layout, char fact, 
            lapack_int n, lapack_int nrhs, const double *d, const double *e, 
            double *df, double *ef, const double *b, lapack_int ldb, double *x,
            lapack_int ldx, double *rcond, double *ferr, double *berr);

    Fptr_NL_LAPACKE_dptsvx DPTSVX;

     /* LAPACKE DPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpttrf) ( lapack_int n, double* d, double* e );

    Fptr_NL_LAPACKE_dpttrf DPTTRF;

    dptsvx_obj = new ptsvx_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dptsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dptsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dptsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dptsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DPTSVX = (Fptr_NL_LAPACKE_dptsvx)dlsym(dptsvx_obj->hModule, "LAPACKE_dptsvx");
    ASSERT_TRUE(DPTSVX != NULL) << "failed to ptt the Netlib LAPACKE_dptsvx symbol";

    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pttrf API to compute the factorized A.  */
    if(dptsvx_obj->fact == 'F') {
        DPTTRF = (Fptr_NL_LAPACKE_dpttrf)dlsym(dptsvx_obj->hModule,"LAPACKE_dpttrf");
        ASSERT_TRUE(DPTTRF != NULL) << "failed to ptt the Netlib LAPACKE_dpttrf symbol";
        
        dptsvx_obj->inforef = DPTTRF( dptsvx_obj->n,
                                      dptsvx_obj->dfref,
                                      dptsvx_obj->efref);
                           
        dptsvx_obj->info = LAPACKE_dpttrf( dptsvx_obj->n,
                                           dptsvx_obj->df,
                                           dptsvx_obj->ef);
    }
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dptsvx_obj->inforef = DPTSVX( dptsvx_obj->matrix_layout, dptsvx_obj->fact,
                                  dptsvx_obj->n,
                                  dptsvx_obj->nrhs,
                                  (const double *)dptsvx_obj->dref,
                                  (const double *)dptsvx_obj->eref, 
                                  dptsvx_obj->dfref, dptsvx_obj->efref,
                                  (const double *)dptsvx_obj->bref, 
                                  dptsvx_obj->ldb,
                                  dptsvx_obj->xref, dptsvx_obj->ldx,
                                  &dptsvx_obj->rcondref, 
                                  dptsvx_obj->ferrref,
                                  dptsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dptsvx_obj->info = LAPACKE_dptsvx( dptsvx_obj->matrix_layout, dptsvx_obj->fact,
                                  dptsvx_obj->n,
                                  dptsvx_obj->nrhs,
                                  (const double *)dptsvx_obj->d,
                                  (const double *)dptsvx_obj->e, 
                                  dptsvx_obj->df, dptsvx_obj->ef,
                                  (const double *)dptsvx_obj->b,
                                  dptsvx_obj->ldb,
                                  dptsvx_obj->x, dptsvx_obj->ldx,
                                  &dptsvx_obj->rcond, 
                                  dptsvx_obj->ferr,
                                  dptsvx_obj->berr);

    if( dptsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dptsvx is wrong\n", dptsvx_obj->info );
    }
    if( dptsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dptsvx is wrong\n", 
        dptsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    if( dptsvx_obj->fact == 'N'){
       dptsvx_obj->diff_df =  computeDiff_d( dptsvx_obj->n, 
                   dptsvx_obj->df, dptsvx_obj->dfref );
   
       dptsvx_obj->diff_ef =  computeDiff_d( dptsvx_obj->n-1, 
                   dptsvx_obj->ef, dptsvx_obj->efref );
    }

    dptsvx_obj->diff_xerr =  computeDiff_d( dptsvx_obj->x_bufsize, 
                dptsvx_obj->x, dptsvx_obj->xref );

    dptsvx_obj->diff_berr =  computeDiff_d( dptsvx_obj->nrhs, 
                dptsvx_obj->berr, dptsvx_obj->berrref );
            
    dptsvx_obj->diff_ferr =  computeDiff_d( dptsvx_obj->nrhs, 
                dptsvx_obj->ferr, dptsvx_obj->ferrref );
}

TEST_F(dptsvx_test, dptsvx1) {
    EXPECT_NEAR(0.0, dptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dptsvx_obj->rcond - dptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dptsvx_test, dptsvx2) {
    EXPECT_NEAR(0.0, dptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dptsvx_obj->rcond - dptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dptsvx_test, dptsvx3) {
    EXPECT_NEAR(0.0, dptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dptsvx_obj->rcond - dptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dptsvx_test, dptsvx4) {
    EXPECT_NEAR(0.0, dptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dptsvx_obj->rcond - dptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin ptsvx_scomplex_parameters  class definition */
class ptsvx_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff_df, diff_ef; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      float *d, *dref; // n diagonal elements of the tridiagonal matrix A
      lapack_complex_float *e, *eref; // (n - 1) subdiagonal elements of the tridiagonal matrix A.
      float *df, *dfref; // n diagonal elements of the diagonal matrix D
      lapack_complex_float* ef, *efref; // (n - 1) subdiagonal elements of the unit bidiagonal factor
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.
  
      /* Output parameters */
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_float* x, *xref; //  solution matrix X to the system of equations
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ptsvx_scomplex_parameters ( int matrix_layout_i, char fact_i,
                               lapack_int n_i, lapack_int nrhs_i);
          
      ~ptsvx_scomplex_parameters (); 
};  /* end of ptsvx_scomplex_parameters  class definition */


/* Constructor ptsvx_scomplex_parameters definition */
ptsvx_scomplex_parameters:: ptsvx_scomplex_parameters ( int matrix_layout_i, 
                char fact_i, lapack_int n_i, lapack_int nrhs_i) {
                                
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    n = n_i;
    nrhs = nrhs_i;

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
   printf(" \n ptsvx lapack_complex_float:  matrix_layout: %d  n: %d, fact: %c ldb: %d nrhs: %d  \
ldx: %d \n",  matrix_layout, n, fact, ldb, nrhs, ldx);
#endif
    diff_df = 0.0;
    diff_ef = 0.0;
    diff_berr = 0.0;
    diff_ferr = 0.0;
    diff_xerr = 0.0;
    hModule = NULL;
    dModule = NULL;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e,  &eref,  n-1);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &ef, &efref, n-1);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_float_buffer_pair( &df, &dfref, n);
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);

    if( (d==NULL) || (dref==NULL) ||  \
        (df==NULL) || (dfref==NULL) || \
        (e==NULL) || (eref==NULL) || \
        (ef==NULL) || (efref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (b==NULL) || (bref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       ptsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand ( d, dref, n);
    memcpy(df, d, (n*sizeof(float)));
    memcpy(dfref, d, (n*sizeof(float)));

    lapacke_gtest_init_scomplex_buffer_pair_rand ( e, eref, n-1);
    memcpy(ef, e, ( (n-1)*sizeof(lapack_complex_float)));
    memcpy(efref, e, ( (n-1)*sizeof(lapack_complex_float)));

    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize); 
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Constructor  */

ptsvx_scomplex_parameters:: ~ptsvx_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptsvx_scomplex_parameters object: destructor invoked. \n");
#endif
   ptsvx_free();
}

//  Test fixture class definition
class cptsvx_test  : public  ::testing::Test {
public:
   ptsvx_scomplex_parameters  *cptsvx_obj;
   void SetUp();  
   void TearDown () { delete cptsvx_obj; }
};


void cptsvx_test::SetUp(){

    /* LAPACKE CPTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_cptsvx) (int matrix_layout, char fact,
                lapack_int n, lapack_int nrhs, const float* d,
                const lapack_complex_float* e, float* df, 
                lapack_complex_float* ef, const lapack_complex_float* b, 
                lapack_int ldb, lapack_complex_float* x, lapack_int ldx,
                float* rcond, float* ferr, float* berr);

    Fptr_NL_LAPACKE_cptsvx CPTSVX;

     /* LAPACKE CPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpttrf) ( lapack_int n, float* d, lapack_complex_float* e );

    Fptr_NL_LAPACKE_cpttrf CPTTRF;

    cptsvx_obj = new ptsvx_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    cptsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cptsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cptsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cptsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CPTSVX = (Fptr_NL_LAPACKE_cptsvx)dlsym(cptsvx_obj->hModule, "LAPACKE_cptsvx");
    ASSERT_TRUE(CPTSVX != NULL) << "failed to ptt the Netlib LAPACKE_cptsvx symbol";

    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pttrf API to compute the factorized A.  */
    if(cptsvx_obj->fact == 'F') {
        CPTTRF = (Fptr_NL_LAPACKE_cpttrf)dlsym(cptsvx_obj->hModule,"LAPACKE_cpttrf");
        ASSERT_TRUE(CPTTRF != NULL) << "failed to ptt the Netlib LAPACKE_cpttrf symbol";
        
        cptsvx_obj->inforef = CPTTRF( cptsvx_obj->n,
                                      cptsvx_obj->dfref,
                                      cptsvx_obj->efref);
                           
        cptsvx_obj->info = LAPACKE_cpttrf( cptsvx_obj->n,
                                           cptsvx_obj->df,
                                           cptsvx_obj->ef);
    }
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cptsvx_obj->inforef = CPTSVX( cptsvx_obj->matrix_layout, cptsvx_obj->fact,
                                  cptsvx_obj->n,
                                  cptsvx_obj->nrhs,
                                  (const float *)cptsvx_obj->dref,
                                  (const lapack_complex_float *)cptsvx_obj->eref, 
                                  cptsvx_obj->dfref, cptsvx_obj->efref,
                                  (const lapack_complex_float *)cptsvx_obj->bref, 
                                  cptsvx_obj->ldb,
                                  cptsvx_obj->xref, cptsvx_obj->ldx,
                                  &cptsvx_obj->rcondref, 
                                  cptsvx_obj->ferrref,
                                  cptsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cptsvx_obj->info = LAPACKE_cptsvx( cptsvx_obj->matrix_layout, cptsvx_obj->fact,
                                  cptsvx_obj->n,
                                  cptsvx_obj->nrhs,
                                  (const float *)cptsvx_obj->d,
                                  (const lapack_complex_float *)cptsvx_obj->e, 
                                  cptsvx_obj->df, cptsvx_obj->ef,
                                  (const lapack_complex_float *)cptsvx_obj->b,
                                  cptsvx_obj->ldb,
                                  cptsvx_obj->x, cptsvx_obj->ldx,
                                  &cptsvx_obj->rcond, 
                                  cptsvx_obj->ferr,
                                  cptsvx_obj->berr);

    if( cptsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cptsvx is wrong\n", cptsvx_obj->info );
    }
    if( cptsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cptsvx is wrong\n", 
        cptsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    if( cptsvx_obj->fact == 'N'){
       cptsvx_obj->diff_df =  computeDiff_s( cptsvx_obj->n, 
                   cptsvx_obj->df, cptsvx_obj->dfref );
   
       cptsvx_obj->diff_ef =  computeDiff_c( cptsvx_obj->n-1, 
                   cptsvx_obj->ef, cptsvx_obj->efref );
    }

    cptsvx_obj->diff_xerr =  computeDiff_c( cptsvx_obj->x_bufsize, 
                cptsvx_obj->x, cptsvx_obj->xref );

    cptsvx_obj->diff_berr =  computeDiff_s( cptsvx_obj->nrhs, 
                cptsvx_obj->berr, cptsvx_obj->berrref );
            
    cptsvx_obj->diff_ferr =  computeDiff_s( cptsvx_obj->nrhs, 
                cptsvx_obj->ferr, cptsvx_obj->ferrref );
}

TEST_F(cptsvx_test, cptsvx1) {
    EXPECT_NEAR(0.0, cptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cptsvx_obj->rcond - cptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cptsvx_test, cptsvx2) {
    EXPECT_NEAR(0.0, cptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cptsvx_obj->rcond - cptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cptsvx_test, cptsvx3) {
    EXPECT_NEAR(0.0, cptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cptsvx_obj->rcond - cptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cptsvx_test, cptsvx4) {
    EXPECT_NEAR(0.0, cptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cptsvx_obj->rcond - cptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin ptsvx_dcomplex_parameters  class definition */
class ptsvx_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff_df, diff_ef; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      double *d, *dref; // n diagonal elements of the tridiagonal matrix A
      lapack_complex_double *e, *eref; // (n - 1) subdiagonal elements of the tridiagonal matrix A.
      double *df, *dfref; // n diagonal elements of the diagonal matrix D
      lapack_complex_double* ef, *efref; // (n - 1) subdiagonal elements of the unit bidiagonal factor
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.
  
      /* Output parameters */
      double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_double* x, *xref; //  solution matrix X to the system of equations
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ptsvx_dcomplex_parameters ( int matrix_layout_i, char fact_i,
                               lapack_int n_i, lapack_int nrhs_i);
          
      ~ptsvx_dcomplex_parameters (); 
};  /* end of ptsvx_dcomplex_parameters  class definition */


/* Constructor ptsvx_dcomplex_parameters definition */
ptsvx_dcomplex_parameters:: ptsvx_dcomplex_parameters ( int matrix_layout_i, 
                char fact_i, lapack_int n_i, lapack_int nrhs_i) {
                                
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    n = n_i;
    nrhs = nrhs_i;

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
   printf(" \n ptsvx lapack_complex_double:  matrix_layout: %d  n: %d, fact: %c ldb: %d nrhs: %d  \
ldx: %d \n",  matrix_layout, n, fact, ldb, nrhs, ldx);
#endif
    diff_df = 0.0;
    diff_ef = 0.0;
    diff_berr = 0.0;
    diff_ferr = 0.0;
    diff_xerr = 0.0;
    hModule = NULL;
    dModule = NULL;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e,  &eref,  n-1);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &ef, &efref, n-1);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_double_buffer_pair( &df, &dfref, n);
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);

    if( (d==NULL) || (dref==NULL) ||  \
        (df==NULL) || (dfref==NULL) || \
        (e==NULL) || (eref==NULL) || \
        (ef==NULL) || (efref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (b==NULL) || (bref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       ptsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand ( d, dref, n);
    memcpy(df, d, (n*sizeof(double)));
    memcpy(dfref, d, (n*sizeof(double)));

    lapacke_gtest_init_dcomplex_buffer_pair_rand ( e, eref, n-1);
    memcpy(ef, e, ( (n-1)*sizeof(lapack_complex_double)));
    memcpy(efref, e, ( (n-1)*sizeof(lapack_complex_double)));

    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize); 
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Constructor  */

ptsvx_dcomplex_parameters:: ~ptsvx_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptsvx_dcomplex_parameters object: destructor invoked. \n");
#endif
   ptsvx_free();
}

//  Test fixture class definition
class zptsvx_test  : public  ::testing::Test {
public:
   ptsvx_dcomplex_parameters  *zptsvx_obj;
   void SetUp();  
   void TearDown () { delete zptsvx_obj; }
};


void zptsvx_test::SetUp(){

    /* LAPACKE ZPTSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zptsvx) (int matrix_layout, char fact,
                lapack_int n, lapack_int nrhs, const double* d,
                const lapack_complex_double* e, double* df, 
                lapack_complex_double* ef, const lapack_complex_double* b, 
                lapack_int ldb, lapack_complex_double* x, lapack_int ldx,
                double* rcond, double* ferr, double* berr);

    Fptr_NL_LAPACKE_zptsvx ZPTSVX;

     /* LAPACKE ZPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpttrf) ( lapack_int n, double* d, lapack_complex_double* e );

    Fptr_NL_LAPACKE_zpttrf ZPTTRF;

    zptsvx_obj = new ptsvx_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_driver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zptsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zptsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zptsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zptsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZPTSVX = (Fptr_NL_LAPACKE_zptsvx)dlsym(zptsvx_obj->hModule, "LAPACKE_zptsvx");
    ASSERT_TRUE(ZPTSVX != NULL) << "failed to ptt the Netlib LAPACKE_zptsvx symbol";

    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pttrf API to compute the factorized A.  */
    if(zptsvx_obj->fact == 'F') {
        ZPTTRF = (Fptr_NL_LAPACKE_zpttrf)dlsym(zptsvx_obj->hModule,"LAPACKE_zpttrf");
        ASSERT_TRUE(ZPTTRF != NULL) << "failed to ptt the Netlib LAPACKE_zpttrf symbol";
        
        zptsvx_obj->inforef = ZPTTRF( zptsvx_obj->n,
                                      zptsvx_obj->dfref,
                                      zptsvx_obj->efref);
                           
        zptsvx_obj->info = LAPACKE_zpttrf( zptsvx_obj->n,
                                           zptsvx_obj->df,
                                           zptsvx_obj->ef);
    }
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zptsvx_obj->inforef = ZPTSVX( zptsvx_obj->matrix_layout, zptsvx_obj->fact,
                                  zptsvx_obj->n,
                                  zptsvx_obj->nrhs,
                                  (const double *)zptsvx_obj->dref,
                                  (const lapack_complex_double *)zptsvx_obj->eref, 
                                  zptsvx_obj->dfref, zptsvx_obj->efref,
                                  (const lapack_complex_double *)zptsvx_obj->bref, 
                                  zptsvx_obj->ldb,
                                  zptsvx_obj->xref, zptsvx_obj->ldx,
                                  &zptsvx_obj->rcondref, 
                                  zptsvx_obj->ferrref,
                                  zptsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zptsvx_obj->info = LAPACKE_zptsvx( zptsvx_obj->matrix_layout, zptsvx_obj->fact,
                                  zptsvx_obj->n,
                                  zptsvx_obj->nrhs,
                                  (const double *)zptsvx_obj->d,
                                  (const lapack_complex_double *)zptsvx_obj->e, 
                                  zptsvx_obj->df, zptsvx_obj->ef,
                                  (const lapack_complex_double *)zptsvx_obj->b,
                                  zptsvx_obj->ldb,
                                  zptsvx_obj->x, zptsvx_obj->ldx,
                                  &zptsvx_obj->rcond, 
                                  zptsvx_obj->ferr,
                                  zptsvx_obj->berr);

    if( zptsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zptsvx is wrong\n", zptsvx_obj->info );
    }
    if( zptsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zptsvx is wrong\n", 
        zptsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    if( zptsvx_obj->fact == 'N'){
       zptsvx_obj->diff_df =  computeDiff_d( zptsvx_obj->n, 
                   zptsvx_obj->df, zptsvx_obj->dfref );
   
       zptsvx_obj->diff_ef =  computeDiff_z( zptsvx_obj->n-1, 
                   zptsvx_obj->ef, zptsvx_obj->efref );
    }

    zptsvx_obj->diff_xerr =  computeDiff_z( zptsvx_obj->x_bufsize, 
                zptsvx_obj->x, zptsvx_obj->xref );

    zptsvx_obj->diff_berr =  computeDiff_d( zptsvx_obj->nrhs, 
                zptsvx_obj->berr, zptsvx_obj->berrref );
            
    zptsvx_obj->diff_ferr =  computeDiff_d( zptsvx_obj->nrhs, 
                zptsvx_obj->ferr, zptsvx_obj->ferrref );
}

TEST_F(zptsvx_test, zptsvx1) {
    EXPECT_NEAR(0.0, zptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zptsvx_obj->rcond - zptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zptsvx_test, zptsvx2) {
    EXPECT_NEAR(0.0, zptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zptsvx_obj->rcond - zptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zptsvx_test, zptsvx3) {
    EXPECT_NEAR(0.0, zptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zptsvx_obj->rcond - zptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zptsvx_test, zptsvx4) {
    EXPECT_NEAR(0.0, zptsvx_obj->diff_df, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_ef, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zptsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zptsvx_obj->rcond - zptsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}


