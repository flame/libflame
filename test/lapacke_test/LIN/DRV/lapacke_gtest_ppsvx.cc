
#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"
#define ppsvx_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (x != NULL)    free (x  ); \
  if (xref != NULL) free (xref); \
  if (af != NULL)    free (af  ); \
  if (afref != NULL) free (afref); \
  if (s != NULL)    free (s  ); \
  if (sref != NULL) free (sref); \
  if (ferr != NULL)    free (ferr  ); \
  if (ferrref != NULL) free (ferrref); \
  if (berr != NULL)    free (berr  ); \
  if (berrref != NULL) free (berrref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ppsvx_float_parameters  class definition */
class ppsvx_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr, diff_s;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.
      float *a, *aref; //The array ab contains the matrix A
      float *af, *afref; //contains the factored form of the matrix A
      float* s, *sref; // the row scale factors for A.
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ppsvx_float_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~ppsvx_float_parameters (); 
};  /* end of ppsvx_float_parameters  class definition */


/* Constructor ppsvx_float_parameters definition */
ppsvx_float_parameters:: ppsvx_float_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
    int a_buf_size =  n_i*(n_i+1)/2;    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
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
    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;
    diff_s = 0;
    diff = 0;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n ppsvx float:  n: %d, fact: %c uplo: %c   \
ldb: %d nrhs: %d   ldx: %d \n",  n, fact, uplo,  
                                          ldb, nrhs,  ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, a_buf_size);
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &af, &afref, a_buf_size);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (s==NULL) || (sref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       ppsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, a_buf_size);
    memcpy(af, a, (a_buf_size*sizeof(float)));
    memcpy(afref, a, (a_buf_size*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_float_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

ppsvx_float_parameters:: ~ppsvx_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppsvx_float_parameters object: destructor invoked. \n");
#endif
   ppsvx_free();
}

//  Test fixture class definition
class sppsvx_test  : public  ::testing::Test {
public:
   ppsvx_float_parameters  *sppsvx_obj;
   void SetUp();  
   void TearDown () { delete sppsvx_obj; }
};


void sppsvx_test::SetUp(){

    /* LAPACKE SPPSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_sppsvx) (int matrix_layout, char fact, 
                    char uplo, lapack_int n, lapack_int nrhs, float* ap,
                    float* afp, char* equed, float* s, float* b, 
                    lapack_int ldb, float* x, lapack_int ldx, 
                    float* rcond, float* ferr, float* berr);

    Fptr_NL_LAPACKE_sppsvx SPPSVX;

     /* LAPACKE SPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spptrf) ( int matrix_layout , char uplo , 
                                                lapack_int n , float * ap );

    Fptr_NL_LAPACKE_spptrf SPPTRF;

    sppsvx_obj = new ppsvx_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    sppsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sppsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sppsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sppsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SPPSVX = (Fptr_NL_LAPACKE_sppsvx)dlsym(sppsvx_obj->hModule, "LAPACKE_sppsvx");
    ASSERT_TRUE(SPPSVX != NULL) << "failed to ppt the Netlib LAPACKE_sppsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pptrf API to compute the factorized A.  */
    if(sppsvx_obj->fact == 'F') {
        SPPTRF = (Fptr_NL_LAPACKE_spptrf)dlsym(sppsvx_obj->hModule,"LAPACKE_spptrf");
        ASSERT_TRUE(SPPTRF != NULL) << "failed to ppt the Netlib LAPACKE_spptrf symbol";
            
        sppsvx_obj->inforef = SPPTRF( sppsvx_obj->matrix_layout,
                                      sppsvx_obj->uplo, sppsvx_obj->n,
                                      sppsvx_obj->afref);
                               
        sppsvx_obj->info = LAPACKE_spptrf( sppsvx_obj->matrix_layout,
                                           sppsvx_obj->uplo, sppsvx_obj->n,
                                           sppsvx_obj->af);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    sppsvx_obj->inforef = SPPSVX( sppsvx_obj->matrix_layout, sppsvx_obj->fact,
                                  sppsvx_obj->uplo, sppsvx_obj->n,
                                  sppsvx_obj->nrhs,
                                  sppsvx_obj->aref, 
                                  sppsvx_obj->afref,
                                  &sppsvx_obj->equedref,
                                  sppsvx_obj->sref,                               
                                  sppsvx_obj->bref, sppsvx_obj->ldb,
                                  sppsvx_obj->xref, sppsvx_obj->ldx,
                                  &sppsvx_obj->rcondref, 
                                  sppsvx_obj->ferrref,
                                  sppsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    sppsvx_obj->info = LAPACKE_sppsvx( sppsvx_obj->matrix_layout, sppsvx_obj->fact,
                                  sppsvx_obj->uplo, sppsvx_obj->n,
                                  sppsvx_obj->nrhs,
                                  sppsvx_obj->a, 
                                  sppsvx_obj->af,
                                  &sppsvx_obj->equed,
                                  sppsvx_obj->s,                                  
                                  sppsvx_obj->b, sppsvx_obj->ldb,
                                  sppsvx_obj->x, sppsvx_obj->ldx,
                                  &sppsvx_obj->rcond, 
                                  sppsvx_obj->ferr,
                                  sppsvx_obj->berr);

    if( sppsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sppsvx is wrong\n", sppsvx_obj->info );
    }
    if( sppsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sppsvx is wrong\n", 
        sppsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sppsvx_obj->diff =  computeDiff_s( sppsvx_obj->b_bufsize, 
                sppsvx_obj->b, sppsvx_obj->bref );

    sppsvx_obj->diff_xerr =  computeDiff_s( sppsvx_obj->x_bufsize, 
                sppsvx_obj->x, sppsvx_obj->xref );

    sppsvx_obj->diff_berr =  computeDiff_s( sppsvx_obj->nrhs, 
                sppsvx_obj->berr, sppsvx_obj->berrref );
                
    sppsvx_obj->diff_ferr =  computeDiff_s( sppsvx_obj->nrhs, 
                sppsvx_obj->ferr, sppsvx_obj->ferrref );
    
    if( sppsvx_obj->fact != 'F') {
       sppsvx_obj->diff_s =  computeDiff_s( sppsvx_obj->n, 
                sppsvx_obj->s, sppsvx_obj->sref );
    }
}

TEST_F(sppsvx_test, sppsvx1) {
    EXPECT_NEAR(0.0, sppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sppsvx_obj->rcond - sppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(sppsvx_obj->equed, sppsvx_obj->equedref);
}

TEST_F(sppsvx_test, sppsvx2) {
    EXPECT_NEAR(0.0, sppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sppsvx_obj->rcond - sppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(sppsvx_obj->equed, sppsvx_obj->equedref);
}

TEST_F(sppsvx_test, sppsvx3) {
    EXPECT_NEAR(0.0, sppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sppsvx_obj->rcond - sppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(sppsvx_obj->equed, sppsvx_obj->equedref);
}

TEST_F(sppsvx_test, sppsvx4) {
    EXPECT_NEAR(0.0, sppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sppsvx_obj->rcond - sppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(sppsvx_obj->equed, sppsvx_obj->equedref);
}

/* Begin ppsvx_double_parameters  class definition */
class ppsvx_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr, diff_s;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.
      double *a, *aref; //The array ab contains the matrix A
      double *af, *afref; //contains the factored form of the matrix A
      double* s, *sref; // the row scale factors for A.
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ppsvx_double_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~ppsvx_double_parameters (); 
};  /* end of ppsvx_double_parameters  class definition */


/* Constructor ppsvx_double_parameters definition */
ppsvx_double_parameters:: ppsvx_double_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
    int a_buf_size =  n_i*(n_i+1)/2;    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
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
    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;
    diff_s = 0;
    diff = 0;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n ppsvx double:  n: %d, fact: %c uplo: %c   \
ldb: %d nrhs: %d   ldx: %d \n",  n, fact, uplo,  
                                          ldb, nrhs,  ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, a_buf_size);
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &af, &afref, a_buf_size);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (s==NULL) || (sref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       ppsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, a_buf_size);
    memcpy(af, a, (a_buf_size*sizeof(double)));
    memcpy(afref, a, (a_buf_size*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_double_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

ppsvx_double_parameters:: ~ppsvx_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppsvx_double_parameters object: destructor invoked. \n");
#endif
   ppsvx_free();
}

//  Test fixture class definition
class dppsvx_test  : public  ::testing::Test {
public:
   ppsvx_double_parameters  *dppsvx_obj;
   void SetUp();  
   void TearDown () { delete dppsvx_obj; }
};


void dppsvx_test::SetUp(){

    /* LAPACKE DPPSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_dppsvx) (int matrix_layout, char fact, 
                    char uplo, lapack_int n, lapack_int nrhs, double* ap,
                    double* afp, char* equed, double* s, double* b, 
                    lapack_int ldb, double* x, lapack_int ldx, 
                    double* rcond, double* ferr, double* berr);

    Fptr_NL_LAPACKE_dppsvx DPPSVX;

     /* LAPACKE DPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpptrf) ( int matrix_layout , char uplo , 
                                                lapack_int n , double * ap );

    Fptr_NL_LAPACKE_dpptrf DPPTRF;

    dppsvx_obj = new ppsvx_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dppsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dppsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dppsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dppsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DPPSVX = (Fptr_NL_LAPACKE_dppsvx)dlsym(dppsvx_obj->hModule, "LAPACKE_dppsvx");
    ASSERT_TRUE(DPPSVX != NULL) << "failed to ppt the Netlib LAPACKE_dppsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pptrf API to compute the factorized A.  */
    if(dppsvx_obj->fact == 'F') {
        DPPTRF = (Fptr_NL_LAPACKE_dpptrf)dlsym(dppsvx_obj->hModule,"LAPACKE_dpptrf");
        ASSERT_TRUE(DPPTRF != NULL) << "failed to ppt the Netlib LAPACKE_dpptrf symbol";
            
        dppsvx_obj->inforef = DPPTRF( dppsvx_obj->matrix_layout,
                                      dppsvx_obj->uplo, dppsvx_obj->n,
                                      dppsvx_obj->afref);
                               
        dppsvx_obj->info = LAPACKE_dpptrf( dppsvx_obj->matrix_layout,
                                           dppsvx_obj->uplo, dppsvx_obj->n,
                                           dppsvx_obj->af);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dppsvx_obj->inforef = DPPSVX( dppsvx_obj->matrix_layout, dppsvx_obj->fact,
                                  dppsvx_obj->uplo, dppsvx_obj->n,
                                  dppsvx_obj->nrhs,
                                  dppsvx_obj->aref, 
                                  dppsvx_obj->afref,
                                  &dppsvx_obj->equedref,
                                  dppsvx_obj->sref,                               
                                  dppsvx_obj->bref, dppsvx_obj->ldb,
                                  dppsvx_obj->xref, dppsvx_obj->ldx,
                                  &dppsvx_obj->rcondref, 
                                  dppsvx_obj->ferrref,
                                  dppsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dppsvx_obj->info = LAPACKE_dppsvx( dppsvx_obj->matrix_layout, dppsvx_obj->fact,
                                  dppsvx_obj->uplo, dppsvx_obj->n,
                                  dppsvx_obj->nrhs,
                                  dppsvx_obj->a, 
                                  dppsvx_obj->af,
                                  &dppsvx_obj->equed,
                                  dppsvx_obj->s,                                  
                                  dppsvx_obj->b, dppsvx_obj->ldb,
                                  dppsvx_obj->x, dppsvx_obj->ldx,
                                  &dppsvx_obj->rcond, 
                                  dppsvx_obj->ferr,
                                  dppsvx_obj->berr);

    if( dppsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dppsvx is wrong\n", dppsvx_obj->info );
    }
    if( dppsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dppsvx is wrong\n", 
        dppsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dppsvx_obj->diff =  computeDiff_d( dppsvx_obj->b_bufsize, 
                dppsvx_obj->b, dppsvx_obj->bref );

    dppsvx_obj->diff_xerr =  computeDiff_d( dppsvx_obj->x_bufsize, 
                dppsvx_obj->x, dppsvx_obj->xref );

    dppsvx_obj->diff_berr =  computeDiff_d( dppsvx_obj->nrhs, 
                dppsvx_obj->berr, dppsvx_obj->berrref );
                
    dppsvx_obj->diff_ferr =  computeDiff_d( dppsvx_obj->nrhs, 
                dppsvx_obj->ferr, dppsvx_obj->ferrref );
    
    if( dppsvx_obj->fact != 'F') {
       dppsvx_obj->diff_s =  computeDiff_d( dppsvx_obj->n, 
                dppsvx_obj->s, dppsvx_obj->sref );
    }
}

TEST_F(dppsvx_test, dppsvx1) {
    EXPECT_NEAR(0.0, dppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dppsvx_obj->rcond - dppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(dppsvx_obj->equed, dppsvx_obj->equedref);
}

TEST_F(dppsvx_test, dppsvx2) {
    EXPECT_NEAR(0.0, dppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dppsvx_obj->rcond - dppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(dppsvx_obj->equed, dppsvx_obj->equedref);
}

TEST_F(dppsvx_test, dppsvx3) {
    EXPECT_NEAR(0.0, dppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dppsvx_obj->rcond - dppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(dppsvx_obj->equed, dppsvx_obj->equedref);
}

TEST_F(dppsvx_test, dppsvx4) {
    EXPECT_NEAR(0.0, dppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dppsvx_obj->rcond - dppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(dppsvx_obj->equed, dppsvx_obj->equedref);
}

/* Begin ppsvx_scomplex_parameters  class definition */
class ppsvx_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr, diff_s;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_float *a, *aref; //The array ab contains the matrix A
      lapack_complex_float *af, *afref; //contains the factored form of the matrix A
      float* s, *sref; // the row scale factors for A.
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ppsvx_scomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~ppsvx_scomplex_parameters (); 
};  /* end of ppsvx_scomplex_parameters  class definition */


/* Constructor ppsvx_scomplex_parameters definition */
ppsvx_scomplex_parameters:: ppsvx_scomplex_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
    int a_buf_size =  n_i*(n_i+1)/2;    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
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
    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;
    diff_s = 0;
    diff = 0;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n ppsvx float:  n: %d, fact: %c uplo: %c   \
ldb: %d nrhs: %d   ldx: %d \n",  n, fact, uplo,  
                                          ldb, nrhs,  ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_buf_size);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &af, &afref, a_buf_size);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (s==NULL) || (sref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       ppsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, a_buf_size);
    memcpy(af, a, (a_buf_size*sizeof(lapack_complex_float)));
    memcpy(afref, a, (a_buf_size*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

ppsvx_scomplex_parameters:: ~ppsvx_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppsvx_scomplex_parameters object: destructor invoked. \n");
#endif
   ppsvx_free();
}

//  Test fixture class definition
class cppsvx_test  : public  ::testing::Test {
public:
   ppsvx_scomplex_parameters  *cppsvx_obj;
   void SetUp();  
   void TearDown () { delete cppsvx_obj; }
};


void cppsvx_test::SetUp(){

    /* LAPACKE CPPSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_cppsvx) (int matrix_layout, char fact, 
                    char uplo, lapack_int n, lapack_int nrhs,
                    lapack_complex_float* ap, lapack_complex_float* afp,
                    char* equed, float* s, lapack_complex_float* b, 
                    lapack_int ldb, lapack_complex_float* x, lapack_int ldx, 
                    float* rcond, float* ferr, float* berr);

    Fptr_NL_LAPACKE_cppsvx CPPSVX;

     /* LAPACKE CPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpptrf) ( int matrix_layout , char uplo , 
                                                lapack_int n , lapack_complex_float * ap );

    Fptr_NL_LAPACKE_cpptrf CPPTRF;

    cppsvx_obj = new ppsvx_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    cppsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cppsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cppsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cppsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CPPSVX = (Fptr_NL_LAPACKE_cppsvx)dlsym(cppsvx_obj->hModule, "LAPACKE_cppsvx");
    ASSERT_TRUE(CPPSVX != NULL) << "failed to ppt the Netlib LAPACKE_cppsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pptrf API to compute the factorized A.  */
    if(cppsvx_obj->fact == 'F') {
        CPPTRF = (Fptr_NL_LAPACKE_cpptrf)dlsym(cppsvx_obj->hModule,"LAPACKE_cpptrf");
        ASSERT_TRUE(CPPTRF != NULL) << "failed to ppt the Netlib LAPACKE_cpptrf symbol";
            
        cppsvx_obj->inforef = CPPTRF( cppsvx_obj->matrix_layout,
                                      cppsvx_obj->uplo, cppsvx_obj->n,
                                      cppsvx_obj->afref);
                               
        cppsvx_obj->info = LAPACKE_cpptrf( cppsvx_obj->matrix_layout,
                                           cppsvx_obj->uplo, cppsvx_obj->n,
                                           cppsvx_obj->af);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cppsvx_obj->inforef = CPPSVX( cppsvx_obj->matrix_layout, cppsvx_obj->fact,
                                  cppsvx_obj->uplo, cppsvx_obj->n,
                                  cppsvx_obj->nrhs,
                                  cppsvx_obj->aref, 
                                  cppsvx_obj->afref,
                                  &cppsvx_obj->equedref,
                                  cppsvx_obj->sref,                               
                                  cppsvx_obj->bref, cppsvx_obj->ldb,
                                  cppsvx_obj->xref, cppsvx_obj->ldx,
                                  &cppsvx_obj->rcondref, 
                                  cppsvx_obj->ferrref,
                                  cppsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cppsvx_obj->info = LAPACKE_cppsvx( cppsvx_obj->matrix_layout, cppsvx_obj->fact,
                                  cppsvx_obj->uplo, cppsvx_obj->n,
                                  cppsvx_obj->nrhs,
                                  cppsvx_obj->a, 
                                  cppsvx_obj->af,
                                  &cppsvx_obj->equed,
                                  cppsvx_obj->s,                                  
                                  cppsvx_obj->b, cppsvx_obj->ldb,
                                  cppsvx_obj->x, cppsvx_obj->ldx,
                                  &cppsvx_obj->rcond, 
                                  cppsvx_obj->ferr,
                                  cppsvx_obj->berr);

    if( cppsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cppsvx is wrong\n", cppsvx_obj->info );
    }
    if( cppsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cppsvx is wrong\n", 
        cppsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cppsvx_obj->diff =  computeDiff_c( cppsvx_obj->b_bufsize, 
                cppsvx_obj->b, cppsvx_obj->bref );

    cppsvx_obj->diff_xerr =  computeDiff_c( cppsvx_obj->x_bufsize, 
                cppsvx_obj->x, cppsvx_obj->xref );

    cppsvx_obj->diff_berr =  computeDiff_s( cppsvx_obj->nrhs, 
                cppsvx_obj->berr, cppsvx_obj->berrref );
                
    cppsvx_obj->diff_ferr =  computeDiff_s( cppsvx_obj->nrhs, 
                cppsvx_obj->ferr, cppsvx_obj->ferrref );
    
    if( cppsvx_obj->fact != 'F') {
       cppsvx_obj->diff_s =  computeDiff_s( cppsvx_obj->n, 
                cppsvx_obj->s, cppsvx_obj->sref );
    }
}

TEST_F(cppsvx_test, cppsvx1) {
    EXPECT_NEAR(0.0, cppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cppsvx_obj->rcond - cppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(cppsvx_obj->equed, cppsvx_obj->equedref);
}

TEST_F(cppsvx_test, cppsvx2) {
    EXPECT_NEAR(0.0, cppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cppsvx_obj->rcond - cppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(cppsvx_obj->equed, cppsvx_obj->equedref);
}

TEST_F(cppsvx_test, cppsvx3) {
    EXPECT_NEAR(0.0, cppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cppsvx_obj->rcond - cppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(cppsvx_obj->equed, cppsvx_obj->equedref);
}

TEST_F(cppsvx_test, cppsvx4) {
    EXPECT_NEAR(0.0, cppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cppsvx_obj->rcond - cppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(cppsvx_obj->equed, cppsvx_obj->equedref);
}

/* Begin ppsvx_dcomplex_parameters  class definition */
class ppsvx_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr, diff_s;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_double *a, *aref; //The array ab contains the matrix A
      lapack_complex_double *af, *afref; //contains the factored form of the matrix A
      double* s, *sref; // the row scale factors for A.
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ppsvx_dcomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~ppsvx_dcomplex_parameters (); 
};  /* end of ppsvx_dcomplex_parameters  class definition */


/* Constructor ppsvx_dcomplex_parameters definition */
ppsvx_dcomplex_parameters:: ppsvx_dcomplex_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
    int a_buf_size =  n_i*(n_i+1)/2;    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
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
    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;
    diff_s = 0;
    diff = 0;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n ppsvx double:  n: %d, fact: %c uplo: %c   \
ldb: %d nrhs: %d   ldx: %d \n",  n, fact, uplo,  
                                          ldb, nrhs,  ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_buf_size);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &af, &afref, a_buf_size);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (s==NULL) || (sref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       ppsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, a_buf_size);
    memcpy(af, a, (a_buf_size*sizeof(lapack_complex_double)));
    memcpy(afref, a, (a_buf_size*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

ppsvx_dcomplex_parameters:: ~ppsvx_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppsvx_dcomplex_parameters object: destructor invoked. \n");
#endif
   ppsvx_free();
}

//  Test fixture class definition
class zppsvx_test  : public  ::testing::Test {
public:
   ppsvx_dcomplex_parameters  *zppsvx_obj;
   void SetUp();  
   void TearDown () { delete zppsvx_obj; }
};


void zppsvx_test::SetUp(){

    /* LAPACKE ZPPSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zppsvx) (int matrix_layout, char fact, 
                    char uplo, lapack_int n, lapack_int nrhs,
                    lapack_complex_double* ap, lapack_complex_double* afp,
                    char* equed, double* s, lapack_complex_double* b, 
                    lapack_int ldb, lapack_complex_double* x, lapack_int ldx, 
                    double* rcond, double* ferr, double* berr);

    Fptr_NL_LAPACKE_zppsvx ZPPSVX;

     /* LAPACKE ZPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpptrf) ( int matrix_layout , char uplo , 
                                                lapack_int n , lapack_complex_double * ap );

    Fptr_NL_LAPACKE_zpptrf ZPPTRF;

    zppsvx_obj = new ppsvx_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zppsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zppsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zppsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zppsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZPPSVX = (Fptr_NL_LAPACKE_zppsvx)dlsym(zppsvx_obj->hModule, "LAPACKE_zppsvx");
    ASSERT_TRUE(ZPPSVX != NULL) << "failed to ppt the Netlib LAPACKE_zppsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pptrf API to compute the factorized A.  */
    if(zppsvx_obj->fact == 'F') {
        ZPPTRF = (Fptr_NL_LAPACKE_zpptrf)dlsym(zppsvx_obj->hModule,"LAPACKE_zpptrf");
        ASSERT_TRUE(ZPPTRF != NULL) << "failed to ppt the Netlib LAPACKE_zpptrf symbol";
            
        zppsvx_obj->inforef = ZPPTRF( zppsvx_obj->matrix_layout,
                                      zppsvx_obj->uplo, zppsvx_obj->n,
                                      zppsvx_obj->afref);
                               
        zppsvx_obj->info = LAPACKE_zpptrf( zppsvx_obj->matrix_layout,
                                           zppsvx_obj->uplo, zppsvx_obj->n,
                                           zppsvx_obj->af);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zppsvx_obj->inforef = ZPPSVX( zppsvx_obj->matrix_layout, zppsvx_obj->fact,
                                  zppsvx_obj->uplo, zppsvx_obj->n,
                                  zppsvx_obj->nrhs,
                                  zppsvx_obj->aref, 
                                  zppsvx_obj->afref,
                                  &zppsvx_obj->equedref,
                                  zppsvx_obj->sref,                               
                                  zppsvx_obj->bref, zppsvx_obj->ldb,
                                  zppsvx_obj->xref, zppsvx_obj->ldx,
                                  &zppsvx_obj->rcondref, 
                                  zppsvx_obj->ferrref,
                                  zppsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zppsvx_obj->info = LAPACKE_zppsvx( zppsvx_obj->matrix_layout, zppsvx_obj->fact,
                                  zppsvx_obj->uplo, zppsvx_obj->n,
                                  zppsvx_obj->nrhs,
                                  zppsvx_obj->a, 
                                  zppsvx_obj->af,
                                  &zppsvx_obj->equed,
                                  zppsvx_obj->s,                                  
                                  zppsvx_obj->b, zppsvx_obj->ldb,
                                  zppsvx_obj->x, zppsvx_obj->ldx,
                                  &zppsvx_obj->rcond, 
                                  zppsvx_obj->ferr,
                                  zppsvx_obj->berr);

    if( zppsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zppsvx is wrong\n", zppsvx_obj->info );
    }
    if( zppsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zppsvx is wrong\n", 
        zppsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zppsvx_obj->diff =  computeDiff_z( zppsvx_obj->b_bufsize, 
                zppsvx_obj->b, zppsvx_obj->bref );

    zppsvx_obj->diff_xerr =  computeDiff_z( zppsvx_obj->x_bufsize, 
                zppsvx_obj->x, zppsvx_obj->xref );

    zppsvx_obj->diff_berr =  computeDiff_d( zppsvx_obj->nrhs, 
                zppsvx_obj->berr, zppsvx_obj->berrref );
                
    zppsvx_obj->diff_ferr =  computeDiff_d( zppsvx_obj->nrhs, 
                zppsvx_obj->ferr, zppsvx_obj->ferrref );
    
    if( zppsvx_obj->fact != 'F') {
       zppsvx_obj->diff_s =  computeDiff_d( zppsvx_obj->n, 
                zppsvx_obj->s, zppsvx_obj->sref );
    }
}

TEST_F(zppsvx_test, zppsvx1) {
    EXPECT_NEAR(0.0, zppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zppsvx_obj->rcond - zppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(zppsvx_obj->equed, zppsvx_obj->equedref);
}

TEST_F(zppsvx_test, zppsvx2) {
    EXPECT_NEAR(0.0, zppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zppsvx_obj->rcond - zppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(zppsvx_obj->equed, zppsvx_obj->equedref);
}

TEST_F(zppsvx_test, zppsvx3) {
    EXPECT_NEAR(0.0, zppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zppsvx_obj->rcond - zppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(zppsvx_obj->equed, zppsvx_obj->equedref);
}

TEST_F(zppsvx_test, zppsvx4) {
    EXPECT_NEAR(0.0, zppsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zppsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zppsvx_obj->rcond - zppsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
    EXPECT_EQ(zppsvx_obj->equed, zppsvx_obj->equedref);
}
