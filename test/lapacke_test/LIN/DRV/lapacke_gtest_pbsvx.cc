#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"
#define pbsvx_free() \
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

/* Begin pbsvx_float_parameters  class definition */
class pbsvx_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr, diff_s;
      float diff_a, diff_af;
      
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
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
      pbsvx_float_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i, lapack_int kd_i);
              
      ~pbsvx_float_parameters (); 
};  /* end of pbsvx_float_parameters  class definition */


/* Constructor pbsvx_float_parameters definition */
pbsvx_float_parameters:: pbsvx_float_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    kd = kd_i;
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
   printf(" \n pbsvx float:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &af, &afref, (n*n));
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
       pbsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand ( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(float)));
    memcpy(afref, a, (n*n*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_float_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

pbsvx_float_parameters:: ~pbsvx_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbsvx_float_parameters object: destructor invoked. \n");
#endif
   pbsvx_free();
}

//  Test fixture class definition
class spbsvx_test  : public  ::testing::Test {
public:
   pbsvx_float_parameters  *spbsvx_obj;
   void SetUp();  
   void TearDown () { delete spbsvx_obj; }
};


void spbsvx_test::SetUp(){

    /* LAPACKE SPBSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_spbsvx) (int matrix_layout,char fact, char uplo,
             lapack_int n, lapack_int kd, lapack_int nrhs, float* a, lapack_int lda, float* af,
             lapack_int ldaf, char* equed, float* s, float* b, lapack_int ldb,
             float* x, lapack_int ldx, float* rcond, float* ferr, float* berr);

    Fptr_NL_LAPACKE_spbsvx SPBSVX;

     /* LAPACKE SPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spbtrf) ( int matrix_layout , char uplo ,
                                lapack_int n , lapack_int kd, float * a , lapack_int lda );

    Fptr_NL_LAPACKE_spbtrf SPBTRF;

    spbsvx_obj = new pbsvx_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd);

    idx = Circular_Increment_Index(idx);

    spbsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spbsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spbsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spbsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SPBSVX = (Fptr_NL_LAPACKE_spbsvx)dlsym(spbsvx_obj->hModule, "LAPACKE_spbsvx");
    ASSERT_TRUE(SPBSVX != NULL) << "failed to pbt the Netlib LAPACKE_spbsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pbtrf API to compute the factorized A.  */
    if(spbsvx_obj->fact == 'F') {
        SPBTRF = (Fptr_NL_LAPACKE_spbtrf)dlsym(spbsvx_obj->hModule,"LAPACKE_spbtrf");
        ASSERT_TRUE(SPBTRF != NULL) << "failed to pbt the Netlib LAPACKE_spbtrf symbol";
            
        spbsvx_obj->inforef = SPBTRF( spbsvx_obj->matrix_layout,
                                      spbsvx_obj->uplo, spbsvx_obj->n,
                                      spbsvx_obj->kd,
                                      spbsvx_obj->afref,
                                      spbsvx_obj->lda);
                               
        spbsvx_obj->info = LAPACKE_spbtrf( spbsvx_obj->matrix_layout,
                                           spbsvx_obj->uplo, spbsvx_obj->n,
                                           spbsvx_obj->kd,
                                           spbsvx_obj->af,
                                           spbsvx_obj->lda);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    spbsvx_obj->inforef = SPBSVX( spbsvx_obj->matrix_layout, spbsvx_obj->fact,
                                  spbsvx_obj->uplo, spbsvx_obj->n,
                                  spbsvx_obj->kd,
                                  spbsvx_obj->nrhs,
                                  spbsvx_obj->aref, spbsvx_obj->lda, 
                                  spbsvx_obj->afref, spbsvx_obj->ldaf,
                                  &spbsvx_obj->equedref,
                                  spbsvx_obj->sref,                               
                                  spbsvx_obj->bref, spbsvx_obj->ldb,
                                  spbsvx_obj->xref, spbsvx_obj->ldx,
                                  &spbsvx_obj->rcondref, 
                                  spbsvx_obj->ferrref,
                                  spbsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    spbsvx_obj->info = LAPACKE_spbsvx( spbsvx_obj->matrix_layout, spbsvx_obj->fact,
                                  spbsvx_obj->uplo, spbsvx_obj->n,
                                  spbsvx_obj->kd,
                                  spbsvx_obj->nrhs,
                                  spbsvx_obj->a, spbsvx_obj->lda, 
                                  spbsvx_obj->af, spbsvx_obj->ldaf,
                                  &spbsvx_obj->equed,
                                  spbsvx_obj->s,                                  
                                  spbsvx_obj->b, spbsvx_obj->ldb,
                                  spbsvx_obj->x, spbsvx_obj->ldx,
                                  &spbsvx_obj->rcond, 
                                  spbsvx_obj->ferr,
                                  spbsvx_obj->berr);

    if( spbsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_spbsvx is wrong\n", spbsvx_obj->info );
    }
    if( spbsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spbsvx is wrong\n", 
        spbsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    spbsvx_obj->diff =  computeDiff_s( spbsvx_obj->b_bufsize, 
                spbsvx_obj->b, spbsvx_obj->bref );

    spbsvx_obj->diff =  computeDiff_s( spbsvx_obj->n * spbsvx_obj->n, 
                spbsvx_obj->a, spbsvx_obj->aref );

    spbsvx_obj->diff =  computeDiff_s( spbsvx_obj->n * spbsvx_obj->n, 
                spbsvx_obj->af, spbsvx_obj->afref );

    spbsvx_obj->diff_xerr =  computeDiff_s( spbsvx_obj->x_bufsize, 
                spbsvx_obj->x, spbsvx_obj->xref );

    spbsvx_obj->diff_berr =  computeDiff_s( spbsvx_obj->nrhs, 
                spbsvx_obj->berr, spbsvx_obj->berrref );
                
    spbsvx_obj->diff_ferr =  computeDiff_s( spbsvx_obj->nrhs, 
                spbsvx_obj->ferr, spbsvx_obj->ferrref );

    if( spbsvx_obj->fact != 'F') {
       spbsvx_obj->diff_s =  computeDiff_s( spbsvx_obj->n, 
                spbsvx_obj->s, spbsvx_obj->sref );
    }

}

TEST_F(spbsvx_test, spbsvx1) {
    EXPECT_NEAR(0.0, spbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (spbsvx_obj->rcond - spbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbsvx_test, spbsvx2) {
    EXPECT_NEAR(0.0, spbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (spbsvx_obj->rcond - spbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbsvx_test, spbsvx3) {
    EXPECT_NEAR(0.0, spbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (spbsvx_obj->rcond - spbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbsvx_test, spbsvx4) {
    EXPECT_NEAR(0.0, spbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, spbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (spbsvx_obj->rcond - spbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}


/* Begin pbsvx_double_parameters  class definition */
class pbsvx_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr, diff_s;
      double diff_a, diff_af;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
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
      pbsvx_double_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i, lapack_int kd_i);
              
      ~pbsvx_double_parameters (); 
};  /* end of pbsvx_double_parameters  class definition */


/* Constructor pbsvx_double_parameters definition */
pbsvx_double_parameters:: pbsvx_double_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    kd = kd_i;
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
   printf(" \n pbsvx Double:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &af, &afref, (n*n));
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
       pbsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand ( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(double)));
    memcpy(afref, a, (n*n*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_double_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

pbsvx_double_parameters:: ~pbsvx_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbsvx_double_parameters object: destructor invoked. \n");
#endif
   pbsvx_free();
}

//  Test fixture class definition
class dpbsvx_test  : public  ::testing::Test {
public:
   pbsvx_double_parameters  *dpbsvx_obj;
   void SetUp();  
   void TearDown () { delete dpbsvx_obj; }
};


void dpbsvx_test::SetUp(){

    /* LAPACKE DPBSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_dpbsvx) (int matrix_layout,char fact, char uplo,
             lapack_int n, lapack_int kd, lapack_int nrhs, double* a, lapack_int lda, double* af,
             lapack_int ldaf, char* equed, double* s, double* b, lapack_int ldb,
             double* x, lapack_int ldx, double* rcond, double* ferr, double* berr);

    Fptr_NL_LAPACKE_dpbsvx DPBSVX;

     /* LAPACKE DPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpbtrf) ( int matrix_layout , char uplo ,
                                lapack_int n , lapack_int kd, double * a , lapack_int lda );

    Fptr_NL_LAPACKE_dpbtrf DPBTRF;

    dpbsvx_obj = new pbsvx_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd);

    idx = Circular_Increment_Index(idx);

    dpbsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpbsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpbsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpbsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DPBSVX = (Fptr_NL_LAPACKE_dpbsvx)dlsym(dpbsvx_obj->hModule, "LAPACKE_dpbsvx");
    ASSERT_TRUE(DPBSVX != NULL) << "failed to pbt the Netlib LAPACKE_dpbsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pbtrf API to compute the factorized A.  */
    if(dpbsvx_obj->fact == 'F') {
        DPBTRF = (Fptr_NL_LAPACKE_dpbtrf)dlsym(dpbsvx_obj->hModule,"LAPACKE_dpbtrf");
        ASSERT_TRUE(DPBTRF != NULL) << "failed to pbt the Netlib LAPACKE_dpbtrf symbol";
            
        dpbsvx_obj->inforef = DPBTRF( dpbsvx_obj->matrix_layout,
                                      dpbsvx_obj->uplo, dpbsvx_obj->n,
                                      dpbsvx_obj->kd,
                                      dpbsvx_obj->afref,
                                      dpbsvx_obj->lda);
                               
        dpbsvx_obj->info = LAPACKE_dpbtrf( dpbsvx_obj->matrix_layout,
                                           dpbsvx_obj->uplo, dpbsvx_obj->n,
                                           dpbsvx_obj->kd,
                                           dpbsvx_obj->af,
                                           dpbsvx_obj->lda);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dpbsvx_obj->inforef = DPBSVX( dpbsvx_obj->matrix_layout, dpbsvx_obj->fact,
                                  dpbsvx_obj->uplo, dpbsvx_obj->n,
                                  dpbsvx_obj->kd,
                                  dpbsvx_obj->nrhs,
                                  dpbsvx_obj->aref, dpbsvx_obj->lda, 
                                  dpbsvx_obj->afref, dpbsvx_obj->ldaf,
                                  &dpbsvx_obj->equedref,
                                  dpbsvx_obj->sref,                               
                                  dpbsvx_obj->bref, dpbsvx_obj->ldb,
                                  dpbsvx_obj->xref, dpbsvx_obj->ldx,
                                  &dpbsvx_obj->rcondref, 
                                  dpbsvx_obj->ferrref,
                                  dpbsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dpbsvx_obj->info = LAPACKE_dpbsvx( dpbsvx_obj->matrix_layout, dpbsvx_obj->fact,
                                  dpbsvx_obj->uplo, dpbsvx_obj->n,
                                  dpbsvx_obj->kd,
                                  dpbsvx_obj->nrhs,
                                  dpbsvx_obj->a, dpbsvx_obj->lda, 
                                  dpbsvx_obj->af, dpbsvx_obj->ldaf,
                                  &dpbsvx_obj->equed,
                                  dpbsvx_obj->s,                                  
                                  dpbsvx_obj->b, dpbsvx_obj->ldb,
                                  dpbsvx_obj->x, dpbsvx_obj->ldx,
                                  &dpbsvx_obj->rcond, 
                                  dpbsvx_obj->ferr,
                                  dpbsvx_obj->berr);

    if( dpbsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dpbsvx is wrong\n", dpbsvx_obj->info );
    }
    if( dpbsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpbsvx is wrong\n", 
        dpbsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dpbsvx_obj->diff =  computeDiff_d( dpbsvx_obj->b_bufsize, 
                dpbsvx_obj->b, dpbsvx_obj->bref );

    dpbsvx_obj->diff_a =  computeDiff_d( dpbsvx_obj->n * dpbsvx_obj->n, 
                dpbsvx_obj->a, dpbsvx_obj->aref );

    dpbsvx_obj->diff_af =  computeDiff_d( dpbsvx_obj->n * dpbsvx_obj->n, 
                dpbsvx_obj->af, dpbsvx_obj->afref );

    dpbsvx_obj->diff_xerr =  computeDiff_d( dpbsvx_obj->x_bufsize, 
                dpbsvx_obj->x, dpbsvx_obj->xref );

    dpbsvx_obj->diff_berr =  computeDiff_d( dpbsvx_obj->nrhs, 
                dpbsvx_obj->berr, dpbsvx_obj->berrref );
                
    dpbsvx_obj->diff_ferr =  computeDiff_d( dpbsvx_obj->nrhs, 
                dpbsvx_obj->ferr, dpbsvx_obj->ferrref );

    if( dpbsvx_obj->fact != 'F') {
       dpbsvx_obj->diff_s =  computeDiff_d( dpbsvx_obj->n, 
                dpbsvx_obj->s, dpbsvx_obj->sref );
    }
}

TEST_F(dpbsvx_test, dpbsvx1) {
    EXPECT_NEAR(0.0, dpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dpbsvx_obj->rcond - dpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbsvx_test, dpbsvx2) {
    EXPECT_NEAR(0.0, dpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dpbsvx_obj->rcond - dpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbsvx_test, dpbsvx3) {
    EXPECT_NEAR(0.0, dpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dpbsvx_obj->rcond - dpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbsvx_test, dpbsvx4) {
    EXPECT_NEAR(0.0, dpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dpbsvx_obj->rcond - dpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin pbsvx_scomplex_parameters  class definition */
class pbsvx_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr, diff_s;
      float diff_a, diff_af;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
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
      pbsvx_scomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i, lapack_int kd_i);
              
      ~pbsvx_scomplex_parameters (); 
};  /* end of pbsvx_scomplex_parameters  class definition */


/* Constructor pbsvx_scomplex_parameters definition */
pbsvx_scomplex_parameters:: pbsvx_scomplex_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    kd = kd_i;
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
   printf(" \n pbsvx scomplex:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &af, &afref, (n*n));
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
       pbsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    //lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    lapacke_gtest_init_scomplex_buffer_pair_rand(a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(lapack_complex_float)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

pbsvx_scomplex_parameters:: ~pbsvx_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbsvx_scomplex_parameters object: destructor invoked. \n");
#endif
   pbsvx_free();
}

//  Test fixture class definition
class cpbsvx_test  : public  ::testing::Test {
public:
   pbsvx_scomplex_parameters  *cpbsvx_obj;
   void SetUp();  
   void TearDown () { delete cpbsvx_obj; }
};


void cpbsvx_test::SetUp(){

    /* LAPACKE CPBSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_cpbsvx) (int matrix_layout,char fact, char uplo,
             lapack_int n, lapack_int kd, lapack_int nrhs, lapack_complex_float* a, lapack_int lda, lapack_complex_float* af,
             lapack_int ldaf, char* equed, float* s, lapack_complex_float* b, lapack_int ldb,
             lapack_complex_float* x, lapack_int ldx, float* rcond, float* ferr, float* berr);

    Fptr_NL_LAPACKE_cpbsvx CPBSVX;

     /* LAPACKE CPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpbtrf) ( int matrix_layout , char uplo ,
                                lapack_int n , lapack_int kd, lapack_complex_float * a , lapack_int lda );

    Fptr_NL_LAPACKE_cpbtrf CPBTRF;

    cpbsvx_obj = new pbsvx_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd);

    idx = Circular_Increment_Index(idx);

    cpbsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpbsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpbsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpbsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CPBSVX = (Fptr_NL_LAPACKE_cpbsvx)dlsym(cpbsvx_obj->hModule, "LAPACKE_cpbsvx");
    ASSERT_TRUE(CPBSVX != NULL) << "failed to pbt the Netlib LAPACKE_cpbsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pbtrf API to compute the factorized A.  */
    if(cpbsvx_obj->fact == 'F') {
        CPBTRF = (Fptr_NL_LAPACKE_cpbtrf)dlsym(cpbsvx_obj->hModule,"LAPACKE_cpbtrf");
        ASSERT_TRUE(CPBTRF != NULL) << "failed to pbt the Netlib LAPACKE_cpbtrf symbol";
            
        cpbsvx_obj->inforef = CPBTRF( cpbsvx_obj->matrix_layout,
                                      cpbsvx_obj->uplo, cpbsvx_obj->n,
                                      cpbsvx_obj->kd,
                                      cpbsvx_obj->afref,
                                      cpbsvx_obj->lda);
                               
        cpbsvx_obj->info = LAPACKE_cpbtrf( cpbsvx_obj->matrix_layout,
                                           cpbsvx_obj->uplo, cpbsvx_obj->n,
                                           cpbsvx_obj->kd,
                                           cpbsvx_obj->af,
                                           cpbsvx_obj->lda);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cpbsvx_obj->inforef = CPBSVX( cpbsvx_obj->matrix_layout, cpbsvx_obj->fact,
                                  cpbsvx_obj->uplo, cpbsvx_obj->n,
                                  cpbsvx_obj->kd,
                                  cpbsvx_obj->nrhs,
                                  cpbsvx_obj->aref, cpbsvx_obj->lda, 
                                  cpbsvx_obj->afref, cpbsvx_obj->ldaf,
                                  &cpbsvx_obj->equedref,
                                  cpbsvx_obj->sref,                               
                                  cpbsvx_obj->bref, cpbsvx_obj->ldb,
                                  cpbsvx_obj->xref, cpbsvx_obj->ldx,
                                  &cpbsvx_obj->rcondref, 
                                  cpbsvx_obj->ferrref,
                                  cpbsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cpbsvx_obj->info = LAPACKE_cpbsvx( cpbsvx_obj->matrix_layout, cpbsvx_obj->fact,
                                  cpbsvx_obj->uplo, cpbsvx_obj->n,
                                  cpbsvx_obj->kd,
                                  cpbsvx_obj->nrhs,
                                  cpbsvx_obj->a, cpbsvx_obj->lda, 
                                  cpbsvx_obj->af, cpbsvx_obj->ldaf,
                                  &cpbsvx_obj->equed,
                                  cpbsvx_obj->s,                                  
                                  cpbsvx_obj->b, cpbsvx_obj->ldb,
                                  cpbsvx_obj->x, cpbsvx_obj->ldx,
                                  &cpbsvx_obj->rcond, 
                                  cpbsvx_obj->ferr,
                                  cpbsvx_obj->berr);

    if( cpbsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cpbsvx is wrong\n", cpbsvx_obj->info );
    }
    if( cpbsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpbsvx is wrong\n", 
        cpbsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cpbsvx_obj->diff =  computeDiff_c( cpbsvx_obj->b_bufsize, 
                cpbsvx_obj->b, cpbsvx_obj->bref );

    cpbsvx_obj->diff_a =  computeDiff_c( cpbsvx_obj->n * cpbsvx_obj->n, 
                cpbsvx_obj->a, cpbsvx_obj->aref );

    cpbsvx_obj->diff_af =  computeDiff_c( cpbsvx_obj->n * cpbsvx_obj->n, 
                cpbsvx_obj->af, cpbsvx_obj->afref );

    cpbsvx_obj->diff_xerr =  computeDiff_c( cpbsvx_obj->x_bufsize, 
                cpbsvx_obj->x, cpbsvx_obj->xref );

    cpbsvx_obj->diff_berr =  computeDiff_s( cpbsvx_obj->nrhs, 
                cpbsvx_obj->berr, cpbsvx_obj->berrref );
                
    cpbsvx_obj->diff_ferr =  computeDiff_s( cpbsvx_obj->nrhs, 
                cpbsvx_obj->ferr, cpbsvx_obj->ferrref );

    if( cpbsvx_obj->fact != 'F') {
       cpbsvx_obj->diff_s =  computeDiff_s( cpbsvx_obj->n, 
                cpbsvx_obj->s, cpbsvx_obj->sref );
    }
}

TEST_F(cpbsvx_test, cpbsvx1) {
    EXPECT_NEAR(0.0, cpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cpbsvx_obj->rcond - cpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbsvx_test, cpbsvx2) {
    EXPECT_NEAR(0.0, cpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cpbsvx_obj->rcond - cpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbsvx_test, cpbsvx3) {
    EXPECT_NEAR(0.0, cpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cpbsvx_obj->rcond - cpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbsvx_test, cpbsvx4) {
    EXPECT_NEAR(0.0, cpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cpbsvx_obj->rcond - cpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin pbsvx_dcomplex_parameters  class definition */
class pbsvx_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr, diff_s;
      double diff_a, diff_af;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char fact; //  Must be 'F', 'N', or 'E'.
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
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
      pbsvx_dcomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i, lapack_int kd_i);
              
      ~pbsvx_dcomplex_parameters (); 
};  /* end of pbsvx_dcomplex_parameters  class definition */


/* Constructor pbsvx_dcomplex_parameters definition */
pbsvx_dcomplex_parameters:: pbsvx_dcomplex_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
    n = n_i;
    nrhs = nrhs_i;
    kd = kd_i;
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
   printf(" \n pbsvx dcomplex:  n: %d, fact: %c uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, fact, uplo, lda, 
                                          ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &af, &afref, (n*n));
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
       pbsvx_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    //lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    lapacke_gtest_init_dcomplex_buffer_pair_rand(a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(lapack_complex_double)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

pbsvx_dcomplex_parameters:: ~pbsvx_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbsvx_dcomplex_parameters object: destructor invoked. \n");
#endif
   pbsvx_free();
}

//  Test fixture class definition
class zpbsvx_test  : public  ::testing::Test {
public:
   pbsvx_dcomplex_parameters  *zpbsvx_obj;
   void SetUp();  
   void TearDown () { delete zpbsvx_obj; }
};


void zpbsvx_test::SetUp(){

    /* LAPACKE ZPBSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zpbsvx) (int matrix_layout,char fact, char uplo,
             lapack_int n, lapack_int kd, lapack_int nrhs, lapack_complex_double* a, lapack_int lda, lapack_complex_double* af,
             lapack_int ldaf, char* equed, double* s, lapack_complex_double* b, lapack_int ldb,
             lapack_complex_double* x, lapack_int ldx, double* rcond, double* ferr, double* berr);

    Fptr_NL_LAPACKE_zpbsvx ZPBSVX;

     /* LAPACKE ZPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpbtrf) ( int matrix_layout , char uplo ,
                                lapack_int n , lapack_int kd, lapack_complex_double * a , lapack_int lda );

    Fptr_NL_LAPACKE_zpbtrf ZPBTRF;

    zpbsvx_obj = new pbsvx_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd);

    idx = Circular_Increment_Index(idx);

    zpbsvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpbsvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpbsvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpbsvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZPBSVX = (Fptr_NL_LAPACKE_zpbsvx)dlsym(zpbsvx_obj->hModule, "LAPACKE_zpbsvx");
    ASSERT_TRUE(ZPBSVX != NULL) << "failed to pbt the Netlib LAPACKE_zpbsvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the pbtrf API to compute the factorized A. */ 
    if(zpbsvx_obj->fact == 'F') {
        ZPBTRF = (Fptr_NL_LAPACKE_zpbtrf)dlsym(zpbsvx_obj->hModule,"LAPACKE_zpbtrf");
        ASSERT_TRUE(ZPBTRF != NULL) << "failed to pbt the Netlib LAPACKE_zpbtrf symbol";
            
        zpbsvx_obj->inforef = ZPBTRF( zpbsvx_obj->matrix_layout,
                                      zpbsvx_obj->uplo, zpbsvx_obj->n,
                                      zpbsvx_obj->kd,
                                      zpbsvx_obj->afref,
                                      zpbsvx_obj->lda);
                               
        zpbsvx_obj->info = LAPACKE_zpbtrf( zpbsvx_obj->matrix_layout,
                                           zpbsvx_obj->uplo, zpbsvx_obj->n,
                                           zpbsvx_obj->kd,
                                           zpbsvx_obj->af,
                                           zpbsvx_obj->lda);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zpbsvx_obj->inforef = ZPBSVX( zpbsvx_obj->matrix_layout, zpbsvx_obj->fact,
                                  zpbsvx_obj->uplo, zpbsvx_obj->n,
                                  zpbsvx_obj->kd,
                                  zpbsvx_obj->nrhs,
                                  zpbsvx_obj->aref, zpbsvx_obj->lda, 
                                  zpbsvx_obj->afref, zpbsvx_obj->ldaf,
                                  &zpbsvx_obj->equedref,
                                  zpbsvx_obj->sref,                               
                                  zpbsvx_obj->bref, zpbsvx_obj->ldb,
                                  zpbsvx_obj->xref, zpbsvx_obj->ldx,
                                  &zpbsvx_obj->rcondref, 
                                  zpbsvx_obj->ferrref,
                                  zpbsvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zpbsvx_obj->info = LAPACKE_zpbsvx( zpbsvx_obj->matrix_layout, zpbsvx_obj->fact,
                                  zpbsvx_obj->uplo, zpbsvx_obj->n,
                                  zpbsvx_obj->kd,
                                  zpbsvx_obj->nrhs,
                                  zpbsvx_obj->a, zpbsvx_obj->lda, 
                                  zpbsvx_obj->af, zpbsvx_obj->ldaf,
                                  &zpbsvx_obj->equed,
                                  zpbsvx_obj->s,                                  
                                  zpbsvx_obj->b, zpbsvx_obj->ldb,
                                  zpbsvx_obj->x, zpbsvx_obj->ldx,
                                  &zpbsvx_obj->rcond, 
                                  zpbsvx_obj->ferr,
                                  zpbsvx_obj->berr);

    if( zpbsvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zpbsvx is wrong\n", zpbsvx_obj->info );
    }
    if( zpbsvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpbsvx is wrong\n", 
        zpbsvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zpbsvx_obj->diff =  computeDiff_z( zpbsvx_obj->b_bufsize, 
                zpbsvx_obj->b, zpbsvx_obj->bref );

    zpbsvx_obj->diff_a =  computeDiff_z( zpbsvx_obj->n *zpbsvx_obj->n, 
                zpbsvx_obj->a, zpbsvx_obj->aref );

    zpbsvx_obj->diff_af =  computeDiff_z( zpbsvx_obj->n *zpbsvx_obj->n, 
                zpbsvx_obj->af, zpbsvx_obj->afref );

    zpbsvx_obj->diff_xerr =  computeDiff_z( zpbsvx_obj->x_bufsize, 
                zpbsvx_obj->x, zpbsvx_obj->xref );

    zpbsvx_obj->diff_berr =  computeDiff_d( zpbsvx_obj->nrhs, 
                zpbsvx_obj->berr, zpbsvx_obj->berrref );
                
    zpbsvx_obj->diff_ferr =  computeDiff_d( zpbsvx_obj->nrhs, 
                zpbsvx_obj->ferr, zpbsvx_obj->ferrref );
    
    if( zpbsvx_obj->fact != 'F') {
       zpbsvx_obj->diff_s =  computeDiff_d( zpbsvx_obj->n, 
                zpbsvx_obj->s, zpbsvx_obj->sref );
    }
}

TEST_F(zpbsvx_test, zpbsvx1) {
    EXPECT_NEAR(0.0, zpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zpbsvx_obj->rcond - zpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbsvx_test, zpbsvx2) {
    EXPECT_NEAR(0.0, zpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zpbsvx_obj->rcond - zpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbsvx_test, zpbsvx3) {
    EXPECT_NEAR(0.0, zpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zpbsvx_obj->rcond - zpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbsvx_test, zpbsvx4) {
    EXPECT_NEAR(0.0, zpbsvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_af, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zpbsvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zpbsvx_obj->rcond - zpbsvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}