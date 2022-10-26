#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define posvx_free() \
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

/* Begin posvx_float_parameters  class definition */
class posvx_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

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
      posvx_float_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~posvx_float_parameters (); 
};  /* end of posvx_float_parameters  class definition */


/* Constructor posvx_float_parameters definition */
posvx_float_parameters:: posvx_float_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
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
   printf(" \n posvx float:  n: %d, fact: %c uplo: %c  lda: %d  \
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
       posvx_free();
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
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

posvx_float_parameters:: ~posvx_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" posvx_float_parameters object: destructor invoked. \n");
#endif
   posvx_free();
}

//  Test fixture class definition
class sposvx_test  : public  ::testing::Test {
public:
   posvx_float_parameters  *sposvx_obj;
   void SetUp();  
   void TearDown () { delete sposvx_obj; }
};


void sposvx_test::SetUp(){

    /* LAPACKE SPOSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_sposvx) (int matrix_layout,char fact, char uplo,
             lapack_int n, lapack_int nrhs, float* a, lapack_int lda, float* af,
             lapack_int ldaf, char* equed, float* s, float* b, lapack_int ldb,
             float* x, lapack_int ldx, float* rcond, float* ferr, float* berr);

    Fptr_NL_LAPACKE_sposvx SPOSVX;

     /* LAPACKE SPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spotrf) ( int matrix_layout , char uplo ,
                                lapack_int n , float * a , lapack_int lda );

    Fptr_NL_LAPACKE_spotrf SPOTRF;

    sposvx_obj = new posvx_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    sposvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sposvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sposvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sposvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SPOSVX = (Fptr_NL_LAPACKE_sposvx)dlsym(sposvx_obj->hModule, "LAPACKE_sposvx");
    ASSERT_TRUE(SPOSVX != NULL) << "failed to pot the Netlib LAPACKE_sposvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the potrf API to compute the factorized A.  */
    if(sposvx_obj->fact == 'F') {
        SPOTRF = (Fptr_NL_LAPACKE_spotrf)dlsym(sposvx_obj->hModule,"LAPACKE_spotrf");
        ASSERT_TRUE(SPOTRF != NULL) << "failed to pot the Netlib LAPACKE_spotrf symbol";
            
        sposvx_obj->inforef = SPOTRF( sposvx_obj->matrix_layout,
                                      sposvx_obj->uplo, sposvx_obj->n,
                                      sposvx_obj->afref,
                                      sposvx_obj->lda);
                               
        sposvx_obj->info = LAPACKE_spotrf( sposvx_obj->matrix_layout,
                                           sposvx_obj->uplo, sposvx_obj->n,
                                           sposvx_obj->af,
                                           sposvx_obj->lda);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    sposvx_obj->inforef = SPOSVX( sposvx_obj->matrix_layout, sposvx_obj->fact,
                                  sposvx_obj->uplo, sposvx_obj->n,
                                  sposvx_obj->nrhs,
                                  sposvx_obj->aref, sposvx_obj->lda, 
                                  sposvx_obj->afref, sposvx_obj->ldaf,
                                  &sposvx_obj->equedref,
                                  sposvx_obj->sref,                               
                                  sposvx_obj->bref, sposvx_obj->ldb,
                                  sposvx_obj->xref, sposvx_obj->ldx,
                                  &sposvx_obj->rcondref, 
                                  sposvx_obj->ferrref,
                                  sposvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    sposvx_obj->info = LAPACKE_sposvx( sposvx_obj->matrix_layout, sposvx_obj->fact,
                                  sposvx_obj->uplo, sposvx_obj->n,
                                  sposvx_obj->nrhs,
                                  sposvx_obj->a, sposvx_obj->lda, 
                                  sposvx_obj->af, sposvx_obj->ldaf,
                                  &sposvx_obj->equed,
                                  sposvx_obj->s,                                  
                                  sposvx_obj->b, sposvx_obj->ldb,
                                  sposvx_obj->x, sposvx_obj->ldx,
                                  &sposvx_obj->rcond, 
                                  sposvx_obj->ferr,
                                  sposvx_obj->berr);

    if( sposvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sposvx is wrong\n", sposvx_obj->info );
    }
    if( sposvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sposvx is wrong\n", 
        sposvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sposvx_obj->diff =  computeDiff_s( sposvx_obj->b_bufsize, 
                sposvx_obj->b, sposvx_obj->bref );

    sposvx_obj->diff_xerr =  computeDiff_s( sposvx_obj->x_bufsize, 
                sposvx_obj->x, sposvx_obj->xref );

    sposvx_obj->diff_berr =  computeDiff_s( sposvx_obj->nrhs, 
                sposvx_obj->berr, sposvx_obj->berrref );
                
    sposvx_obj->diff_ferr =  computeDiff_s( sposvx_obj->nrhs, 
                sposvx_obj->ferr, sposvx_obj->ferrref );
}

TEST_F(sposvx_test, sposvx1) {
    EXPECT_NEAR(0.0, sposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sposvx_obj->rcond - sposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sposvx_test, sposvx2) {
    EXPECT_NEAR(0.0, sposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sposvx_obj->rcond - sposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sposvx_test, sposvx3) {
    EXPECT_NEAR(0.0, sposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sposvx_obj->rcond - sposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sposvx_test, sposvx4) {
    EXPECT_NEAR(0.0, sposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, sposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (sposvx_obj->rcond - sposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}


/* Begin posvx_double_parameters  class definition */
class posvx_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

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
      posvx_double_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~posvx_double_parameters (); 
};  /* end of posvx_double_parameters  class definition */


/* Constructor posvx_double_parameters definition */
posvx_double_parameters:: posvx_double_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
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
   printf(" \n posvx Double:  n: %d, fact: %c uplo: %c  lda: %d  \
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
       posvx_free();
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
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

posvx_double_parameters:: ~posvx_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" posvx_double_parameters object: destructor invoked. \n");
#endif
   posvx_free();
}

//  Test fixture class definition
class dposvx_test  : public  ::testing::Test {
public:
   posvx_double_parameters  *dposvx_obj;
   void SetUp();  
   void TearDown () { delete dposvx_obj; }
};


void dposvx_test::SetUp(){

    /* LAPACKE DPOSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_dposvx) (int matrix_layout,char fact, char uplo,
             lapack_int n, lapack_int nrhs, double* a, lapack_int lda, double* af,
             lapack_int ldaf, char* equed, double* s, double* b, lapack_int ldb,
             double* x, lapack_int ldx, double* rcond, double* ferr, double* berr);

    Fptr_NL_LAPACKE_dposvx DPOSVX;

     /* LAPACKE DPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpotrf) ( int matrix_layout , char uplo ,
                                lapack_int n , double * a , lapack_int lda );

    Fptr_NL_LAPACKE_dpotrf DPOTRF;

    dposvx_obj = new posvx_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    dposvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dposvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dposvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dposvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DPOSVX = (Fptr_NL_LAPACKE_dposvx)dlsym(dposvx_obj->hModule, "LAPACKE_dposvx");
    ASSERT_TRUE(DPOSVX != NULL) << "failed to pot the Netlib LAPACKE_dposvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the potrf API to compute the factorized A.  */
    if(dposvx_obj->fact == 'F') {
        DPOTRF = (Fptr_NL_LAPACKE_dpotrf)dlsym(dposvx_obj->hModule,"LAPACKE_dpotrf");
        ASSERT_TRUE(DPOTRF != NULL) << "failed to pot the Netlib LAPACKE_dpotrf symbol";
            
        dposvx_obj->inforef = DPOTRF( dposvx_obj->matrix_layout,
                                      dposvx_obj->uplo, dposvx_obj->n,
                                      dposvx_obj->afref,
                                      dposvx_obj->lda);
                               
        dposvx_obj->info = LAPACKE_dpotrf( dposvx_obj->matrix_layout,
                                           dposvx_obj->uplo, dposvx_obj->n,
                                           dposvx_obj->af,
                                           dposvx_obj->lda);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dposvx_obj->inforef = DPOSVX( dposvx_obj->matrix_layout, dposvx_obj->fact,
                                  dposvx_obj->uplo, dposvx_obj->n,
                                  dposvx_obj->nrhs,
                                  dposvx_obj->aref, dposvx_obj->lda, 
                                  dposvx_obj->afref, dposvx_obj->ldaf,
                                  &dposvx_obj->equedref,
                                  dposvx_obj->sref,                               
                                  dposvx_obj->bref, dposvx_obj->ldb,
                                  dposvx_obj->xref, dposvx_obj->ldx,
                                  &dposvx_obj->rcondref, 
                                  dposvx_obj->ferrref,
                                  dposvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dposvx_obj->info = LAPACKE_dposvx( dposvx_obj->matrix_layout, dposvx_obj->fact,
                                  dposvx_obj->uplo, dposvx_obj->n,
                                  dposvx_obj->nrhs,
                                  dposvx_obj->a, dposvx_obj->lda, 
                                  dposvx_obj->af, dposvx_obj->ldaf,
                                  &dposvx_obj->equed,
                                  dposvx_obj->s,                                  
                                  dposvx_obj->b, dposvx_obj->ldb,
                                  dposvx_obj->x, dposvx_obj->ldx,
                                  &dposvx_obj->rcond, 
                                  dposvx_obj->ferr,
                                  dposvx_obj->berr);

    if( dposvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dposvx is wrong\n", dposvx_obj->info );
    }
    if( dposvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dposvx is wrong\n", 
        dposvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dposvx_obj->diff =  computeDiff_d( dposvx_obj->b_bufsize, 
                dposvx_obj->b, dposvx_obj->bref );

    dposvx_obj->diff_xerr =  computeDiff_d( dposvx_obj->x_bufsize, 
                dposvx_obj->x, dposvx_obj->xref );

    dposvx_obj->diff_berr =  computeDiff_d( dposvx_obj->nrhs, 
                dposvx_obj->berr, dposvx_obj->berrref );
                
    dposvx_obj->diff_ferr =  computeDiff_d( dposvx_obj->nrhs, 
                dposvx_obj->ferr, dposvx_obj->ferrref );
}

TEST_F(dposvx_test, dposvx1) {
    EXPECT_NEAR(0.0, dposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dposvx_obj->rcond - dposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dposvx_test, dposvx2) {
    EXPECT_NEAR(0.0, dposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dposvx_obj->rcond - dposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dposvx_test, dposvx3) {
    EXPECT_NEAR(0.0, dposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dposvx_obj->rcond - dposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dposvx_test, dposvx4) {
    EXPECT_NEAR(0.0, dposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, dposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (dposvx_obj->rcond - dposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin posvx_scomplex_parameters  class definition */
class posvx_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

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
      posvx_scomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~posvx_scomplex_parameters (); 
};  /* end of posvx_scomplex_parameters  class definition */


/* Constructor posvx_scomplex_parameters definition */
posvx_scomplex_parameters:: posvx_scomplex_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
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
   printf(" \n posvx scomplex:  n: %d, fact: %c uplo: %c  lda: %d  \
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
       posvx_free();
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
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

posvx_scomplex_parameters:: ~posvx_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" posvx_scomplex_parameters object: destructor invoked. \n");
#endif
   posvx_free();
}

//  Test fixture class definition
class cposvx_test  : public  ::testing::Test {
public:
   posvx_scomplex_parameters  *cposvx_obj;
   void SetUp();  
   void TearDown () { delete cposvx_obj; }
};


void cposvx_test::SetUp(){

    /* LAPACKE CPOSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_cposvx) (int matrix_layout,char fact, char uplo,
             lapack_int n, lapack_int nrhs, lapack_complex_float* a, lapack_int lda, lapack_complex_float* af,
             lapack_int ldaf, char* equed, float* s, lapack_complex_float* b, lapack_int ldb,
             lapack_complex_float* x, lapack_int ldx, float* rcond, float* ferr, float* berr);

    Fptr_NL_LAPACKE_cposvx CPOSVX;

     /* LAPACKE CPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpotrf) ( int matrix_layout , char uplo ,
                                lapack_int n , lapack_complex_float * a , lapack_int lda );

    Fptr_NL_LAPACKE_cpotrf CPOTRF;

    cposvx_obj = new posvx_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    cposvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cposvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cposvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cposvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CPOSVX = (Fptr_NL_LAPACKE_cposvx)dlsym(cposvx_obj->hModule, "LAPACKE_cposvx");
    ASSERT_TRUE(CPOSVX != NULL) << "failed to pot the Netlib LAPACKE_cposvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the potrf API to compute the factorized A.  */
    if(cposvx_obj->fact == 'F') {
        CPOTRF = (Fptr_NL_LAPACKE_cpotrf)dlsym(cposvx_obj->hModule,"LAPACKE_cpotrf");
        ASSERT_TRUE(CPOTRF != NULL) << "failed to pot the Netlib LAPACKE_cpotrf symbol";
            
        cposvx_obj->inforef = CPOTRF( cposvx_obj->matrix_layout,
                                      cposvx_obj->uplo, cposvx_obj->n,
                                      cposvx_obj->afref,
                                      cposvx_obj->lda);
                               
        cposvx_obj->info = LAPACKE_cpotrf( cposvx_obj->matrix_layout,
                                           cposvx_obj->uplo, cposvx_obj->n,
                                           cposvx_obj->af,
                                           cposvx_obj->lda);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cposvx_obj->inforef = CPOSVX( cposvx_obj->matrix_layout, cposvx_obj->fact,
                                  cposvx_obj->uplo, cposvx_obj->n,
                                  cposvx_obj->nrhs,
                                  cposvx_obj->aref, cposvx_obj->lda, 
                                  cposvx_obj->afref, cposvx_obj->ldaf,
                                  &cposvx_obj->equedref,
                                  cposvx_obj->sref,                               
                                  cposvx_obj->bref, cposvx_obj->ldb,
                                  cposvx_obj->xref, cposvx_obj->ldx,
                                  &cposvx_obj->rcondref, 
                                  cposvx_obj->ferrref,
                                  cposvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cposvx_obj->info = LAPACKE_cposvx( cposvx_obj->matrix_layout, cposvx_obj->fact,
                                  cposvx_obj->uplo, cposvx_obj->n,
                                  cposvx_obj->nrhs,
                                  cposvx_obj->a, cposvx_obj->lda, 
                                  cposvx_obj->af, cposvx_obj->ldaf,
                                  &cposvx_obj->equed,
                                  cposvx_obj->s,                                  
                                  cposvx_obj->b, cposvx_obj->ldb,
                                  cposvx_obj->x, cposvx_obj->ldx,
                                  &cposvx_obj->rcond, 
                                  cposvx_obj->ferr,
                                  cposvx_obj->berr);

    if( cposvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cposvx is wrong\n", cposvx_obj->info );
    }
    if( cposvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cposvx is wrong\n", 
        cposvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cposvx_obj->diff =  computeDiff_c( cposvx_obj->b_bufsize, 
                cposvx_obj->b, cposvx_obj->bref );

    cposvx_obj->diff_xerr =  computeDiff_c( cposvx_obj->x_bufsize, 
                cposvx_obj->x, cposvx_obj->xref );

    cposvx_obj->diff_berr =  computeDiff_s( cposvx_obj->nrhs, 
                cposvx_obj->berr, cposvx_obj->berrref );
                
    cposvx_obj->diff_ferr =  computeDiff_s( cposvx_obj->nrhs, 
                cposvx_obj->ferr, cposvx_obj->ferrref );
}

TEST_F(cposvx_test, cposvx1) {
    EXPECT_NEAR(0.0, cposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cposvx_obj->rcond - cposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cposvx_test, cposvx2) {
    EXPECT_NEAR(0.0, cposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cposvx_obj->rcond - cposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cposvx_test, cposvx3) {
    EXPECT_NEAR(0.0, cposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cposvx_obj->rcond - cposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cposvx_test, cposvx4) {
    EXPECT_NEAR(0.0, cposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, cposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (cposvx_obj->rcond - cposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

/* Begin posvx_dcomplex_parameters  class definition */
class posvx_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

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
      posvx_dcomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~posvx_dcomplex_parameters (); 
};  /* end of posvx_dcomplex_parameters  class definition */


/* Constructor posvx_dcomplex_parameters definition */
posvx_dcomplex_parameters:: posvx_dcomplex_parameters ( int matrix_layout_i, 
                char fact_i, char uplo_i, char equed_i, lapack_int n_i,
                  lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    fact= fact_i;
    uplo = uplo_i;
    // equed acts as input when " fact = 'F' " else output
    equed = equed_i;
    equedref = equed_i;
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
   printf(" \n posvx dcomplex:  n: %d, fact: %c uplo: %c  lda: %d  \
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
       posvx_free();
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
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

posvx_dcomplex_parameters:: ~posvx_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" posvx_dcomplex_parameters object: destructor invoked. \n");
#endif
   posvx_free();
}

//  Test fixture class definition
class zposvx_test  : public  ::testing::Test {
public:
   posvx_dcomplex_parameters  *zposvx_obj;
   void SetUp();  
   void TearDown () { delete zposvx_obj; }
};


void zposvx_test::SetUp(){

    /* LAPACKE ZPOSVX prototype */
    typedef int (*Fptr_NL_LAPACKE_zposvx) (int matrix_layout,char fact, char uplo,
             lapack_int n, lapack_int nrhs, lapack_complex_double* a, lapack_int lda, lapack_complex_double* af,
             lapack_int ldaf, char* equed, double* s, lapack_complex_double* b, lapack_int ldb,
             lapack_complex_double* x, lapack_int ldx, double* rcond, double* ferr, double* berr);

    Fptr_NL_LAPACKE_zposvx ZPOSVX;

     /* LAPACKE ZPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpotrf) ( int matrix_layout , char uplo ,
                                lapack_int n , lapack_complex_double * a , lapack_int lda );

    Fptr_NL_LAPACKE_zpotrf ZPOTRF;

    zposvx_obj = new posvx_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    idx = Circular_Increment_Index(idx);

    zposvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zposvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zposvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zposvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZPOSVX = (Fptr_NL_LAPACKE_zposvx)dlsym(zposvx_obj->hModule, "LAPACKE_zposvx");
    ASSERT_TRUE(ZPOSVX != NULL) << "failed to pot the Netlib LAPACKE_zposvx symbol";
    
    /* if  'fact' == 'F', then the i/p matrix A should be in 
       factored form.
       so, invoke the potrf API to compute the factorized A.  */
    if(zposvx_obj->fact == 'F') {
        ZPOTRF = (Fptr_NL_LAPACKE_zpotrf)dlsym(zposvx_obj->hModule,"LAPACKE_zpotrf");
        ASSERT_TRUE(ZPOTRF != NULL) << "failed to pot the Netlib LAPACKE_zpotrf symbol";
            
        zposvx_obj->inforef = ZPOTRF( zposvx_obj->matrix_layout,
                                      zposvx_obj->uplo, zposvx_obj->n,
                                      zposvx_obj->afref,
                                      zposvx_obj->lda);
                               
        zposvx_obj->info = LAPACKE_zpotrf( zposvx_obj->matrix_layout,
                                           zposvx_obj->uplo, zposvx_obj->n,
                                           zposvx_obj->af,
                                           zposvx_obj->lda);
    }

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zposvx_obj->inforef = ZPOSVX( zposvx_obj->matrix_layout, zposvx_obj->fact,
                                  zposvx_obj->uplo, zposvx_obj->n,
                                  zposvx_obj->nrhs,
                                  zposvx_obj->aref, zposvx_obj->lda, 
                                  zposvx_obj->afref, zposvx_obj->ldaf,
                                  &zposvx_obj->equedref,
                                  zposvx_obj->sref,                               
                                  zposvx_obj->bref, zposvx_obj->ldb,
                                  zposvx_obj->xref, zposvx_obj->ldx,
                                  &zposvx_obj->rcondref, 
                                  zposvx_obj->ferrref,
                                  zposvx_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zposvx_obj->info = LAPACKE_zposvx( zposvx_obj->matrix_layout, zposvx_obj->fact,
                                  zposvx_obj->uplo, zposvx_obj->n,
                                  zposvx_obj->nrhs,
                                  zposvx_obj->a, zposvx_obj->lda, 
                                  zposvx_obj->af, zposvx_obj->ldaf,
                                  &zposvx_obj->equed,
                                  zposvx_obj->s,                                  
                                  zposvx_obj->b, zposvx_obj->ldb,
                                  zposvx_obj->x, zposvx_obj->ldx,
                                  &zposvx_obj->rcond, 
                                  zposvx_obj->ferr,
                                  zposvx_obj->berr);

    if( zposvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zposvx is wrong\n", zposvx_obj->info );
    }
    if( zposvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zposvx is wrong\n", 
        zposvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zposvx_obj->diff =  computeDiff_z( zposvx_obj->b_bufsize, 
                zposvx_obj->b, zposvx_obj->bref );

    zposvx_obj->diff_xerr =  computeDiff_z( zposvx_obj->x_bufsize, 
                zposvx_obj->x, zposvx_obj->xref );

    zposvx_obj->diff_berr =  computeDiff_d( zposvx_obj->nrhs, 
                zposvx_obj->berr, zposvx_obj->berrref );
                
    zposvx_obj->diff_ferr =  computeDiff_d( zposvx_obj->nrhs, 
                zposvx_obj->ferr, zposvx_obj->ferrref );
}

TEST_F(zposvx_test, zposvx1) {
    EXPECT_NEAR(0.0, zposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zposvx_obj->rcond - zposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zposvx_test, zposvx2) {
    EXPECT_NEAR(0.0, zposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zposvx_obj->rcond - zposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zposvx_test, zposvx3) {
    EXPECT_NEAR(0.0, zposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zposvx_obj->rcond - zposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zposvx_test, zposvx4) {
    EXPECT_NEAR(0.0, zposvx_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zposvx_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, (zposvx_obj->rcond - zposvx_obj->rcondref),
                LAPACKE_GTEST_THRESHOLD);
}