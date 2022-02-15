#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define porfs_free() \
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
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin porfs_float_parameters  class definition */
class porfs_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
	  float threshold;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      porfs_float_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~porfs_float_parameters (); 
};  /* end of porfs_float_parameters  class definition */


/* Constructor porfs_float_parameters definition */
porfs_float_parameters:: porfs_float_parameters ( int matrix_layout_i, 
                 char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
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
   printf(" \n porfs float:  matrix_layout_i: %d n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  threshold: %f \n", matrix_layout_i, n, uplo, lda,
ldb, nrhs, lin_solver_paramslist[idx].solver_threhold );
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       porfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(float)));
    memcpy(afref, a, (n*n*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(float)));
    memcpy(xref, b, (b_bufsize*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

porfs_float_parameters:: ~porfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" porfs_float_parameters object: destructor invoked. \n");
#endif
   porfs_free();
}

//  Test fixture class definition
class sporfs_test  : public  ::testing::Test {
public:
   porfs_float_parameters  *sporfs_obj;
   void SetUp();  
   void TearDown () { delete sporfs_obj; }
};


void sporfs_test::SetUp(){

    /* LAPACKE SPORFS prototype */
    typedef int (*Fptr_NL_LAPACKE_sporfs) (int matrix_layout, char uplo,
        lapack_int n, lapack_int nrhs, const float* a, lapack_int lda,
        const float* af, lapack_int ldaf, const float* b, lapack_int ldb,
        float* x, lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_sporfs SPORFS;

     /* LAPACKE SPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spotrf) ( int matrix_layout , char uplo ,
                                lapack_int n , float * a , lapack_int lda );

    Fptr_NL_LAPACKE_spotrf SPOTRF;

    /* LAPACKE SPOTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_spotrs) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs, const float * a,
                          lapack_int lda, float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_spotrs SPOTRS;

    sporfs_obj = new porfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

	sporfs_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
    idx = Circular_Increment_Index(idx);

    sporfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sporfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sporfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sporfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SPORFS = (Fptr_NL_LAPACKE_sporfs)dlsym(sporfs_obj->hModule, "LAPACKE_sporfs");
    ASSERT_TRUE(SPORFS != NULL) << "failed to pot the Netlib LAPACKE_sporfs symbol";

    SPOTRS = (Fptr_NL_LAPACKE_spotrs)dlsym(sporfs_obj->hModule, "LAPACKE_spotrs");
    ASSERT_TRUE(SPOTRS != NULL) << "failed to get the Netlib LAPACKE_spotrs symbol";
    
    SPOTRF = (Fptr_NL_LAPACKE_spotrf)dlsym(sporfs_obj->hModule,"LAPACKE_spotrf");
    ASSERT_TRUE(SPOTRF != NULL) << "failed to pot the Netlib LAPACKE_spotrf symbol";
    /*  Generate the i/ps for porfs by calling the sequence of potrf, potrs APIs  */
            
    sporfs_obj->inforef = SPOTRF( sporfs_obj->matrix_layout,
                                  sporfs_obj->uplo, sporfs_obj->n,
                                  sporfs_obj->afref,
                                  sporfs_obj->lda);
                                  
    sporfs_obj->inforef = SPOTRS( sporfs_obj->matrix_layout,
                                  sporfs_obj->uplo, 
                                  sporfs_obj->n,
                                  sporfs_obj->nrhs,
                                  (const float *)sporfs_obj->afref, 
                                  sporfs_obj->lda,
                                  sporfs_obj->xref,
                                  sporfs_obj->ldb);
                          
    // invoking netlib potrf for generation of libflame i/p as well.                       
    sporfs_obj->info = SPOTRF( sporfs_obj->matrix_layout,
                               sporfs_obj->uplo, sporfs_obj->n,
                               sporfs_obj->af,
                               sporfs_obj->lda);

    sporfs_obj->info = SPOTRS( sporfs_obj->matrix_layout,
                               sporfs_obj->uplo, sporfs_obj->n,
                               sporfs_obj->nrhs,
                               (const float *)sporfs_obj->af, 
                                sporfs_obj->lda,
                                sporfs_obj->x,
                                sporfs_obj->ldb );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    sporfs_obj->inforef = SPORFS( sporfs_obj->matrix_layout,
                                  sporfs_obj->uplo, sporfs_obj->n,
                                  sporfs_obj->nrhs,
                                  sporfs_obj->aref, sporfs_obj->lda, 
                                  sporfs_obj->afref, sporfs_obj->ldaf,
                                  sporfs_obj->bref, sporfs_obj->ldb,
                                  sporfs_obj->xref, sporfs_obj->ldx,
                                  sporfs_obj->ferrref,
                                  sporfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    sporfs_obj->info = LAPACKE_sporfs( sporfs_obj->matrix_layout,
                                  sporfs_obj->uplo, sporfs_obj->n,
                                  sporfs_obj->nrhs,
                                  sporfs_obj->a, sporfs_obj->lda, 
                                  sporfs_obj->af, sporfs_obj->ldaf,
                                  sporfs_obj->b, sporfs_obj->ldb,
                                  sporfs_obj->x, sporfs_obj->ldx,
                                  sporfs_obj->ferr,
                                  sporfs_obj->berr);

    if( sporfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sporfs is wrong\n", sporfs_obj->info );
    }
    if( sporfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sporfs is wrong\n", 
        sporfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sporfs_obj->diff_xerr =  computeDiff_s( sporfs_obj->x_bufsize, 
                sporfs_obj->x, sporfs_obj->xref );

    sporfs_obj->diff_berr =  computeDiff_s( sporfs_obj->nrhs, 
                sporfs_obj->berr, sporfs_obj->berrref );
                
    sporfs_obj->diff_ferr =  computeDiff_s( sporfs_obj->nrhs, 
                sporfs_obj->ferr, sporfs_obj->ferrref );
}

TEST_F(sporfs_test, sporfs1) {
    EXPECT_NEAR(0.0, sporfs_obj->diff_xerr, sporfs_obj->threshold);
    EXPECT_NEAR(0.0, sporfs_obj->diff_berr, sporfs_obj->threshold);
    EXPECT_NEAR(0.0, sporfs_obj->diff_ferr, sporfs_obj->threshold);
}

TEST_F(sporfs_test, sporfs2) {
    EXPECT_NEAR(0.0, sporfs_obj->diff_xerr, sporfs_obj->threshold);
    EXPECT_NEAR(0.0, sporfs_obj->diff_berr, sporfs_obj->threshold);
    EXPECT_NEAR(0.0, sporfs_obj->diff_ferr, sporfs_obj->threshold);
}

TEST_F(sporfs_test, sporfs3) {
    EXPECT_NEAR(0.0, sporfs_obj->diff_xerr, sporfs_obj->threshold);
    EXPECT_NEAR(0.0, sporfs_obj->diff_berr, sporfs_obj->threshold);
    EXPECT_NEAR(0.0, sporfs_obj->diff_ferr, sporfs_obj->threshold);
}

TEST_F(sporfs_test, sporfs4) {
    EXPECT_NEAR(0.0, sporfs_obj->diff_xerr, sporfs_obj->threshold);
    EXPECT_NEAR(0.0, sporfs_obj->diff_berr, sporfs_obj->threshold);
    EXPECT_NEAR(0.0, sporfs_obj->diff_ferr, sporfs_obj->threshold);
}

/* Begin porfs_double_parameters  class definition */
class porfs_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
	  float threshold;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      porfs_double_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~porfs_double_parameters (); 
};  /* end of porfs_double_parameters  class definition */


/* Constructor porfs_double_parameters definition */
porfs_double_parameters:: porfs_double_parameters ( int matrix_layout_i, 
                 char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
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
   printf(" \n porfs double:  matrix_layout_i: %d n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  threshold: %f \n", matrix_layout_i, n, uplo, lda,
ldb, nrhs, lin_solver_paramslist[idx].solver_threhold );
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       porfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(double)));
    memcpy(afref, a, (n*n*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(double)));
    memcpy(xref, b, (b_bufsize*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

porfs_double_parameters:: ~porfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" porfs_double_parameters object: destructor invoked. \n");
#endif
   porfs_free();
}

//  Test fixture class definition
class dporfs_test  : public  ::testing::Test {
public:
   porfs_double_parameters  *dporfs_obj;
   void SetUp();  
   void TearDown () { delete dporfs_obj; }
};


void dporfs_test::SetUp(){

    /* LAPACKE DPORFS prototype */
    typedef int (*Fptr_NL_LAPACKE_dporfs) (int matrix_layout, char uplo,
        lapack_int n, lapack_int nrhs, const double* a, lapack_int lda,
        const double* af, lapack_int ldaf, const double* b, lapack_int ldb,
        double* x, lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_dporfs DPORFS;

     /* LAPACKE DPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpotrf) ( int matrix_layout , char uplo ,
                                lapack_int n , double * a , lapack_int lda );

    Fptr_NL_LAPACKE_dpotrf DPOTRF;

    /* LAPACKE DPOTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dpotrs) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs, const double * a,
                          lapack_int lda, double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dpotrs DPOTRS;

    dporfs_obj = new porfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

	dporfs_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
    idx = Circular_Increment_Index(idx);

    dporfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dporfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dporfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dporfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DPORFS = (Fptr_NL_LAPACKE_dporfs)dlsym(dporfs_obj->hModule, "LAPACKE_dporfs");
    ASSERT_TRUE(DPORFS != NULL) << "failed to pot the Netlib LAPACKE_dporfs symbol";

    DPOTRS = (Fptr_NL_LAPACKE_dpotrs)dlsym(dporfs_obj->hModule, "LAPACKE_dpotrs");
    ASSERT_TRUE(DPOTRS != NULL) << "failed to get the Netlib LAPACKE_dpotrs symbol";
    
    DPOTRF = (Fptr_NL_LAPACKE_dpotrf)dlsym(dporfs_obj->hModule,"LAPACKE_dpotrf");
    ASSERT_TRUE(DPOTRF != NULL) << "failed to pot the Netlib LAPACKE_dpotrf symbol";
    /*  Generate the i/ps for porfs by calling the sequence of potrf, potrs APIs  */
            
    dporfs_obj->inforef = DPOTRF( dporfs_obj->matrix_layout,
                                  dporfs_obj->uplo, dporfs_obj->n,
                                  dporfs_obj->afref,
                                  dporfs_obj->lda);
                                  
    dporfs_obj->inforef = DPOTRS( dporfs_obj->matrix_layout,
                                  dporfs_obj->uplo, 
                                  dporfs_obj->n,
                                  dporfs_obj->nrhs,
                                  (const double *)dporfs_obj->afref, 
                                  dporfs_obj->lda,
                                  dporfs_obj->xref,
                                  dporfs_obj->ldb);
                          
    // invoking netlib potrf for generation of libflame i/p as well.                       
    dporfs_obj->info = DPOTRF( dporfs_obj->matrix_layout,
                               dporfs_obj->uplo, dporfs_obj->n,
                               dporfs_obj->af,
                               dporfs_obj->lda);

    dporfs_obj->info = DPOTRS( dporfs_obj->matrix_layout,
                               dporfs_obj->uplo, dporfs_obj->n,
                               dporfs_obj->nrhs,
                               (const double *)dporfs_obj->af, 
                                dporfs_obj->lda,
                                dporfs_obj->x,
                                dporfs_obj->ldb );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dporfs_obj->inforef = DPORFS( dporfs_obj->matrix_layout,
                                  dporfs_obj->uplo, dporfs_obj->n,
                                  dporfs_obj->nrhs,
                                  dporfs_obj->aref, dporfs_obj->lda, 
                                  dporfs_obj->afref, dporfs_obj->ldaf,
                                  dporfs_obj->bref, dporfs_obj->ldb,
                                  dporfs_obj->xref, dporfs_obj->ldx,
                                  dporfs_obj->ferrref,
                                  dporfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dporfs_obj->info = LAPACKE_dporfs( dporfs_obj->matrix_layout,
                                  dporfs_obj->uplo, dporfs_obj->n,
                                  dporfs_obj->nrhs,
                                  dporfs_obj->a, dporfs_obj->lda, 
                                  dporfs_obj->af, dporfs_obj->ldaf,
                                  dporfs_obj->b, dporfs_obj->ldb,
                                  dporfs_obj->x, dporfs_obj->ldx,
                                  dporfs_obj->ferr,
                                  dporfs_obj->berr);

    if( dporfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dporfs is wrong\n", dporfs_obj->info );
    }
    if( dporfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dporfs is wrong\n", 
        dporfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dporfs_obj->diff_xerr =  computeDiff_d( dporfs_obj->x_bufsize, 
                dporfs_obj->x, dporfs_obj->xref );

    dporfs_obj->diff_berr =  computeDiff_d( dporfs_obj->nrhs, 
                dporfs_obj->berr, dporfs_obj->berrref );
                
    dporfs_obj->diff_ferr =  computeDiff_d( dporfs_obj->nrhs, 
                dporfs_obj->ferr, dporfs_obj->ferrref );
}

TEST_F(dporfs_test, dporfs1) {
    EXPECT_NEAR(0.0, dporfs_obj->diff_xerr, dporfs_obj->threshold);
    EXPECT_NEAR(0.0, dporfs_obj->diff_berr, dporfs_obj->threshold);
    EXPECT_NEAR(0.0, dporfs_obj->diff_ferr, dporfs_obj->threshold);
}

TEST_F(dporfs_test, dporfs2) {
    EXPECT_NEAR(0.0, dporfs_obj->diff_xerr, dporfs_obj->threshold);
    EXPECT_NEAR(0.0, dporfs_obj->diff_berr, dporfs_obj->threshold);
    EXPECT_NEAR(0.0, dporfs_obj->diff_ferr, dporfs_obj->threshold);
}

TEST_F(dporfs_test, dporfs3) {
    EXPECT_NEAR(0.0, dporfs_obj->diff_xerr, dporfs_obj->threshold);
    EXPECT_NEAR(0.0, dporfs_obj->diff_berr, dporfs_obj->threshold);
    EXPECT_NEAR(0.0, dporfs_obj->diff_ferr, dporfs_obj->threshold);
}

TEST_F(dporfs_test, dporfs4) {
    EXPECT_NEAR(0.0, dporfs_obj->diff_xerr, dporfs_obj->threshold);
    EXPECT_NEAR(0.0, dporfs_obj->diff_berr, dporfs_obj->threshold);
    EXPECT_NEAR(0.0, dporfs_obj->diff_ferr, dporfs_obj->threshold);
}

/* Begin porfs_scomplex_parameters  class definition */
class porfs_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
	  float threshold;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      porfs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~porfs_scomplex_parameters (); 
};  /* end of porfs_scomplex_parameters  class definition */


/* Constructor porfs_scomplex_parameters definition */
porfs_scomplex_parameters:: porfs_scomplex_parameters ( int matrix_layout_i, 
                 char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
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
   printf(" \n porfs lapack_complex_float:  matrix_layout_i: %d n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  threshold: %f \n", matrix_layout_i, n, uplo, lda,
ldb, nrhs, lin_solver_paramslist[idx].solver_threhold );
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       porfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(lapack_complex_float)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, (b_bufsize*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

porfs_scomplex_parameters:: ~porfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" porfs_scomplex_parameters object: destructor invoked. \n");
#endif
   porfs_free();
}

//  Test fixture class definition
class cporfs_test  : public  ::testing::Test {
public:
   porfs_scomplex_parameters  *cporfs_obj;
   void SetUp();  
   void TearDown () { delete cporfs_obj; }
};


void cporfs_test::SetUp(){

    /* LAPACKE CPORFS prototype */
    typedef int (*Fptr_NL_LAPACKE_cporfs) (int matrix_layout, char uplo,
        lapack_int n, lapack_int nrhs, const lapack_complex_float* a, lapack_int lda,
        const lapack_complex_float* af, lapack_int ldaf, const lapack_complex_float* b, lapack_int ldb,
        lapack_complex_float* x, lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_cporfs CPORFS;

     /* LAPACKE CPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpotrf) ( int matrix_layout , char uplo ,
                                lapack_int n , lapack_complex_float * a , lapack_int lda );

    Fptr_NL_LAPACKE_cpotrf CPOTRF;

    /* LAPACKE CPOTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_cpotrs) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs, const lapack_complex_float * a,
                          lapack_int lda, lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_cpotrs CPOTRS;

    cporfs_obj = new porfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

	cporfs_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
    idx = Circular_Increment_Index(idx);

    cporfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cporfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cporfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cporfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CPORFS = (Fptr_NL_LAPACKE_cporfs)dlsym(cporfs_obj->hModule, "LAPACKE_cporfs");
    ASSERT_TRUE(CPORFS != NULL) << "failed to pot the Netlib LAPACKE_cporfs symbol";

    CPOTRS = (Fptr_NL_LAPACKE_cpotrs)dlsym(cporfs_obj->hModule, "LAPACKE_cpotrs");
    ASSERT_TRUE(CPOTRS != NULL) << "failed to get the Netlib LAPACKE_cpotrs symbol";
    
    CPOTRF = (Fptr_NL_LAPACKE_cpotrf)dlsym(cporfs_obj->hModule,"LAPACKE_cpotrf");
    ASSERT_TRUE(CPOTRF != NULL) << "failed to pot the Netlib LAPACKE_cpotrf symbol";
    /*  Generate the i/ps for porfs by calling the sequence of potrf, potrs APIs  */
            
    cporfs_obj->inforef = CPOTRF( cporfs_obj->matrix_layout,
                                  cporfs_obj->uplo, cporfs_obj->n,
                                  cporfs_obj->afref,
                                  cporfs_obj->lda);
                                  
    cporfs_obj->inforef = CPOTRS( cporfs_obj->matrix_layout,
                                  cporfs_obj->uplo, 
                                  cporfs_obj->n,
                                  cporfs_obj->nrhs,
                                  (const lapack_complex_float *)cporfs_obj->afref, 
                                  cporfs_obj->lda,
                                  cporfs_obj->xref,
                                  cporfs_obj->ldb);
                          
    // invoking netlib potrf for generation of libflame i/p as well.                       
    cporfs_obj->info = CPOTRF( cporfs_obj->matrix_layout,
                               cporfs_obj->uplo, cporfs_obj->n,
                               cporfs_obj->af,
                               cporfs_obj->lda);

    cporfs_obj->info = CPOTRS( cporfs_obj->matrix_layout,
                               cporfs_obj->uplo, cporfs_obj->n,
                               cporfs_obj->nrhs,
                               (const lapack_complex_float *)cporfs_obj->af, 
                                cporfs_obj->lda,
                                cporfs_obj->x,
                                cporfs_obj->ldb );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cporfs_obj->inforef = CPORFS( cporfs_obj->matrix_layout,
                                  cporfs_obj->uplo, cporfs_obj->n,
                                  cporfs_obj->nrhs,
                                  cporfs_obj->aref, cporfs_obj->lda, 
                                  cporfs_obj->afref, cporfs_obj->ldaf,
                                  cporfs_obj->bref, cporfs_obj->ldb,
                                  cporfs_obj->xref, cporfs_obj->ldx,
                                  cporfs_obj->ferrref,
                                  cporfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cporfs_obj->info = LAPACKE_cporfs( cporfs_obj->matrix_layout,
                                  cporfs_obj->uplo, cporfs_obj->n,
                                  cporfs_obj->nrhs,
                                  cporfs_obj->a, cporfs_obj->lda, 
                                  cporfs_obj->af, cporfs_obj->ldaf,
                                  cporfs_obj->b, cporfs_obj->ldb,
                                  cporfs_obj->x, cporfs_obj->ldx,
                                  cporfs_obj->ferr,
                                  cporfs_obj->berr);

    if( cporfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cporfs is wrong\n", cporfs_obj->info );
    }
    if( cporfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cporfs is wrong\n", 
        cporfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cporfs_obj->diff_xerr =  computeDiff_c( cporfs_obj->x_bufsize, 
                cporfs_obj->x, cporfs_obj->xref );

    cporfs_obj->diff_berr =  computeDiff_s( cporfs_obj->nrhs, 
                cporfs_obj->berr, cporfs_obj->berrref );
                
    cporfs_obj->diff_ferr =  computeDiff_s( cporfs_obj->nrhs, 
                cporfs_obj->ferr, cporfs_obj->ferrref );
}

TEST_F(cporfs_test, cporfs1) {
    EXPECT_NEAR(0.0, cporfs_obj->diff_xerr, cporfs_obj->threshold);
    EXPECT_NEAR(0.0, cporfs_obj->diff_berr, cporfs_obj->threshold);
    EXPECT_NEAR(0.0, cporfs_obj->diff_ferr, cporfs_obj->threshold);
}

TEST_F(cporfs_test, cporfs2) {
    EXPECT_NEAR(0.0, cporfs_obj->diff_xerr, cporfs_obj->threshold);
    EXPECT_NEAR(0.0, cporfs_obj->diff_berr, cporfs_obj->threshold);
    EXPECT_NEAR(0.0, cporfs_obj->diff_ferr, cporfs_obj->threshold);
}

TEST_F(cporfs_test, cporfs3) {
    EXPECT_NEAR(0.0, cporfs_obj->diff_xerr, cporfs_obj->threshold);
    EXPECT_NEAR(0.0, cporfs_obj->diff_berr, cporfs_obj->threshold);
    EXPECT_NEAR(0.0, cporfs_obj->diff_ferr, cporfs_obj->threshold);
}

TEST_F(cporfs_test, cporfs4) {
    EXPECT_NEAR(0.0, cporfs_obj->diff_xerr, cporfs_obj->threshold);
    EXPECT_NEAR(0.0, cporfs_obj->diff_berr, cporfs_obj->threshold);
    EXPECT_NEAR(0.0, cporfs_obj->diff_ferr, cporfs_obj->threshold);
}

/* Begin porfs_dcomplex_parameters  class definition */
class porfs_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
	  double threshold;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      porfs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~porfs_dcomplex_parameters (); 
};  /* end of porfs_dcomplex_parameters  class definition */


/* Constructor porfs_dcomplex_parameters definition */
porfs_dcomplex_parameters:: porfs_dcomplex_parameters ( int matrix_layout_i, 
                 char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
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
   printf(" \n porfs lapack_complex_double:  matrix_layout_i: %d n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  threshold: %f \n", matrix_layout_i, n, uplo, lda,
ldb, nrhs, lin_solver_paramslist[idx].solver_threhold );
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       porfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix ( a, aref, n, n, uplo);
    memcpy(af, a, (n*n*sizeof(lapack_complex_double)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, (b_bufsize*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

porfs_dcomplex_parameters:: ~porfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" porfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   porfs_free();
}

//  Test fixture class definition
class zporfs_test  : public  ::testing::Test {
public:
   porfs_dcomplex_parameters  *zporfs_obj;
   void SetUp();  
   void TearDown () { delete zporfs_obj; }
};


void zporfs_test::SetUp(){

    /* LAPACKE ZPORFS prototype */
    typedef int (*Fptr_NL_LAPACKE_zporfs) (int matrix_layout, char uplo,
        lapack_int n, lapack_int nrhs, const lapack_complex_double* a, lapack_int lda,
        const lapack_complex_double* af, lapack_int ldaf, const lapack_complex_double* b, lapack_int ldb,
        lapack_complex_double* x, lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_zporfs ZPORFS;

     /* LAPACKE ZPOTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpotrf) ( int matrix_layout , char uplo ,
                                lapack_int n , lapack_complex_double * a , lapack_int lda );

    Fptr_NL_LAPACKE_zpotrf ZPOTRF;

    /* LAPACKE ZPOTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zpotrs) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs, const lapack_complex_double * a,
                          lapack_int lda, lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zpotrs ZPOTRS;

    zporfs_obj = new porfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

	zporfs_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
    idx = Circular_Increment_Index(idx);

    zporfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zporfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zporfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zporfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZPORFS = (Fptr_NL_LAPACKE_zporfs)dlsym(zporfs_obj->hModule, "LAPACKE_zporfs");
    ASSERT_TRUE(ZPORFS != NULL) << "failed to pot the Netlib LAPACKE_zporfs symbol";

    ZPOTRS = (Fptr_NL_LAPACKE_zpotrs)dlsym(zporfs_obj->hModule, "LAPACKE_zpotrs");
    ASSERT_TRUE(ZPOTRS != NULL) << "failed to get the Netlib LAPACKE_zpotrs symbol";
    
    ZPOTRF = (Fptr_NL_LAPACKE_zpotrf)dlsym(zporfs_obj->hModule,"LAPACKE_zpotrf");
    ASSERT_TRUE(ZPOTRF != NULL) << "failed to pot the Netlib LAPACKE_zpotrf symbol";
    /*  Generate the i/ps for porfs by calling the sequence of potrf, potrs APIs  */
            
    zporfs_obj->inforef = ZPOTRF( zporfs_obj->matrix_layout,
                                  zporfs_obj->uplo, zporfs_obj->n,
                                  zporfs_obj->afref,
                                  zporfs_obj->lda);
                                  
    zporfs_obj->inforef = ZPOTRS( zporfs_obj->matrix_layout,
                                  zporfs_obj->uplo, 
                                  zporfs_obj->n,
                                  zporfs_obj->nrhs,
                                  (const lapack_complex_double *)zporfs_obj->afref, 
                                  zporfs_obj->lda,
                                  zporfs_obj->xref,
                                  zporfs_obj->ldb);
                          
    // invoking netlib potrf for generation of libflame i/p as well.                       
    zporfs_obj->info = ZPOTRF( zporfs_obj->matrix_layout,
                               zporfs_obj->uplo, zporfs_obj->n,
                               zporfs_obj->af,
                               zporfs_obj->lda);

    zporfs_obj->info = ZPOTRS( zporfs_obj->matrix_layout,
                               zporfs_obj->uplo, zporfs_obj->n,
                               zporfs_obj->nrhs,
                               (const lapack_complex_double *)zporfs_obj->af, 
                                zporfs_obj->lda,
                                zporfs_obj->x,
                                zporfs_obj->ldb );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zporfs_obj->inforef = ZPORFS( zporfs_obj->matrix_layout,
                                  zporfs_obj->uplo, zporfs_obj->n,
                                  zporfs_obj->nrhs,
                                  zporfs_obj->aref, zporfs_obj->lda, 
                                  zporfs_obj->afref, zporfs_obj->ldaf,
                                  zporfs_obj->bref, zporfs_obj->ldb,
                                  zporfs_obj->xref, zporfs_obj->ldx,
                                  zporfs_obj->ferrref,
                                  zporfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zporfs_obj->info = LAPACKE_zporfs( zporfs_obj->matrix_layout,
                                  zporfs_obj->uplo, zporfs_obj->n,
                                  zporfs_obj->nrhs,
                                  zporfs_obj->a, zporfs_obj->lda, 
                                  zporfs_obj->af, zporfs_obj->ldaf,
                                  zporfs_obj->b, zporfs_obj->ldb,
                                  zporfs_obj->x, zporfs_obj->ldx,
                                  zporfs_obj->ferr,
                                  zporfs_obj->berr);

    if( zporfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zporfs is wrong\n", zporfs_obj->info );
    }
    if( zporfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zporfs is wrong\n", 
        zporfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zporfs_obj->diff_xerr =  computeDiff_z( zporfs_obj->x_bufsize, 
                zporfs_obj->x, zporfs_obj->xref );

    zporfs_obj->diff_berr =  computeDiff_d( zporfs_obj->nrhs, 
                zporfs_obj->berr, zporfs_obj->berrref );
                
    zporfs_obj->diff_ferr =  computeDiff_d( zporfs_obj->nrhs, 
                zporfs_obj->ferr, zporfs_obj->ferrref );
}

TEST_F(zporfs_test, zporfs1) {
    EXPECT_NEAR(0.0, zporfs_obj->diff_xerr, zporfs_obj->threshold);
    EXPECT_NEAR(0.0, zporfs_obj->diff_berr, zporfs_obj->threshold);
    EXPECT_NEAR(0.0, zporfs_obj->diff_ferr, zporfs_obj->threshold);
}

TEST_F(zporfs_test, zporfs2) {
    EXPECT_NEAR(0.0, zporfs_obj->diff_xerr, zporfs_obj->threshold);
    EXPECT_NEAR(0.0, zporfs_obj->diff_berr, zporfs_obj->threshold);
    EXPECT_NEAR(0.0, zporfs_obj->diff_ferr, zporfs_obj->threshold);
}

TEST_F(zporfs_test, zporfs3) {
    EXPECT_NEAR(0.0, zporfs_obj->diff_xerr, zporfs_obj->threshold);
    EXPECT_NEAR(0.0, zporfs_obj->diff_berr, zporfs_obj->threshold);
    EXPECT_NEAR(0.0, zporfs_obj->diff_ferr, zporfs_obj->threshold);
}

TEST_F(zporfs_test, zporfs4) {
    EXPECT_NEAR(0.0, zporfs_obj->diff_xerr, zporfs_obj->threshold);
    EXPECT_NEAR(0.0, zporfs_obj->diff_berr, zporfs_obj->threshold);
    EXPECT_NEAR(0.0, zporfs_obj->diff_ferr, zporfs_obj->threshold);
}