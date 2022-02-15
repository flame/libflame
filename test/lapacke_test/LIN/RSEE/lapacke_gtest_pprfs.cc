
#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define pprfs_free() \
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

/* Begin pprfs_float_parameters  class definition */
class pprfs_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr, diff_s;
      void *hModule, *dModule;
      float threshold;

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
      pprfs_float_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~pprfs_float_parameters (); 
};  /* end of pprfs_float_parameters  class definition */


/* Constructor pprfs_float_parameters definition */
pprfs_float_parameters:: pprfs_float_parameters ( int matrix_layout_i, 
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
   printf(" \n pprfs float:  n: %d, fact: %c uplo: %c   \
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
       pprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, a_buf_size);
    memcpy(af, a, (a_buf_size*sizeof(float)));
    memcpy(afref, a, (a_buf_size*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(float)));
    memcpy(xref, b, (b_bufsize*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

pprfs_float_parameters:: ~pprfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pprfs_float_parameters object: destructor invoked. \n");
#endif
   pprfs_free();
}

//  Test fixture class definition
class spprfs_test  : public  ::testing::Test {
public:
   pprfs_float_parameters  *spprfs_obj;
   void SetUp();  
   void TearDown () { delete spprfs_obj; }
};


void spprfs_test::SetUp(){

    /* LAPACKE SPPRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_spprfs) (int matrix_layout, 
		char uplo, lapack_int n, lapack_int nrhs, const float* ap,
		const float* afp, const float* b, lapack_int ldb,
		float* x, lapack_int ldx, float* ferr, float* berr );

    Fptr_NL_LAPACKE_spprfs SPPRFS;

    /* LAPACKE SPPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_spptrs) ( int matrix_layout, 
					char uplo, lapack_int n, lapack_int nrhs,
              const float * a, float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_spptrs SPPTRS;

     /* LAPACKE SPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spptrf) ( int matrix_layout, 
           char uplo , lapack_int n , float * a );

    Fptr_NL_LAPACKE_spptrf SPPTRF;

    spprfs_obj = new pprfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);
    spprfs_obj->threshold = lin_solver_paramslist[idx].solver_threhold;

    idx = Circular_Increment_Index(idx);

    spprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SPPRFS = (Fptr_NL_LAPACKE_spprfs)dlsym(spprfs_obj->hModule, "LAPACKE_spprfs");
    ASSERT_TRUE(SPPRFS != NULL) << "failed to ppt the Netlib LAPACKE_spprfs symbol";

    SPPTRS = (Fptr_NL_LAPACKE_spptrs)dlsym(spprfs_obj->hModule, "LAPACKE_spptrs");
    ASSERT_TRUE(SPPTRS != NULL) << "failed to get the Netlib LAPACKE_spptrs symbol";

    SPPTRF = (Fptr_NL_LAPACKE_spptrf)dlsym(spprfs_obj->hModule,"LAPACKE_spptrf");
    ASSERT_TRUE(SPPTRF != NULL) << "failed to get the Netlib LAPACKE_spptrf symbol";
    
    /* invoke the pptrf API to compute the factorized A.  */
        spprfs_obj->inforef = SPPTRF( spprfs_obj->matrix_layout,
                                      spprfs_obj->uplo, spprfs_obj->n,
                                      spprfs_obj->afref);
                               
        spprfs_obj->info = LAPACKE_spptrf( spprfs_obj->matrix_layout,
                                           spprfs_obj->uplo, spprfs_obj->n,
                                           spprfs_obj->af);

    spprfs_obj->inforef = SPPTRS( spprfs_obj->matrix_layout, 
          spprfs_obj->uplo, spprfs_obj->n,
         spprfs_obj->nrhs, (const float *)spprfs_obj->afref, 
                          spprfs_obj->xref, spprfs_obj->ldx);

    spprfs_obj->info = LAPACKE_spptrs( spprfs_obj->matrix_layout, 
               spprfs_obj->uplo, spprfs_obj->n,
                 spprfs_obj->nrhs, (const float *)spprfs_obj->af, 
                                 spprfs_obj->x, spprfs_obj->ldx );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
	
    spprfs_obj->inforef = SPPRFS( spprfs_obj->matrix_layout,
                                  spprfs_obj->uplo, spprfs_obj->n,
                                  spprfs_obj->nrhs,
                                  spprfs_obj->aref, 
                                  spprfs_obj->afref,
                                  spprfs_obj->bref, spprfs_obj->ldb,
                                  spprfs_obj->xref, spprfs_obj->ldx,
                                  spprfs_obj->ferrref,
                                  spprfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    spprfs_obj->info = LAPACKE_spprfs( spprfs_obj->matrix_layout,
                                  spprfs_obj->uplo, spprfs_obj->n,
                                  spprfs_obj->nrhs,
                                  spprfs_obj->a, 
                                  spprfs_obj->af,
                                  spprfs_obj->b, spprfs_obj->ldb,
                                  spprfs_obj->x, spprfs_obj->ldx,
                                  spprfs_obj->ferr,
                                  spprfs_obj->berr);

    if( spprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_spprfs is wrong\n", spprfs_obj->info );
    }
    if( spprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spprfs is wrong\n", 
        spprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    spprfs_obj->diff_xerr =  computeDiff_s( spprfs_obj->x_bufsize, 
                spprfs_obj->x, spprfs_obj->xref );

    spprfs_obj->diff_berr =  computeDiff_s( spprfs_obj->nrhs, 
                spprfs_obj->berr, spprfs_obj->berrref );
                
    spprfs_obj->diff_ferr =  computeDiff_s( spprfs_obj->nrhs, 
                spprfs_obj->ferr, spprfs_obj->ferrref );
    
}

TEST_F(spprfs_test, spprfs1) {
    EXPECT_NEAR(0.0, spprfs_obj->diff_xerr, spprfs_obj->threshold);
    EXPECT_NEAR(0.0, spprfs_obj->diff_berr, spprfs_obj->threshold);
    EXPECT_NEAR(0.0, spprfs_obj->diff_ferr, spprfs_obj->threshold);
}

TEST_F(spprfs_test, spprfs2) {
    EXPECT_NEAR(0.0, spprfs_obj->diff_xerr, spprfs_obj->threshold);
    EXPECT_NEAR(0.0, spprfs_obj->diff_berr, spprfs_obj->threshold);
    EXPECT_NEAR(0.0, spprfs_obj->diff_ferr, spprfs_obj->threshold);
}

TEST_F(spprfs_test, spprfs3) {
    EXPECT_NEAR(0.0, spprfs_obj->diff_xerr, spprfs_obj->threshold);
    EXPECT_NEAR(0.0, spprfs_obj->diff_berr, spprfs_obj->threshold);
    EXPECT_NEAR(0.0, spprfs_obj->diff_ferr, spprfs_obj->threshold);
}

TEST_F(spprfs_test, spprfs4) {
    EXPECT_NEAR(0.0, spprfs_obj->diff_xerr, spprfs_obj->threshold);
    EXPECT_NEAR(0.0, spprfs_obj->diff_berr, spprfs_obj->threshold);
    EXPECT_NEAR(0.0, spprfs_obj->diff_ferr, spprfs_obj->threshold);
}

/* Begin pprfs_double_parameters  class definition */
class pprfs_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr, diff_s;
      void *hModule, *dModule;
      float threshold;

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
      pprfs_double_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~pprfs_double_parameters (); 
};  /* end of pprfs_double_parameters  class definition */


/* Constructor pprfs_double_parameters definition */
pprfs_double_parameters:: pprfs_double_parameters ( int matrix_layout_i, 
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
   printf(" \n pprfs double:  n: %d, fact: %c uplo: %c   \
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
       pprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, a_buf_size);
    memcpy(af, a, (a_buf_size*sizeof(double)));
    memcpy(afref, a, (a_buf_size*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(double)));
    memcpy(xref, b, (b_bufsize*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

pprfs_double_parameters:: ~pprfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pprfs_double_parameters object: destructor invoked. \n");
#endif
   pprfs_free();
}

//  Test fixture class definition
class dpprfs_test  : public  ::testing::Test {
public:
   pprfs_double_parameters  *dpprfs_obj;
   void SetUp();  
   void TearDown () { delete dpprfs_obj; }
};


void dpprfs_test::SetUp(){

    /* LAPACKE DPPRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_dpprfs) (int matrix_layout, 
		char uplo, lapack_int n, lapack_int nrhs, const double* ap,
		const double* afp, const double* b, lapack_int ldb,
		double* x, lapack_int ldx, double* ferr, double* berr );

    Fptr_NL_LAPACKE_dpprfs DPPRFS;

    /* LAPACKE DPPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dpptrs) ( int matrix_layout, 
					char uplo, lapack_int n, lapack_int nrhs,
              const double * a, double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dpptrs DPPTRS;

     /* LAPACKE DPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpptrf) ( int matrix_layout, 
           char uplo , lapack_int n , double * a );

    Fptr_NL_LAPACKE_dpptrf DPPTRF;

    dpprfs_obj = new pprfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);
    dpprfs_obj->threshold = lin_solver_paramslist[idx].solver_threhold;

    idx = Circular_Increment_Index(idx);

    dpprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DPPRFS = (Fptr_NL_LAPACKE_dpprfs)dlsym(dpprfs_obj->hModule, "LAPACKE_dpprfs");
    ASSERT_TRUE(DPPRFS != NULL) << "failed to ppt the Netlib LAPACKE_dpprfs symbol";

    DPPTRS = (Fptr_NL_LAPACKE_dpptrs)dlsym(dpprfs_obj->hModule, "LAPACKE_dpptrs");
    ASSERT_TRUE(DPPTRS != NULL) << "failed to get the Netlib LAPACKE_dpptrs symbol";

    DPPTRF = (Fptr_NL_LAPACKE_dpptrf)dlsym(dpprfs_obj->hModule,"LAPACKE_dpptrf");
    ASSERT_TRUE(DPPTRF != NULL) << "failed to get the Netlib LAPACKE_dpptrf symbol";
    
    /* invoke the pptrf API to compute the factorized A.  */
        dpprfs_obj->inforef = DPPTRF( dpprfs_obj->matrix_layout,
                                      dpprfs_obj->uplo, dpprfs_obj->n,
                                      dpprfs_obj->afref);
                               
        dpprfs_obj->info = LAPACKE_dpptrf( dpprfs_obj->matrix_layout,
                                           dpprfs_obj->uplo, dpprfs_obj->n,
                                           dpprfs_obj->af);

    dpprfs_obj->inforef = DPPTRS( dpprfs_obj->matrix_layout, 
          dpprfs_obj->uplo, dpprfs_obj->n,
         dpprfs_obj->nrhs, (const double *)dpprfs_obj->afref, 
                          dpprfs_obj->xref, dpprfs_obj->ldx);

    dpprfs_obj->info = LAPACKE_dpptrs( dpprfs_obj->matrix_layout, 
               dpprfs_obj->uplo, dpprfs_obj->n,
                 dpprfs_obj->nrhs, (const double *)dpprfs_obj->af, 
                                 dpprfs_obj->x, dpprfs_obj->ldx );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
	
    dpprfs_obj->inforef = DPPRFS( dpprfs_obj->matrix_layout,
                                  dpprfs_obj->uplo, dpprfs_obj->n,
                                  dpprfs_obj->nrhs,
                                  dpprfs_obj->aref, 
                                  dpprfs_obj->afref,
                                  dpprfs_obj->bref, dpprfs_obj->ldb,
                                  dpprfs_obj->xref, dpprfs_obj->ldx,
                                  dpprfs_obj->ferrref,
                                  dpprfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dpprfs_obj->info = LAPACKE_dpprfs( dpprfs_obj->matrix_layout,
                                  dpprfs_obj->uplo, dpprfs_obj->n,
                                  dpprfs_obj->nrhs,
                                  dpprfs_obj->a, 
                                  dpprfs_obj->af,
                                  dpprfs_obj->b, dpprfs_obj->ldb,
                                  dpprfs_obj->x, dpprfs_obj->ldx,
                                  dpprfs_obj->ferr,
                                  dpprfs_obj->berr);

    if( dpprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dpprfs is wrong\n", dpprfs_obj->info );
    }
    if( dpprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpprfs is wrong\n", 
        dpprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dpprfs_obj->diff_xerr =  computeDiff_d( dpprfs_obj->x_bufsize, 
                dpprfs_obj->x, dpprfs_obj->xref );

    dpprfs_obj->diff_berr =  computeDiff_d( dpprfs_obj->nrhs, 
                dpprfs_obj->berr, dpprfs_obj->berrref );
                
    dpprfs_obj->diff_ferr =  computeDiff_d( dpprfs_obj->nrhs, 
                dpprfs_obj->ferr, dpprfs_obj->ferrref );
    
}

TEST_F(dpprfs_test, dpprfs1) {
    EXPECT_NEAR(0.0, dpprfs_obj->diff_xerr, dpprfs_obj->threshold);
    EXPECT_NEAR(0.0, dpprfs_obj->diff_berr, dpprfs_obj->threshold);
    EXPECT_NEAR(0.0, dpprfs_obj->diff_ferr, dpprfs_obj->threshold);
}

TEST_F(dpprfs_test, dpprfs2) {
    EXPECT_NEAR(0.0, dpprfs_obj->diff_xerr, dpprfs_obj->threshold);
    EXPECT_NEAR(0.0, dpprfs_obj->diff_berr, dpprfs_obj->threshold);
    EXPECT_NEAR(0.0, dpprfs_obj->diff_ferr, dpprfs_obj->threshold);
}

TEST_F(dpprfs_test, dpprfs3) {
    EXPECT_NEAR(0.0, dpprfs_obj->diff_xerr, dpprfs_obj->threshold);
    EXPECT_NEAR(0.0, dpprfs_obj->diff_berr, dpprfs_obj->threshold);
    EXPECT_NEAR(0.0, dpprfs_obj->diff_ferr, dpprfs_obj->threshold);
}

TEST_F(dpprfs_test, dpprfs4) {
    EXPECT_NEAR(0.0, dpprfs_obj->diff_xerr, dpprfs_obj->threshold);
    EXPECT_NEAR(0.0, dpprfs_obj->diff_berr, dpprfs_obj->threshold);
    EXPECT_NEAR(0.0, dpprfs_obj->diff_ferr, dpprfs_obj->threshold);
}

/* Begin pprfs_scomplex_parameters  class definition */
class pprfs_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr, diff_s;
      void *hModule, *dModule;
      float threshold;

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
      lapack_complex_float* s, *sref; // the row scale factors for A.
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      lapack_complex_float rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pprfs_scomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~pprfs_scomplex_parameters (); 
};  /* end of pprfs_scomplex_parameters  class definition */


/* Constructor pprfs_scomplex_parameters definition */
pprfs_scomplex_parameters:: pprfs_scomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n pprfs lapack_complex_float:  n: %d, fact: %c uplo: %c   \
ldb: %d nrhs: %d   ldx: %d \n",  n, fact, uplo,  
                                          ldb, nrhs,  ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_buf_size);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &af, &afref, a_buf_size);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &s, &sref, n);
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
       pprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, a_buf_size);
    memcpy(af, a, (a_buf_size*sizeof(lapack_complex_float)));
    memcpy(afref, a, (a_buf_size*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, (b_bufsize*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

pprfs_scomplex_parameters:: ~pprfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pprfs_scomplex_parameters object: destructor invoked. \n");
#endif
   pprfs_free();
}

//  Test fixture class definition
class cpprfs_test  : public  ::testing::Test {
public:
   pprfs_scomplex_parameters  *cpprfs_obj;
   void SetUp();  
   void TearDown () { delete cpprfs_obj; }
};


void cpprfs_test::SetUp(){

    /* LAPACKE CPPRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_cpprfs) (int matrix_layout, char uplo, 
	lapack_int n, lapack_int nrhs, const lapack_complex_float* ap, 
	const lapack_complex_float* afp, const lapack_complex_float* b, 
	lapack_int ldb, lapack_complex_float* x, lapack_int ldx, 
	float* ferr, float* berr );

    Fptr_NL_LAPACKE_cpprfs CPPRFS;

    /* LAPACKE CPPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_cpptrs) ( int matrix_layout , 
	char uplo , lapack_int n , lapack_int nrhs , 
	const lapack_complex_float * ap , 
	lapack_complex_float * b , lapack_int ldb  );

    Fptr_NL_LAPACKE_cpptrs CPPTRS;

     /* LAPACKE CPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpptrf) ( int matrix_layout , 
		char uplo , lapack_int n , lapack_complex_float * ap );

    Fptr_NL_LAPACKE_cpptrf CPPTRF;

    cpprfs_obj = new pprfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);
    cpprfs_obj->threshold = lin_solver_paramslist[idx].solver_threhold;

    idx = Circular_Increment_Index(idx);

    cpprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CPPRFS = (Fptr_NL_LAPACKE_cpprfs)dlsym(cpprfs_obj->hModule, "LAPACKE_cpprfs");
    ASSERT_TRUE(CPPRFS != NULL) << "failed to ppt the Netlib LAPACKE_cpprfs symbol";

    CPPTRS = (Fptr_NL_LAPACKE_cpptrs)dlsym(cpprfs_obj->hModule, "LAPACKE_cpptrs");
    ASSERT_TRUE(CPPTRS != NULL) << "failed to get the Netlib LAPACKE_cpptrs symbol";

    CPPTRF = (Fptr_NL_LAPACKE_cpptrf)dlsym(cpprfs_obj->hModule,"LAPACKE_cpptrf");
    ASSERT_TRUE(CPPTRF != NULL) << "failed to get the Netlib LAPACKE_cpptrf symbol";
    
    /* invoke the pptrf API to compute the factorized A.  */
        cpprfs_obj->inforef = CPPTRF( cpprfs_obj->matrix_layout,
                                      cpprfs_obj->uplo, cpprfs_obj->n,
                                      cpprfs_obj->afref);
                               
        cpprfs_obj->info = LAPACKE_cpptrf( cpprfs_obj->matrix_layout,
                                           cpprfs_obj->uplo, cpprfs_obj->n,
                                           cpprfs_obj->af);

    cpprfs_obj->inforef = CPPTRS( cpprfs_obj->matrix_layout, 
          cpprfs_obj->uplo, cpprfs_obj->n,
         cpprfs_obj->nrhs, (const lapack_complex_float *)cpprfs_obj->afref, 
                          cpprfs_obj->xref, cpprfs_obj->ldx);

    cpprfs_obj->info = LAPACKE_cpptrs( cpprfs_obj->matrix_layout, 
               cpprfs_obj->uplo, cpprfs_obj->n,
                 cpprfs_obj->nrhs, (const lapack_complex_float *)cpprfs_obj->af, 
                                 cpprfs_obj->x, cpprfs_obj->ldx );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
	
    cpprfs_obj->inforef = CPPRFS( cpprfs_obj->matrix_layout,
                                  cpprfs_obj->uplo, cpprfs_obj->n,
                                  cpprfs_obj->nrhs,
                                  cpprfs_obj->aref, 
                                  cpprfs_obj->afref,
                                  cpprfs_obj->bref, cpprfs_obj->ldb,
                                  cpprfs_obj->xref, cpprfs_obj->ldx,
                                  cpprfs_obj->ferrref,
                                  cpprfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cpprfs_obj->info = LAPACKE_cpprfs( cpprfs_obj->matrix_layout,
                                  cpprfs_obj->uplo, cpprfs_obj->n,
                                  cpprfs_obj->nrhs,
                                  cpprfs_obj->a, 
                                  cpprfs_obj->af,
                                  cpprfs_obj->b, cpprfs_obj->ldb,
                                  cpprfs_obj->x, cpprfs_obj->ldx,
                                  cpprfs_obj->ferr,
                                  cpprfs_obj->berr);

    if( cpprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cpprfs is wrong\n", cpprfs_obj->info );
    }
    if( cpprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpprfs is wrong\n", 
        cpprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cpprfs_obj->diff_xerr =  computeDiff_c( cpprfs_obj->x_bufsize, 
                cpprfs_obj->x, cpprfs_obj->xref );

    cpprfs_obj->diff_berr =  computeDiff_s( cpprfs_obj->nrhs, 
                cpprfs_obj->berr, cpprfs_obj->berrref );
                
    cpprfs_obj->diff_ferr =  computeDiff_s( cpprfs_obj->nrhs, 
                cpprfs_obj->ferr, cpprfs_obj->ferrref );
    
}

TEST_F(cpprfs_test, cpprfs1) {
    EXPECT_NEAR(0.0, cpprfs_obj->diff_xerr, cpprfs_obj->threshold);
    EXPECT_NEAR(0.0, cpprfs_obj->diff_berr, cpprfs_obj->threshold);
    EXPECT_NEAR(0.0, cpprfs_obj->diff_ferr, cpprfs_obj->threshold);
}

TEST_F(cpprfs_test, cpprfs2) {
    EXPECT_NEAR(0.0, cpprfs_obj->diff_xerr, cpprfs_obj->threshold);
    EXPECT_NEAR(0.0, cpprfs_obj->diff_berr, cpprfs_obj->threshold);
    EXPECT_NEAR(0.0, cpprfs_obj->diff_ferr, cpprfs_obj->threshold);
}

TEST_F(cpprfs_test, cpprfs3) {
    EXPECT_NEAR(0.0, cpprfs_obj->diff_xerr, cpprfs_obj->threshold);
    EXPECT_NEAR(0.0, cpprfs_obj->diff_berr, cpprfs_obj->threshold);
    EXPECT_NEAR(0.0, cpprfs_obj->diff_ferr, cpprfs_obj->threshold);
}

TEST_F(cpprfs_test, cpprfs4) {
    EXPECT_NEAR(0.0, cpprfs_obj->diff_xerr, cpprfs_obj->threshold);
    EXPECT_NEAR(0.0, cpprfs_obj->diff_berr, cpprfs_obj->threshold);
    EXPECT_NEAR(0.0, cpprfs_obj->diff_ferr, cpprfs_obj->threshold);
}

/* Begin pprfs_dcomplex_parameters  class definition */
class pprfs_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr, diff_s;
      void *hModule, *dModule;
      float threshold;

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
      lapack_complex_double* s, *sref; // the row scale factors for A.
      
      /* Output parameters */
      char   equed, equedref; //  Must be 'N', 'R', 'C', or 'B'.
      lapack_complex_double rcond, rcondref; // estimate of the reciprocal condition number of the matrix A.
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pprfs_dcomplex_parameters ( int matrix_layout_i, char fact_i, char uplo_i,
                               char equed_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~pprfs_dcomplex_parameters (); 
};  /* end of pprfs_dcomplex_parameters  class definition */


/* Constructor pprfs_dcomplex_parameters definition */
pprfs_dcomplex_parameters:: pprfs_dcomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n pprfs lapack_complex_double:  n: %d, fact: %c uplo: %c   \
ldb: %d nrhs: %d   ldx: %d \n",  n, fact, uplo,  
                                          ldb, nrhs,  ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_buf_size);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &af, &afref, a_buf_size);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &s, &sref, n);
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
       pprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, a_buf_size);
    memcpy(af, a, (a_buf_size*sizeof(lapack_complex_double)));
    memcpy(afref, a, (a_buf_size*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, (b_bufsize*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(s, sref, n, 0.0);

    
   } /* end of Constructor  */

pprfs_dcomplex_parameters:: ~pprfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pprfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   pprfs_free();
}

//  Test fixture class definition
class zpprfs_test  : public  ::testing::Test {
public:
   pprfs_dcomplex_parameters  *zpprfs_obj;
   void SetUp();  
   void TearDown () { delete zpprfs_obj; }
};


void zpprfs_test::SetUp(){

    /* LAPACKE ZPPRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_zpprfs) (int matrix_layout, char uplo, 
	lapack_int n, lapack_int nrhs, const lapack_complex_double* ap, 
	const lapack_complex_double* afp, const lapack_complex_double* b, 
	lapack_int ldb, lapack_complex_double* x, lapack_int ldx, 
	double* ferr, double* berr );

    Fptr_NL_LAPACKE_zpprfs ZPPRFS;

    /* LAPACKE ZPPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zpptrs) ( int matrix_layout , 
	char uplo , lapack_int n , lapack_int nrhs , 
	const lapack_complex_double * ap , 
	lapack_complex_double * b , lapack_int ldb  );

    Fptr_NL_LAPACKE_zpptrs ZPPTRS;

     /* LAPACKE ZPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpptrf) ( int matrix_layout , 
		char uplo , lapack_int n , lapack_complex_double * ap );

    Fptr_NL_LAPACKE_zpptrf ZPPTRF;

    zpprfs_obj = new pprfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].fact,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].equed,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);
    zpprfs_obj->threshold = lin_solver_paramslist[idx].solver_threhold;

    idx = Circular_Increment_Index(idx);

    zpprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZPPRFS = (Fptr_NL_LAPACKE_zpprfs)dlsym(zpprfs_obj->hModule, "LAPACKE_zpprfs");
    ASSERT_TRUE(ZPPRFS != NULL) << "failed to ppt the Netlib LAPACKE_zpprfs symbol";

    ZPPTRS = (Fptr_NL_LAPACKE_zpptrs)dlsym(zpprfs_obj->hModule, "LAPACKE_zpptrs");
    ASSERT_TRUE(ZPPTRS != NULL) << "failed to get the Netlib LAPACKE_zpptrs symbol";

    ZPPTRF = (Fptr_NL_LAPACKE_zpptrf)dlsym(zpprfs_obj->hModule,"LAPACKE_zpptrf");
    ASSERT_TRUE(ZPPTRF != NULL) << "failed to get the Netlib LAPACKE_zpptrf symbol";
    
    /* invoke the pptrf API to compute the factorized A.  */
        zpprfs_obj->inforef = ZPPTRF( zpprfs_obj->matrix_layout,
                                      zpprfs_obj->uplo, zpprfs_obj->n,
                                      zpprfs_obj->afref);
                               
        zpprfs_obj->info = LAPACKE_zpptrf( zpprfs_obj->matrix_layout,
                                           zpprfs_obj->uplo, zpprfs_obj->n,
                                           zpprfs_obj->af);

    zpprfs_obj->inforef = ZPPTRS( zpprfs_obj->matrix_layout, 
          zpprfs_obj->uplo, zpprfs_obj->n,
         zpprfs_obj->nrhs, (const lapack_complex_double *)zpprfs_obj->afref, 
                          zpprfs_obj->xref, zpprfs_obj->ldx);

    zpprfs_obj->info = LAPACKE_zpptrs( zpprfs_obj->matrix_layout, 
               zpprfs_obj->uplo, zpprfs_obj->n,
                 zpprfs_obj->nrhs, (const lapack_complex_double *)zpprfs_obj->af, 
                                 zpprfs_obj->x, zpprfs_obj->ldx );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */
	
    zpprfs_obj->inforef = ZPPRFS( zpprfs_obj->matrix_layout,
                                  zpprfs_obj->uplo, zpprfs_obj->n,
                                  zpprfs_obj->nrhs,
                                  zpprfs_obj->aref, 
                                  zpprfs_obj->afref,
                                  zpprfs_obj->bref, zpprfs_obj->ldb,
                                  zpprfs_obj->xref, zpprfs_obj->ldx,
                                  zpprfs_obj->ferrref,
                                  zpprfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zpprfs_obj->info = LAPACKE_zpprfs( zpprfs_obj->matrix_layout,
                                  zpprfs_obj->uplo, zpprfs_obj->n,
                                  zpprfs_obj->nrhs,
                                  zpprfs_obj->a, 
                                  zpprfs_obj->af,
                                  zpprfs_obj->b, zpprfs_obj->ldb,
                                  zpprfs_obj->x, zpprfs_obj->ldx,
                                  zpprfs_obj->ferr,
                                  zpprfs_obj->berr);

    if( zpprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zpprfs is wrong\n", zpprfs_obj->info );
    }
    if( zpprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpprfs is wrong\n", 
        zpprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zpprfs_obj->diff_xerr =  computeDiff_z( zpprfs_obj->x_bufsize, 
                zpprfs_obj->x, zpprfs_obj->xref );

    zpprfs_obj->diff_berr =  computeDiff_d( zpprfs_obj->nrhs, 
                zpprfs_obj->berr, zpprfs_obj->berrref );
                
    zpprfs_obj->diff_ferr =  computeDiff_d( zpprfs_obj->nrhs, 
                zpprfs_obj->ferr, zpprfs_obj->ferrref );
    
}

TEST_F(zpprfs_test, zpprfs1) {
    EXPECT_NEAR(0.0, zpprfs_obj->diff_xerr, zpprfs_obj->threshold);
    EXPECT_NEAR(0.0, zpprfs_obj->diff_berr, zpprfs_obj->threshold);
    EXPECT_NEAR(0.0, zpprfs_obj->diff_ferr, zpprfs_obj->threshold);
}

TEST_F(zpprfs_test, zpprfs2) {
    EXPECT_NEAR(0.0, zpprfs_obj->diff_xerr, zpprfs_obj->threshold);
    EXPECT_NEAR(0.0, zpprfs_obj->diff_berr, zpprfs_obj->threshold);
    EXPECT_NEAR(0.0, zpprfs_obj->diff_ferr, zpprfs_obj->threshold);
}

TEST_F(zpprfs_test, zpprfs3) {
    EXPECT_NEAR(0.0, zpprfs_obj->diff_xerr, zpprfs_obj->threshold);
    EXPECT_NEAR(0.0, zpprfs_obj->diff_berr, zpprfs_obj->threshold);
    EXPECT_NEAR(0.0, zpprfs_obj->diff_ferr, zpprfs_obj->threshold);
}

TEST_F(zpprfs_test, zpprfs4) {
    EXPECT_NEAR(0.0, zpprfs_obj->diff_xerr, zpprfs_obj->threshold);
    EXPECT_NEAR(0.0, zpprfs_obj->diff_berr, zpprfs_obj->threshold);
    EXPECT_NEAR(0.0, zpprfs_obj->diff_ferr, zpprfs_obj->threshold);
}