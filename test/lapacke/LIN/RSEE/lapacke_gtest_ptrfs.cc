#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define ptrfs_free() \
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

/* Begin ptrfs_float_parameters  class definition */
class ptrfs_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR

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
      float* x, *xref; //  solution matrix X to the system of equations
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ptrfs_float_parameters ( int matrix_layout_i,
                               lapack_int n_i, lapack_int nrhs_i);
          
      ~ptrfs_float_parameters (); 
};  /* end of ptrfs_float_parameters  class definition */


/* Constructor ptrfs_float_parameters definition */
ptrfs_float_parameters:: ptrfs_float_parameters ( int matrix_layout_i, 
                                 lapack_int n_i, lapack_int nrhs_i) {
                                
    matrix_layout = matrix_layout_i;
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
   printf(" \n ptrfs float:  matrix_layout: %d  n: %d, ldb: %d nrhs: %d  \
ldx: %d \n",  matrix_layout, n, ldb, nrhs, ldx);
#endif

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
       ptrfs_free();
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
    memcpy(x, b, ( b_bufsize*sizeof(float)));
    memcpy(xref, b, ( b_bufsize*sizeof(float)));

    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Constructor  */

ptrfs_float_parameters:: ~ptrfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptrfs_float_parameters object: destructor invoked. \n");
#endif
   ptrfs_free();
}

//  Test fixture class definition
class sptrfs_test  : public  ::testing::Test {
public:
   ptrfs_float_parameters  *sptrfs_obj;
   void SetUp();  
   void TearDown () { delete sptrfs_obj; }
};


void sptrfs_test::SetUp(){

    /* LAPACKE SPTRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_sptrfs) (int matrix_layout,
		lapack_int n, lapack_int nrhs, const float* d,
		const float* e, const float* df, const float* ef,
		const float* b, lapack_int ldb, float* x,
		lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_sptrfs SPTRFS;

     /* LAPACKE SPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spttrf) ( lapack_int n, float* d, float* e );
    Fptr_NL_LAPACKE_spttrf SPTTRF;
	
	typedef int (*Fptr_NL_LAPACKE_spttrs)( int matrix_layout,
		lapack_int n, lapack_int nrhs, const float* d,
		const float* e, float* b, lapack_int ldb );
	
	Fptr_NL_LAPACKE_spttrs  SPTTRS;

    sptrfs_obj = new ptrfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    sptrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sptrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sptrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sptrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SPTRFS = (Fptr_NL_LAPACKE_sptrfs)dlsym(sptrfs_obj->hModule, "LAPACKE_sptrfs");
    ASSERT_TRUE(SPTRFS != NULL) << "failed to ptt the Netlib LAPACKE_sptrfs symbol";

	SPTTRF = (Fptr_NL_LAPACKE_spttrf)dlsym(sptrfs_obj->hModule,"LAPACKE_spttrf");
	ASSERT_TRUE(SPTTRF != NULL) << "failed to ptt the Netlib LAPACKE_spttrf symbol";
	
	SPTTRS = (Fptr_NL_LAPACKE_spttrs)dlsym(sptrfs_obj->hModule,"LAPACKE_spttrs");
	ASSERT_TRUE(SPTTRS != NULL) << "failed to ptt the Netlib LAPACKE_spttrs symbol";

	sptrfs_obj->inforef = SPTTRF( sptrfs_obj->n,
								  sptrfs_obj->dfref,
								  sptrfs_obj->efref);

	sptrfs_obj->info = LAPACKE_spttrf( sptrfs_obj->n,
									  sptrfs_obj->df,
									  sptrfs_obj->ef);

	sptrfs_obj->inforef = SPTTRS( sptrfs_obj->matrix_layout,
                                  sptrfs_obj->n,
                                  sptrfs_obj->nrhs,
                                  (const float *)sptrfs_obj->dfref,
                                  (const float *)sptrfs_obj->efref,
                                  sptrfs_obj->xref,
								  sptrfs_obj->ldx);

	sptrfs_obj->info = SPTTRS( sptrfs_obj->matrix_layout,
                                  sptrfs_obj->n,
                                  sptrfs_obj->nrhs,
                                  (const float *)sptrfs_obj->df,
                                  (const float *)sptrfs_obj->ef,
                                  sptrfs_obj->x,
								  sptrfs_obj->ldx);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    sptrfs_obj->inforef = SPTRFS( sptrfs_obj->matrix_layout,
                                  sptrfs_obj->n,
                                  sptrfs_obj->nrhs,
                                  (const float *)sptrfs_obj->dref,
                                  (const float *)sptrfs_obj->eref, 
                                  (const float *)sptrfs_obj->dfref,
								  (const float *)sptrfs_obj->efref,
                                  (const float *)sptrfs_obj->bref, 
                                  sptrfs_obj->ldb,
                                  sptrfs_obj->xref, sptrfs_obj->ldx,
                                  sptrfs_obj->ferrref,
                                  sptrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    sptrfs_obj->info = LAPACKE_sptrfs( sptrfs_obj->matrix_layout,
                                  sptrfs_obj->n,
                                  sptrfs_obj->nrhs,
                                  (const float *)sptrfs_obj->d,
                                  (const float *)sptrfs_obj->e, 
                                  (const float *)sptrfs_obj->df,
								  (const float *)sptrfs_obj->ef,
                                  (const float *)sptrfs_obj->b,
                                  sptrfs_obj->ldb,
                                  sptrfs_obj->x, sptrfs_obj->ldx,
                                  sptrfs_obj->ferr,
                                  sptrfs_obj->berr);

    if( sptrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sptrfs is wrong\n", sptrfs_obj->info );
    }
    if( sptrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sptrfs is wrong\n", 
        sptrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sptrfs_obj->diff_xerr =  computeDiff_s( sptrfs_obj->x_bufsize, 
                sptrfs_obj->x, sptrfs_obj->xref );

    sptrfs_obj->diff_berr =  computeDiff_s( sptrfs_obj->nrhs, 
                sptrfs_obj->berr, sptrfs_obj->berrref );
            
    sptrfs_obj->diff_ferr =  computeDiff_s( sptrfs_obj->nrhs, 
                sptrfs_obj->ferr, sptrfs_obj->ferrref );
}

TEST_F(sptrfs_test, sptrfs1) {
    EXPECT_NEAR(0.0, sptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(sptrfs_test, sptrfs2) {
    EXPECT_NEAR(0.0, sptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(sptrfs_test, sptrfs3) {
    EXPECT_NEAR(0.0, sptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(sptrfs_test, sptrfs4) {
    EXPECT_NEAR(0.0, sptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, sptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin ptrfs_double_parameters  class definition */
class ptrfs_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR

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
      double* x, *xref; //  solution matrix X to the system of equations
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ptrfs_double_parameters ( int matrix_layout_i,
                               lapack_int n_i, lapack_int nrhs_i);
          
      ~ptrfs_double_parameters (); 
};  /* end of ptrfs_double_parameters  class definition */


/* Constructor ptrfs_double_parameters definition */
ptrfs_double_parameters:: ptrfs_double_parameters ( int matrix_layout_i, 
                                 lapack_int n_i, lapack_int nrhs_i) {
                                
    matrix_layout = matrix_layout_i;
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
   printf(" \n ptrfs double:  matrix_layout: %d  n: %d, ldb: %d nrhs: %d  \
ldx: %d \n",  matrix_layout, n, ldb, nrhs, ldx);
#endif

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
       ptrfs_free();
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
    memcpy(x, b, ( b_bufsize*sizeof(double)));
    memcpy(xref, b, ( b_bufsize*sizeof(double)));

    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Constructor  */

ptrfs_double_parameters:: ~ptrfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptrfs_double_parameters object: destructor invoked. \n");
#endif
   ptrfs_free();
}

//  Test fixture class definition
class dptrfs_test  : public  ::testing::Test {
public:
   ptrfs_double_parameters  *dptrfs_obj;
   void SetUp();  
   void TearDown () { delete dptrfs_obj; }
};


void dptrfs_test::SetUp(){

    /* LAPACKE DPTRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_dptrfs) (int matrix_layout,
		lapack_int n, lapack_int nrhs, const double* d,
		const double* e, const double* df, const double* ef,
		const double* b, lapack_int ldb, double* x,
		lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_dptrfs DPTRFS;

     /* LAPACKE DPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpttrf) ( lapack_int n, double* d, double* e );
    Fptr_NL_LAPACKE_dpttrf DPTTRF;
	
	typedef int (*Fptr_NL_LAPACKE_dpttrs)( int matrix_layout,
		lapack_int n, lapack_int nrhs, const double* d,
		const double* e, double* b, lapack_int ldb );
	
	Fptr_NL_LAPACKE_dpttrs  DPTTRS;

    dptrfs_obj = new ptrfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    dptrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dptrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dptrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dptrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DPTRFS = (Fptr_NL_LAPACKE_dptrfs)dlsym(dptrfs_obj->hModule, "LAPACKE_dptrfs");
    ASSERT_TRUE(DPTRFS != NULL) << "failed to ptt the Netlib LAPACKE_dptrfs symbol";

	DPTTRF = (Fptr_NL_LAPACKE_dpttrf)dlsym(dptrfs_obj->hModule,"LAPACKE_dpttrf");
	ASSERT_TRUE(DPTTRF != NULL) << "failed to ptt the Netlib LAPACKE_dpttrf symbol";
	
	DPTTRS = (Fptr_NL_LAPACKE_dpttrs)dlsym(dptrfs_obj->hModule,"LAPACKE_dpttrs");
	ASSERT_TRUE(DPTTRS != NULL) << "failed to ptt the Netlib LAPACKE_dpttrs symbol";

	dptrfs_obj->inforef = DPTTRF( dptrfs_obj->n,
								  dptrfs_obj->dfref,
								  dptrfs_obj->efref);

	dptrfs_obj->info = LAPACKE_dpttrf( dptrfs_obj->n,
									  dptrfs_obj->df,
									  dptrfs_obj->ef);

	dptrfs_obj->inforef = DPTTRS( dptrfs_obj->matrix_layout,
                                  dptrfs_obj->n,
                                  dptrfs_obj->nrhs,
                                  (const double *)dptrfs_obj->dfref,
                                  (const double *)dptrfs_obj->efref,
                                  dptrfs_obj->xref,
								  dptrfs_obj->ldx);

	dptrfs_obj->info = DPTTRS( dptrfs_obj->matrix_layout,
                                  dptrfs_obj->n,
                                  dptrfs_obj->nrhs,
                                  (const double *)dptrfs_obj->df,
                                  (const double *)dptrfs_obj->ef,
                                  dptrfs_obj->x,
								  dptrfs_obj->ldx);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dptrfs_obj->inforef = DPTRFS( dptrfs_obj->matrix_layout,
                                  dptrfs_obj->n,
                                  dptrfs_obj->nrhs,
                                  (const double *)dptrfs_obj->dref,
                                  (const double *)dptrfs_obj->eref, 
                                  (const double *)dptrfs_obj->dfref,
								  (const double *)dptrfs_obj->efref,
                                  (const double *)dptrfs_obj->bref, 
                                  dptrfs_obj->ldb,
                                  dptrfs_obj->xref, dptrfs_obj->ldx,
                                  dptrfs_obj->ferrref,
                                  dptrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dptrfs_obj->info = LAPACKE_dptrfs( dptrfs_obj->matrix_layout,
                                  dptrfs_obj->n,
                                  dptrfs_obj->nrhs,
                                  (const double *)dptrfs_obj->d,
                                  (const double *)dptrfs_obj->e, 
                                  (const double *)dptrfs_obj->df,
								  (const double *)dptrfs_obj->ef,
                                  (const double *)dptrfs_obj->b,
                                  dptrfs_obj->ldb,
                                  dptrfs_obj->x, dptrfs_obj->ldx,
                                  dptrfs_obj->ferr,
                                  dptrfs_obj->berr);

    if( dptrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dptrfs is wrong\n", dptrfs_obj->info );
    }
    if( dptrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dptrfs is wrong\n", 
        dptrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dptrfs_obj->diff_xerr =  computeDiff_d( dptrfs_obj->x_bufsize, 
                dptrfs_obj->x, dptrfs_obj->xref );

    dptrfs_obj->diff_berr =  computeDiff_d( dptrfs_obj->nrhs, 
                dptrfs_obj->berr, dptrfs_obj->berrref );
            
    dptrfs_obj->diff_ferr =  computeDiff_d( dptrfs_obj->nrhs, 
                dptrfs_obj->ferr, dptrfs_obj->ferrref );
}

TEST_F(dptrfs_test, dptrfs1) {
    EXPECT_NEAR(0.0, dptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dptrfs_test, dptrfs2) {
    EXPECT_NEAR(0.0, dptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dptrfs_test, dptrfs3) {
    EXPECT_NEAR(0.0, dptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dptrfs_test, dptrfs4) {
    EXPECT_NEAR(0.0, dptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin ptrfs_scomplex_parameters  class definition */
class ptrfs_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // U or L
      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      float *d, *dref; // n diagonal elements of the tridiagonal matrix A
      float *df, *dfref; // n diagonal elements of the diagonal matrix D
      lapack_complex_float *e, *eref; // (n - 1) subdiagonal elements of the tridiagonal matrix A.
      lapack_complex_float* ef, *efref; // (n - 1) subdiagonal elements of the unit bidiagonal factor
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.
  
      /* Output parameters */
      lapack_complex_float* x, *xref; //  solution matrix X to the system of equations
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ptrfs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
          
      ~ptrfs_scomplex_parameters (); 
};  /* end of ptrfs_scomplex_parameters  class definition */


/* Constructor ptrfs_scomplex_parameters definition */
ptrfs_scomplex_parameters:: ptrfs_scomplex_parameters ( int matrix_layout_i, 
                          char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                
    matrix_layout = matrix_layout_i;
    n = n_i;
    nrhs = nrhs_i;
	uplo = uplo_i;

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
   printf(" \n ptrfs lapack_complex_float:  matrix_layout: %d  n: %d, ldb: %d nrhs: %d  \
ldx: %d \n",  matrix_layout, n, ldb, nrhs, ldx);
#endif

    diff_berr = 0.0;
    diff_ferr = 0.0;
    diff_xerr = 0.0;
    hModule = NULL;
    dModule = NULL;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_float_buffer_pair( &df, &dfref, n);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e,  &eref,  n-1);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &ef, &efref, n-1);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
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
       ptrfs_free();
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
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_float)));

    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Constructor  */

ptrfs_scomplex_parameters:: ~ptrfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptrfs_scomplex_parameters object: destructor invoked. \n");
#endif
   ptrfs_free();
}

//  Test fixture class definition
class cptrfs_test  : public  ::testing::Test {
public:
   ptrfs_scomplex_parameters  *cptrfs_obj;
   void SetUp();  
   void TearDown () { delete cptrfs_obj; }
};


void cptrfs_test::SetUp(){

    /* LAPACKE CPTRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_cptrfs) (int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const float* d,
		const lapack_complex_float* e, const float* df,
		const lapack_complex_float* ef, const lapack_complex_float* b,
		lapack_int ldb, lapack_complex_float* x, lapack_int ldx,
		float* ferr, float* berr);

    Fptr_NL_LAPACKE_cptrfs CPTRFS;

     /* LAPACKE CPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpttrf) ( lapack_int n, 
		float* d, lapack_complex_float* e );

    Fptr_NL_LAPACKE_cpttrf CPTTRF;
	
	typedef int (*Fptr_NL_LAPACKE_cpttrs)( int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const float* d,
		const lapack_complex_float* e, lapack_complex_float* b,
		lapack_int ldb );
	
	Fptr_NL_LAPACKE_cpttrs  CPTTRS;

    cptrfs_obj = new ptrfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    cptrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cptrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cptrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cptrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CPTRFS = (Fptr_NL_LAPACKE_cptrfs)dlsym(cptrfs_obj->hModule, "LAPACKE_cptrfs");
    ASSERT_TRUE(CPTRFS != NULL) << "failed to ptt the Netlib LAPACKE_cptrfs symbol";

	CPTTRF = (Fptr_NL_LAPACKE_cpttrf)dlsym(cptrfs_obj->hModule,"LAPACKE_cpttrf");
	ASSERT_TRUE(CPTTRF != NULL) << "failed to ptt the Netlib LAPACKE_cpttrf symbol";
	
	CPTTRS = (Fptr_NL_LAPACKE_cpttrs)dlsym(cptrfs_obj->hModule,"LAPACKE_cpttrs");
	ASSERT_TRUE(CPTTRS != NULL) << "failed to ptt the Netlib LAPACKE_cpttrs symbol";

	cptrfs_obj->inforef = CPTTRF( cptrfs_obj->n,
								  cptrfs_obj->dfref,
								  cptrfs_obj->efref);

	cptrfs_obj->info = LAPACKE_cpttrf( cptrfs_obj->n,
									  cptrfs_obj->df,
									  cptrfs_obj->ef);

	cptrfs_obj->inforef = CPTTRS( cptrfs_obj->matrix_layout,
                                  cptrfs_obj->uplo,
                                  cptrfs_obj->n,
                                  cptrfs_obj->nrhs,
                                  (const float *)cptrfs_obj->dfref,
                                  (const lapack_complex_float *)cptrfs_obj->efref,
                                  cptrfs_obj->xref,
								  cptrfs_obj->ldx);

	cptrfs_obj->info = LAPACKE_cpttrs( cptrfs_obj->matrix_layout,
                                  cptrfs_obj->uplo,
                                  cptrfs_obj->n,
                                  cptrfs_obj->nrhs,
                                  (const float *)cptrfs_obj->df,
                                  (const lapack_complex_float *)cptrfs_obj->ef,
                                  cptrfs_obj->x,
								  cptrfs_obj->ldx);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    cptrfs_obj->inforef = CPTRFS( cptrfs_obj->matrix_layout,
                                  cptrfs_obj->uplo,
                                  cptrfs_obj->n,
                                  cptrfs_obj->nrhs,
                                  (const float *)cptrfs_obj->dref,
                                  (const lapack_complex_float *)cptrfs_obj->eref, 
                                  (const float *)cptrfs_obj->dfref,
								  (const lapack_complex_float *)cptrfs_obj->efref,
                                  (const lapack_complex_float *)cptrfs_obj->bref, 
                                  cptrfs_obj->ldb,
                                  cptrfs_obj->xref, cptrfs_obj->ldx,
                                  cptrfs_obj->ferrref,
                                  cptrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cptrfs_obj->info = LAPACKE_cptrfs( cptrfs_obj->matrix_layout,
                                  cptrfs_obj->uplo,
                                  cptrfs_obj->n,
                                  cptrfs_obj->nrhs,
                                  (const float *)cptrfs_obj->d,
                                  (const lapack_complex_float *)cptrfs_obj->e, 
                                  (const float *)cptrfs_obj->df,
								  (const lapack_complex_float *)cptrfs_obj->ef,
                                  (const lapack_complex_float *)cptrfs_obj->b,
                                  cptrfs_obj->ldb,
                                  cptrfs_obj->x, cptrfs_obj->ldx,
                                  cptrfs_obj->ferr,
                                  cptrfs_obj->berr);

    if( cptrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cptrfs is wrong\n", cptrfs_obj->info );
    }
    if( cptrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cptrfs is wrong\n", 
        cptrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cptrfs_obj->diff_xerr =  computeDiff_c( cptrfs_obj->x_bufsize, 
                cptrfs_obj->x, cptrfs_obj->xref );

    cptrfs_obj->diff_berr =  computeDiff_s( cptrfs_obj->nrhs, 
                cptrfs_obj->berr, cptrfs_obj->berrref );
            
    cptrfs_obj->diff_ferr =  computeDiff_s( cptrfs_obj->nrhs, 
                cptrfs_obj->ferr, cptrfs_obj->ferrref );
}

TEST_F(cptrfs_test, cptrfs1) {
    EXPECT_NEAR(0.0, cptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(cptrfs_test, cptrfs2) {
    EXPECT_NEAR(0.0, cptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(cptrfs_test, cptrfs3) {
    EXPECT_NEAR(0.0, cptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(cptrfs_test, cptrfs4) {
    EXPECT_NEAR(0.0, cptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin ptrfs_dcomplex_parameters  class definition */
class ptrfs_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // U or L
      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      double *d, *dref; // n diagonal elements of the tridiagonal matrix A
      double *df, *dfref; // n diagonal elements of the diagonal matrix D
      lapack_complex_double *e, *eref; // (n - 1) subdiagonal elements of the tridiagonal matrix A.
      lapack_complex_double* ef, *efref; // (n - 1) subdiagonal elements of the unit bidiagonal factor
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.
  
      /* Output parameters */
      lapack_complex_double* x, *xref; //  solution matrix X to the system of equations
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ptrfs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
          
      ~ptrfs_dcomplex_parameters (); 
};  /* end of ptrfs_dcomplex_parameters  class definition */


/* Constructor ptrfs_dcomplex_parameters definition */
ptrfs_dcomplex_parameters:: ptrfs_dcomplex_parameters ( int matrix_layout_i, 
                          char uplo_i, lapack_int n_i, lapack_int nrhs_i) {
                                
    matrix_layout = matrix_layout_i;
    n = n_i;
    nrhs = nrhs_i;
	uplo = uplo_i;

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
   printf(" \n ptrfs lapack_complex_double:  matrix_layout: %d  n: %d, ldb: %d nrhs: %d  \
ldx: %d \n",  matrix_layout, n, ldb, nrhs, ldx);
#endif

    diff_berr = 0.0;
    diff_ferr = 0.0;
    diff_xerr = 0.0;
    hModule = NULL;
    dModule = NULL;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_double_buffer_pair( &df, &dfref, n);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e,  &eref,  n-1);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &ef, &efref, n-1);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
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
       ptrfs_free();
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
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_double)));

    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Constructor  */

ptrfs_dcomplex_parameters:: ~ptrfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptrfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   ptrfs_free();
}

//  Test fixture class definition
class zptrfs_test  : public  ::testing::Test {
public:
   ptrfs_dcomplex_parameters  *zptrfs_obj;
   void SetUp();  
   void TearDown () { delete zptrfs_obj; }
};


void zptrfs_test::SetUp(){

    /* LAPACKE ZPTRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_zptrfs) (int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const double* d,
		const lapack_complex_double* e, const double* df,
		const lapack_complex_double* ef, const lapack_complex_double* b,
		lapack_int ldb, lapack_complex_double* x, lapack_int ldx,
		double* ferr, double* berr);

    Fptr_NL_LAPACKE_zptrfs ZPTRFS;

     /* LAPACKE ZPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpttrf) ( lapack_int n, 
		double* d, lapack_complex_double* e );

    Fptr_NL_LAPACKE_zpttrf ZPTTRF;
	
	typedef int (*Fptr_NL_LAPACKE_zpttrs)( int matrix_layout,
		char uplo, lapack_int n, lapack_int nrhs, const double* d,
		const lapack_complex_double* e, lapack_complex_double* b,
		lapack_int ldb );
	
	Fptr_NL_LAPACKE_zpttrs  ZPTTRS;

    zptrfs_obj = new ptrfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    zptrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zptrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zptrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zptrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZPTRFS = (Fptr_NL_LAPACKE_zptrfs)dlsym(zptrfs_obj->hModule, "LAPACKE_zptrfs");
    ASSERT_TRUE(ZPTRFS != NULL) << "failed to ptt the Netlib LAPACKE_zptrfs symbol";

	ZPTTRF = (Fptr_NL_LAPACKE_zpttrf)dlsym(zptrfs_obj->hModule,"LAPACKE_zpttrf");
	ASSERT_TRUE(ZPTTRF != NULL) << "failed to ptt the Netlib LAPACKE_zpttrf symbol";
	
	ZPTTRS = (Fptr_NL_LAPACKE_zpttrs)dlsym(zptrfs_obj->hModule,"LAPACKE_zpttrs");
	ASSERT_TRUE(ZPTTRS != NULL) << "failed to ptt the Netlib LAPACKE_zpttrs symbol";

	zptrfs_obj->inforef = ZPTTRF( zptrfs_obj->n,
								  zptrfs_obj->dfref,
								  zptrfs_obj->efref);

	zptrfs_obj->info = LAPACKE_zpttrf( zptrfs_obj->n,
									  zptrfs_obj->df,
									  zptrfs_obj->ef);

	zptrfs_obj->inforef = ZPTTRS( zptrfs_obj->matrix_layout,
                                  zptrfs_obj->uplo,
                                  zptrfs_obj->n,
                                  zptrfs_obj->nrhs,
                                  (const double *)zptrfs_obj->dfref,
                                  (const lapack_complex_double *)zptrfs_obj->efref,
                                  zptrfs_obj->xref,
								  zptrfs_obj->ldx);

	zptrfs_obj->info = LAPACKE_zpttrs( zptrfs_obj->matrix_layout,
                                  zptrfs_obj->uplo,
                                  zptrfs_obj->n,
                                  zptrfs_obj->nrhs,
                                  (const double *)zptrfs_obj->df,
                                  (const lapack_complex_double *)zptrfs_obj->ef,
                                  zptrfs_obj->x,
								  zptrfs_obj->ldx);

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zptrfs_obj->inforef = ZPTRFS( zptrfs_obj->matrix_layout,
                                  zptrfs_obj->uplo,
                                  zptrfs_obj->n,
                                  zptrfs_obj->nrhs,
                                  (const double *)zptrfs_obj->dref,
                                  (const lapack_complex_double *)zptrfs_obj->eref, 
                                  (const double *)zptrfs_obj->dfref,
								  (const lapack_complex_double *)zptrfs_obj->efref,
                                  (const lapack_complex_double *)zptrfs_obj->bref, 
                                  zptrfs_obj->ldb,
                                  zptrfs_obj->xref, zptrfs_obj->ldx,
                                  zptrfs_obj->ferrref,
                                  zptrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zptrfs_obj->info = LAPACKE_zptrfs( zptrfs_obj->matrix_layout,
                                  zptrfs_obj->uplo,
                                  zptrfs_obj->n,
                                  zptrfs_obj->nrhs,
                                  (const double *)zptrfs_obj->d,
                                  (const lapack_complex_double *)zptrfs_obj->e, 
                                  (const double *)zptrfs_obj->df,
								  (const lapack_complex_double *)zptrfs_obj->ef,
                                  (const lapack_complex_double *)zptrfs_obj->b,
                                  zptrfs_obj->ldb,
                                  zptrfs_obj->x, zptrfs_obj->ldx,
                                  zptrfs_obj->ferr,
                                  zptrfs_obj->berr);

    if( zptrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zptrfs is wrong\n", zptrfs_obj->info );
    }
    if( zptrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zptrfs is wrong\n", 
        zptrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zptrfs_obj->diff_xerr =  computeDiff_z( zptrfs_obj->x_bufsize, 
                zptrfs_obj->x, zptrfs_obj->xref );

    zptrfs_obj->diff_berr =  computeDiff_d( zptrfs_obj->nrhs, 
                zptrfs_obj->berr, zptrfs_obj->berrref );
            
    zptrfs_obj->diff_ferr =  computeDiff_d( zptrfs_obj->nrhs, 
                zptrfs_obj->ferr, zptrfs_obj->ferrref );
}

TEST_F(zptrfs_test, zptrfs1) {
    EXPECT_NEAR(0.0, zptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zptrfs_test, zptrfs2) {
    EXPECT_NEAR(0.0, zptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zptrfs_test, zptrfs3) {
    EXPECT_NEAR(0.0, zptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zptrfs_test, zptrfs4) {
    EXPECT_NEAR(0.0, zptrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zptrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zptrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}


