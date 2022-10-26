#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define tprfs_free() \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
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

/* Begin tprfs_float_parameters  class definition */
class tprfs_float_parameters{
   public:
      int b_bufsize;
      void *hModule, *dModule;
      float diff_berr, diff_ferr, diff_xerr;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx;  //  leading dimension of 'x'
      float *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.
      float *x, *xref; //right-hand sides for the systems of equations.

      /* Output parameters */
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tprfs_float_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~tprfs_float_parameters (); 
};  /* end of tprfs_float_parameters  class definition */


/* Constructor tprfs_float_parameters definition */
tprfs_float_parameters:: tprfs_float_parameters ( int matrix_layout_i, 
                               char uplo_i, char trans_i, char diag_i, 
                                    lapack_int n_i, lapack_int nrhs_i)
{
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = n;
    nrhs = nrhs_i;

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
   printf(" \n tprfs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
#endif

    b_bufsize = n*nrhs;
	
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*(n+1)/2) );
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       tprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, (n*(n+1)/2) );
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(float)));
    memcpy(xref, b, ( b_bufsize*sizeof(float)));

    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Constructor  */

tprfs_float_parameters:: ~tprfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tprfs_float_parameters object: destructor invoked. \n");
#endif
   tprfs_free();
}

//  Test fixture class definition
class stprfs_test  : public  ::testing::Test {
public:
   tprfs_float_parameters  *stprfs_obj;
   void SetUp();  
   void TearDown () { delete stprfs_obj; }
};


void stprfs_test::SetUp(){

    /* LAPACKE STPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_stptrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const float *a, 
                                    float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_stptrs STPTRS;

    typedef int (*Fptr_NL_LAPACKE_stprfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int nrhs, const float* a,
			const float* b, lapack_int ldb, const float* x,
			lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_stprfs STPRFS;

    stprfs_obj = new tprfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    stprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    STPTRS = (Fptr_NL_LAPACKE_stptrs)dlsym(stprfs_obj->hModule, "LAPACKE_stptrs");
    ASSERT_TRUE(STPTRS != NULL) << "failed to get the Netlib LAPACKE_stptrs symbol";

    STPRFS = (Fptr_NL_LAPACKE_stprfs)dlsym(stprfs_obj->hModule, "LAPACKE_stprfs");
    ASSERT_TRUE(STPRFS != NULL) << "failed to get the Netlib LAPACKE_stprfs symbol";

    /* Generate i/ps to tprfs API from 'tptrs' API */
    stprfs_obj->inforef = STPTRS( stprfs_obj->matrix_layout, stprfs_obj->uplo,
                                  stprfs_obj->trans, stprfs_obj->diag, 
                                  stprfs_obj->n, stprfs_obj->nrhs,
                                  (const float *)stprfs_obj->aref, 
                                  stprfs_obj->xref,
                                  stprfs_obj->ldb );

    stprfs_obj->info = LAPACKE_stptrs( stprfs_obj->matrix_layout, stprfs_obj->uplo,
                                  stprfs_obj->trans, stprfs_obj->diag, 
                                  stprfs_obj->n, stprfs_obj->nrhs,
                                  (const float *)stprfs_obj->a, 
                                  stprfs_obj->x,
                                  stprfs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    stprfs_obj->inforef = STPRFS( stprfs_obj->matrix_layout, stprfs_obj->uplo,
                                  stprfs_obj->trans, stprfs_obj->diag, 
                                  stprfs_obj->n, stprfs_obj->nrhs,
                                  (const float *)stprfs_obj->aref, 
                                  stprfs_obj->bref, stprfs_obj->ldb,
								  stprfs_obj->xref, stprfs_obj->ldx,
                                  stprfs_obj->ferrref,
                                  stprfs_obj->berrref
								  );

    stprfs_obj->info = LAPACKE_stprfs( stprfs_obj->matrix_layout, stprfs_obj->uplo,
                                  stprfs_obj->trans, stprfs_obj->diag, 
                                  stprfs_obj->n, stprfs_obj->nrhs,
                                  (const float *)stprfs_obj->a, 
                                  stprfs_obj->b, stprfs_obj->ldb,
								  stprfs_obj->x, stprfs_obj->ldx,
                                  stprfs_obj->ferr,
                                  stprfs_obj->berr
								  );
								  
    if( stprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stprfs is wrong\n", stprfs_obj->info );
    }
    if( stprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stprfs is wrong\n", 
        stprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stprfs_obj->diff_berr =  computeDiff_s( stprfs_obj->nrhs, 
                stprfs_obj->berr, stprfs_obj->berrref );
                
    stprfs_obj->diff_ferr =  computeDiff_s( stprfs_obj->nrhs, 
                stprfs_obj->ferr, stprfs_obj->ferrref );

}

TEST_F(stprfs_test, stprfs1) {
    EXPECT_NEAR(0.0, stprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, stprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(stprfs_test, stprfs2) {
    EXPECT_NEAR(0.0, stprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, stprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(stprfs_test, stprfs3) {
    EXPECT_NEAR(0.0, stprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, stprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(stprfs_test, stprfs4) {
    EXPECT_NEAR(0.0, stprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, stprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin tprfs_double_parameters  class definition */
class tprfs_double_parameters{
   public:
      int b_bufsize;
      void *hModule, *dModule;
      double diff_berr, diff_ferr, diff_xerr;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx;  //  leading dimension of 'x'
      double *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.
      double *x, *xref; //right-hand sides for the systems of equations.

      /* Output parameters */
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tprfs_double_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~tprfs_double_parameters (); 
};  /* end of tprfs_double_parameters  class definition */


/* Condtructor tprfs_double_parameters definition */
tprfs_double_parameters:: tprfs_double_parameters ( int matrix_layout_i, 
                               char uplo_i, char trans_i, char diag_i, 
                                    lapack_int n_i, lapack_int nrhs_i)
{
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = n;
    nrhs = nrhs_i;

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
   printf(" \n tprfs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
#endif

    b_bufsize = n*nrhs;
	
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*(n+1)/2) );
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       tprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, (n*(n+1)/2) );
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(double)));
    memcpy(xref, b, ( b_bufsize*sizeof(double)));

    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Condtructor  */

tprfs_double_parameters:: ~tprfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tprfs_double_parameters object: dedtructor invoked. \n");
#endif
   tprfs_free();
}

//  Test fixture class definition
class dtprfs_test  : public  ::testing::Test {
public:
   tprfs_double_parameters  *dtprfs_obj;
   void SetUp();  
   void TearDown () { delete dtprfs_obj; }
};


void dtprfs_test::SetUp(){

    /* LAPACKE DTPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dtptrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const double *a, 
                                     double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_dtptrs DTPTRS;

    typedef int (*Fptr_NL_LAPACKE_dtprfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int nrhs, const double* a,
			 const double* b, lapack_int ldb, const double* x,
			lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_dtprfs DTPRFS;

    dtprfs_obj = new tprfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    dtprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DTPTRS = (Fptr_NL_LAPACKE_dtptrs)dlsym(dtprfs_obj->hModule, "LAPACKE_dtptrs");
    ASSERT_TRUE(DTPTRS != NULL) << "failed to get the Netlib LAPACKE_dtptrs symbol";

    DTPRFS = (Fptr_NL_LAPACKE_dtprfs)dlsym(dtprfs_obj->hModule, "LAPACKE_dtprfs");
    ASSERT_TRUE(DTPRFS != NULL) << "failed to get the Netlib LAPACKE_dtprfs symbol";

    /* Generate i/ps to tprfs API from 'tptrs' API */
    dtprfs_obj->inforef = DTPTRS( dtprfs_obj->matrix_layout, dtprfs_obj->uplo,
                                  dtprfs_obj->trans, dtprfs_obj->diag, 
                                  dtprfs_obj->n, dtprfs_obj->nrhs,
                                  (const double *)dtprfs_obj->aref, 
                                  dtprfs_obj->xref,
                                  dtprfs_obj->ldb );

    dtprfs_obj->info = LAPACKE_dtptrs( dtprfs_obj->matrix_layout, dtprfs_obj->uplo,
                                  dtprfs_obj->trans, dtprfs_obj->diag, 
                                  dtprfs_obj->n, dtprfs_obj->nrhs,
                                  (const double *)dtprfs_obj->a, 
                                  dtprfs_obj->x,
                                  dtprfs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    dtprfs_obj->inforef = DTPRFS( dtprfs_obj->matrix_layout, dtprfs_obj->uplo,
                                  dtprfs_obj->trans, dtprfs_obj->diag, 
                                  dtprfs_obj->n, dtprfs_obj->nrhs,
                                  (const double *)dtprfs_obj->aref, 
                                  dtprfs_obj->bref, dtprfs_obj->ldb,
								  dtprfs_obj->xref, dtprfs_obj->ldx,
                                  dtprfs_obj->ferrref,
                                  dtprfs_obj->berrref
								  );

    dtprfs_obj->info = LAPACKE_dtprfs( dtprfs_obj->matrix_layout, dtprfs_obj->uplo,
                                  dtprfs_obj->trans, dtprfs_obj->diag, 
                                  dtprfs_obj->n, dtprfs_obj->nrhs,
                                  (const double *)dtprfs_obj->a, 
                                  dtprfs_obj->b, dtprfs_obj->ldb,
								  dtprfs_obj->x, dtprfs_obj->ldx,
                                  dtprfs_obj->ferr,
                                  dtprfs_obj->berr
								  );
								  
    if( dtprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtprfs is wrong\n", dtprfs_obj->info );
    }
    if( dtprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtprfs is wrong\n", 
        dtprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtprfs_obj->diff_berr =  computeDiff_d( dtprfs_obj->nrhs, 
                dtprfs_obj->berr, dtprfs_obj->berrref );
                
    dtprfs_obj->diff_ferr =  computeDiff_d( dtprfs_obj->nrhs, 
                dtprfs_obj->ferr, dtprfs_obj->ferrref );

}

TEST_F(dtprfs_test, dtprfs1) {
    EXPECT_NEAR(0.0, dtprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dtprfs_test, dtprfs2) {
    EXPECT_NEAR(0.0, dtprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dtprfs_test, dtprfs3) {
    EXPECT_NEAR(0.0, dtprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dtprfs_test, dtprfs4) {
    EXPECT_NEAR(0.0, dtprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin tprfs_scomplex_parameters  class definition */
class tprfs_scomplex_parameters{
   public:
      int b_bufsize;
      void *hModule, *dModule;
      float diff_berr, diff_ferr, diff_xerr;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx;  //  leading dimension of 'x'
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_float *x, *xref; //right-hand sides for the systems of equations.

      /* Output parameters */
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tprfs_scomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~tprfs_scomplex_parameters (); 
};  /* end of tprfs_scomplex_parameters  class definition */


/* Conctructor tprfs_scomplex_parameters definition */
tprfs_scomplex_parameters:: tprfs_scomplex_parameters ( int matrix_layout_i, 
                               char uplo_i, char trans_i, char diag_i, 
                                    lapack_int n_i, lapack_int nrhs_i)
{
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = n;
    nrhs = nrhs_i;

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
   printf(" \n tprfs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
#endif

    b_bufsize = n*nrhs;
	
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*(n+1)/2) );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       tprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (n*(n+1)/2) );
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_float)));

    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Conctructor  */

tprfs_scomplex_parameters:: ~tprfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tprfs_scomplex_parameters object: dectructor invoked. \n");
#endif
   tprfs_free();
}

//  Test fixture class definition
class ctprfs_test  : public  ::testing::Test {
public:
   tprfs_scomplex_parameters  *ctprfs_obj;
   void SetUp();  
   void TearDown () { delete ctprfs_obj; }
};


void ctprfs_test::SetUp(){

    /* LAPACKE CTPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ctptrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const lapack_complex_float *a, 
                                     lapack_complex_float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ctptrs CTPTRS;

    typedef int (*Fptr_NL_LAPACKE_ctprfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int nrhs, const lapack_complex_float* a,
			 const lapack_complex_float* b, lapack_int ldb, const lapack_complex_float* x,
			lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_ctprfs CTPRFS;

    ctprfs_obj = new tprfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    ctprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CTPTRS = (Fptr_NL_LAPACKE_ctptrs)dlsym(ctprfs_obj->hModule, "LAPACKE_ctptrs");
    ASSERT_TRUE(CTPTRS != NULL) << "failed to get the Netlib LAPACKE_ctptrs symbol";

    CTPRFS = (Fptr_NL_LAPACKE_ctprfs)dlsym(ctprfs_obj->hModule, "LAPACKE_ctprfs");
    ASSERT_TRUE(CTPRFS != NULL) << "failed to get the Netlib LAPACKE_ctprfs symbol";

    /* Generate i/ps to tprfs API from 'tptrs' API */
    ctprfs_obj->inforef = CTPTRS( ctprfs_obj->matrix_layout, ctprfs_obj->uplo,
                                  ctprfs_obj->trans, ctprfs_obj->diag, 
                                  ctprfs_obj->n, ctprfs_obj->nrhs,
                                  (const lapack_complex_float *)ctprfs_obj->aref, 
                                   ctprfs_obj->xref,
                                  ctprfs_obj->ldb );

    ctprfs_obj->info = LAPACKE_ctptrs( ctprfs_obj->matrix_layout, ctprfs_obj->uplo,
                                  ctprfs_obj->trans, ctprfs_obj->diag, 
                                  ctprfs_obj->n, ctprfs_obj->nrhs,
                                  (const lapack_complex_float *)ctprfs_obj->a, 
                                   ctprfs_obj->x,
                                  ctprfs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    ctprfs_obj->inforef = CTPRFS( ctprfs_obj->matrix_layout, ctprfs_obj->uplo,
                                  ctprfs_obj->trans, ctprfs_obj->diag, 
                                  ctprfs_obj->n, ctprfs_obj->nrhs,
                                  (const lapack_complex_float *)ctprfs_obj->aref, 
                                  
                                  ctprfs_obj->bref, ctprfs_obj->ldb,
								  ctprfs_obj->xref, ctprfs_obj->ldx,
                                  ctprfs_obj->ferrref,
                                  ctprfs_obj->berrref
								  );

    ctprfs_obj->info = LAPACKE_ctprfs( ctprfs_obj->matrix_layout, ctprfs_obj->uplo,
                                  ctprfs_obj->trans, ctprfs_obj->diag, 
                                  ctprfs_obj->n, ctprfs_obj->nrhs,
                                  (const lapack_complex_float *)ctprfs_obj->a, 
                                  
                                  ctprfs_obj->b, ctprfs_obj->ldb,
								  ctprfs_obj->x, ctprfs_obj->ldx,
                                  ctprfs_obj->ferr,
                                  ctprfs_obj->berr
								  );
								  
    if( ctprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctprfs is wrong\n", ctprfs_obj->info );
    }
    if( ctprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctprfs is wrong\n", 
        ctprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctprfs_obj->diff_berr =  computeDiff_s( ctprfs_obj->nrhs, 
                ctprfs_obj->berr, ctprfs_obj->berrref );
                
    ctprfs_obj->diff_ferr =  computeDiff_s( ctprfs_obj->nrhs, 
                ctprfs_obj->ferr, ctprfs_obj->ferrref );

}

TEST_F(ctprfs_test, ctprfs1) {
    EXPECT_NEAR(0.0, ctprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ctprfs_test, ctprfs2) {
    EXPECT_NEAR(0.0, ctprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ctprfs_test, ctprfs3) {
    EXPECT_NEAR(0.0, ctprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ctprfs_test, ctprfs4) {
    EXPECT_NEAR(0.0, ctprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin tprfs_dcomplex_parameters  class definition */
class tprfs_dcomplex_parameters{
   public:
      int b_bufsize;
      void *hModule, *dModule;
      double diff_berr, diff_ferr, diff_xerr;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx;  //  leading dimension of 'x'
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_double *x, *xref; //right-hand sides for the systems of equations.

      /* Output parameters */
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tprfs_dcomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~tprfs_dcomplex_parameters (); 
};  /* end of tprfs_dcomplex_parameters  class definition */


/* Conztructor tprfs_dcomplex_parameters definition */
tprfs_dcomplex_parameters:: tprfs_dcomplex_parameters ( int matrix_layout_i, 
                               char uplo_i, char trans_i, char diag_i, 
                                    lapack_int n_i, lapack_int nrhs_i)
{
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = n;
    nrhs = nrhs_i;

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
   printf(" \n tprfs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
#endif

    b_bufsize = n*nrhs;
	
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*(n+1)/2) );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       tprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref,(n*(n+1)/2) );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_double)));

    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Conztructor  */

tprfs_dcomplex_parameters:: ~tprfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tprfs_dcomplex_parameters object: deztructor invoked. \n");
#endif
   tprfs_free();
}

//  Test fixture class definition
class ztprfs_test  : public  ::testing::Test {
public:
   tprfs_dcomplex_parameters  *ztprfs_obj;
   void SetUp();  
   void TearDown () { delete ztprfs_obj; }
};


void ztprfs_test::SetUp(){

    /* LAPACKE CTPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ztptrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const lapack_complex_double *a, 
                                    lapack_complex_double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ztptrs CTPTRS;

    typedef int (*Fptr_NL_LAPACKE_ztprfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int nrhs, const lapack_complex_double* a,
			 const lapack_complex_double* b, lapack_int ldb, const lapack_complex_double* x,
			lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_ztprfs CTPRFS;

    ztprfs_obj = new tprfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    ztprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CTPTRS = (Fptr_NL_LAPACKE_ztptrs)dlsym(ztprfs_obj->hModule, "LAPACKE_ztptrs");
    ASSERT_TRUE(CTPTRS != NULL) << "failed to get the Netlib LAPACKE_ztptrs symbol";

    CTPRFS = (Fptr_NL_LAPACKE_ztprfs)dlsym(ztprfs_obj->hModule, "LAPACKE_ztprfs");
    ASSERT_TRUE(CTPRFS != NULL) << "failed to get the Netlib LAPACKE_ztprfs symbol";

    /* Generate i/ps to tprfs API from 'tptrs' API */
    ztprfs_obj->inforef = CTPTRS( ztprfs_obj->matrix_layout, ztprfs_obj->uplo,
                                  ztprfs_obj->trans, ztprfs_obj->diag, 
                                  ztprfs_obj->n, ztprfs_obj->nrhs,
                                  (const lapack_complex_double *)ztprfs_obj->aref, 
                                  ztprfs_obj->xref,
                                  ztprfs_obj->ldb );

    ztprfs_obj->info = LAPACKE_ztptrs( ztprfs_obj->matrix_layout, ztprfs_obj->uplo,
                                  ztprfs_obj->trans, ztprfs_obj->diag, 
                                  ztprfs_obj->n, ztprfs_obj->nrhs,
                                  (const lapack_complex_double *)ztprfs_obj->a, 
                                  ztprfs_obj->x,
                                  ztprfs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    ztprfs_obj->inforef = CTPRFS( ztprfs_obj->matrix_layout, ztprfs_obj->uplo,
                                  ztprfs_obj->trans, ztprfs_obj->diag, 
                                  ztprfs_obj->n, ztprfs_obj->nrhs,
                                  (const lapack_complex_double *)ztprfs_obj->aref, 
                                  ztprfs_obj->bref, ztprfs_obj->ldb,
								  ztprfs_obj->xref, ztprfs_obj->ldx,
                                  ztprfs_obj->ferrref,
                                  ztprfs_obj->berrref
								  );

    ztprfs_obj->info = LAPACKE_ztprfs( ztprfs_obj->matrix_layout, ztprfs_obj->uplo,
                                  ztprfs_obj->trans, ztprfs_obj->diag, 
                                  ztprfs_obj->n, ztprfs_obj->nrhs,
                                  (const lapack_complex_double *)ztprfs_obj->a, 
                                  ztprfs_obj->b, ztprfs_obj->ldb,
								  ztprfs_obj->x, ztprfs_obj->ldx,
                                  ztprfs_obj->ferr,
                                  ztprfs_obj->berr
								  );
								  
    if( ztprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztprfs is wrong\n", ztprfs_obj->info );
    }
    if( ztprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztprfs is wrong\n", 
        ztprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztprfs_obj->diff_berr =  computeDiff_d( ztprfs_obj->nrhs, 
                ztprfs_obj->berr, ztprfs_obj->berrref );
                
    ztprfs_obj->diff_ferr =  computeDiff_d( ztprfs_obj->nrhs, 
                ztprfs_obj->ferr, ztprfs_obj->ferrref );

}

TEST_F(ztprfs_test, ztprfs1) {
    EXPECT_NEAR(0.0, ztprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ztprfs_test, ztprfs2) {
    EXPECT_NEAR(0.0, ztprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ztprfs_test, ztprfs3) {
    EXPECT_NEAR(0.0, ztprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ztprfs_test, ztprfs4) {
    EXPECT_NEAR(0.0, ztprfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztprfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

