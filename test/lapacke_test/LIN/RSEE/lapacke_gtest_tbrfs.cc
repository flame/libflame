#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"
#define LAPACKE_TEST_VERBOSE (1)

#define tbrfs_free() \
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

/* Begin tbrfs_float_parameters  class definition */
class tbrfs_float_parameters{
   public:
      int b_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      void *hModule, *dModule;
      float diff_berr, diff_ferr, diff_xerr;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd; // The number of superdiagonals or subdiagonals in A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
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
      tbrfs_float_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i,
								lapack_int kd_i);
              
      ~tbrfs_float_parameters (); 
};  /* end of tbrfs_float_parameters  class definition */


/* Constructor tbrfs_float_parameters definition */
tbrfs_float_parameters:: tbrfs_float_parameters ( int matrix_layout_i, char uplo_i,
					char trans_i, char diag_i, lapack_int n_i, lapack_int nrhs_i,
								lapack_int kd_i)
{
    matrix_layout = matrix_layout_i;
    n = n_i;
	kd = kd_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = kd+1;
    nrhs = nrhs_i;
    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
		lda = n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n stbrfs : matrix_layout_i: %d  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d kd: %d \n",  matrix_layout, n, uplo, trans, diag, lda, ldb, nrhs, kd);
#endif

    b_bufsize = n*nrhs;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*lda));
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
       tbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(float)));
    memcpy(xref, b, ( b_bufsize*sizeof(float)));

    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    
   } /* end of Constructor  */

tbrfs_float_parameters:: ~tbrfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbrfs_float_parameters object: destructor invoked. \n");
#endif
   tbrfs_free();
}

//  Test fixture class definition
class stbrfs_test  : public  ::testing::Test {
public:
   tbrfs_float_parameters  *stbrfs_obj;
   void SetUp();  
   void TearDown () { delete stbrfs_obj; }
};


void stbrfs_test::SetUp(){

    /* LAPACKE STBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_stbtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int kd, lapack_int nrhs, const float *a, 
                                    lapack_int lda, float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_stbtrs STBTRS;

    typedef int (*Fptr_NL_LAPACKE_stbrfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int kd,
			lapack_int nrhs, const float* ab, lapack_int ldab,
			const float* b, lapack_int ldb, const float* x,
			lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_stbrfs STBRFS;

    stbrfs_obj = new tbrfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd );


    stbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    STBTRS = (Fptr_NL_LAPACKE_stbtrs)dlsym(stbrfs_obj->hModule, "LAPACKE_stbtrs");
    ASSERT_TRUE(STBTRS != NULL) << "failed to get the Netlib LAPACKE_stbtrs symbol";

    STBRFS = (Fptr_NL_LAPACKE_stbrfs)dlsym(stbrfs_obj->hModule, "LAPACKE_stbrfs");
    ASSERT_TRUE(STBRFS != NULL) << "failed to get the Netlib LAPACKE_stbrfs symbol";

    /* Generate i/ps to tbrfs API from 'tbrfs' API */
    stbrfs_obj->inforef = STBTRS( stbrfs_obj->matrix_layout, stbrfs_obj->uplo,
                                  stbrfs_obj->trans, stbrfs_obj->diag, 
                                  stbrfs_obj->n, stbrfs_obj->kd, stbrfs_obj->nrhs,
                                  (const float *)stbrfs_obj->aref, 
                                  stbrfs_obj->lda, stbrfs_obj->xref,
                                  stbrfs_obj->ldb );

    stbrfs_obj->info = LAPACKE_stbtrs( stbrfs_obj->matrix_layout, stbrfs_obj->uplo,
                                  stbrfs_obj->trans, stbrfs_obj->diag, 
                                  stbrfs_obj->n, stbrfs_obj->kd, stbrfs_obj->nrhs,
                                  (const float *)stbrfs_obj->a, 
                                  stbrfs_obj->lda, stbrfs_obj->x,
                                  stbrfs_obj->ldb );

    stbrfs_obj->inforef = STBRFS( stbrfs_obj->matrix_layout, stbrfs_obj->uplo,
                                  stbrfs_obj->trans, stbrfs_obj->diag, 
                                  stbrfs_obj->n, stbrfs_obj->kd, stbrfs_obj->nrhs,
                                  (const float *)stbrfs_obj->aref, 
                                  stbrfs_obj->lda, 
								  stbrfs_obj->bref, stbrfs_obj->ldb,
								  stbrfs_obj->xref, stbrfs_obj->ldb,
								  stbrfs_obj->ferrref, stbrfs_obj->berrref
								  );

    stbrfs_obj->info = LAPACKE_stbrfs( stbrfs_obj->matrix_layout, stbrfs_obj->uplo,
                                  stbrfs_obj->trans, stbrfs_obj->diag, 
                                  stbrfs_obj->n, stbrfs_obj->kd, stbrfs_obj->nrhs,
                                  (const float *)stbrfs_obj->a, 
                                  stbrfs_obj->lda, 
								  stbrfs_obj->b, stbrfs_obj->ldb,
								  stbrfs_obj->x, stbrfs_obj->ldb,
								  stbrfs_obj->ferr, stbrfs_obj->berr
								  );

    if( stbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stbrfs is wrong\n", stbrfs_obj->info );
    }
    if( stbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stbrfs is wrong\n", 
        stbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stbrfs_obj->diff_berr =  computeDiff_s( stbrfs_obj->nrhs, 
                stbrfs_obj->berr, stbrfs_obj->berrref );
                
    stbrfs_obj->diff_ferr =  computeDiff_s( stbrfs_obj->nrhs, 
                stbrfs_obj->ferr, stbrfs_obj->ferrref );
}

TEST_F(stbrfs_test, stbrfs1) {
    EXPECT_NEAR(0.0, stbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, stbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(stbrfs_test, stbrfs2) {
    EXPECT_NEAR(0.0, stbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, stbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(stbrfs_test, stbrfs3) {
    EXPECT_NEAR(0.0, stbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, stbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(stbrfs_test, stbrfs4) {
    EXPECT_NEAR(0.0, stbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, stbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin tbrfs_double_parameters  class definition */
class tbrfs_double_parameters{
   public:
      int b_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      void *hModule, *dModule;
      double diff_berr, diff_ferr, diff_xerr;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd; // The number of superdiagonals or subdiagonals in A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
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
      tbrfs_double_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i,
								lapack_int kd_i);
              
      ~tbrfs_double_parameters (); 
};  /* end of tbrfs_double_parameters  class definition */


/* Constructor tbrfs_double_parameters definition */
tbrfs_double_parameters:: tbrfs_double_parameters ( int matrix_layout_i, char uplo_i,
					char trans_i, char diag_i, lapack_int n_i, lapack_int nrhs_i,
								lapack_int kd_i)
{
    matrix_layout = matrix_layout_i;
    n = n_i;
	kd = kd_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = kd+1;
    nrhs = nrhs_i;
    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
		lda = n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n dtbrfs : matrix_layout_i: %d  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d kd: %d \n",  matrix_layout, n, uplo, trans, diag, lda, ldb, nrhs, kd);
#endif

    b_bufsize = n*nrhs;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*lda));
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
       tbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(double)));
    memcpy(xref, b, ( b_bufsize*sizeof(double)));

    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    
   } /* end of Constructor  */

tbrfs_double_parameters:: ~tbrfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbrfs_double_parameters object: destructor invoked. \n");
#endif
   tbrfs_free();
}

//  Test fixture class definition
class dtbrfs_test  : public  ::testing::Test {
public:
   tbrfs_double_parameters  *dtbrfs_obj;
   void SetUp();  
   void TearDown () { delete dtbrfs_obj; }
};


void dtbrfs_test::SetUp(){

    /* LAPACKE DTBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dtbtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int kd, lapack_int nrhs, const double *a, 
                                    lapack_int lda, double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_dtbtrs DTBTRS;

    typedef int (*Fptr_NL_LAPACKE_dtbrfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int kd,
			lapack_int nrhs, const double* ab, lapack_int ldab,
			const double* b, lapack_int ldb, const double* x,
			lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_dtbrfs DTBRFS;

    dtbrfs_obj = new tbrfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd );


    dtbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DTBTRS = (Fptr_NL_LAPACKE_dtbtrs)dlsym(dtbrfs_obj->hModule, "LAPACKE_dtbtrs");
    ASSERT_TRUE(DTBTRS != NULL) << "failed to get the Netlib LAPACKE_dtbtrs symbol";

    DTBRFS = (Fptr_NL_LAPACKE_dtbrfs)dlsym(dtbrfs_obj->hModule, "LAPACKE_dtbrfs");
    ASSERT_TRUE(DTBRFS != NULL) << "failed to get the Netlib LAPACKE_dtbrfs symbol";

    /* Generate i/ps to tbrfs API from 'tbrfs' API */
    dtbrfs_obj->inforef = DTBTRS( dtbrfs_obj->matrix_layout, dtbrfs_obj->uplo,
                                  dtbrfs_obj->trans, dtbrfs_obj->diag, 
                                  dtbrfs_obj->n, dtbrfs_obj->kd, dtbrfs_obj->nrhs,
                                  (const double *)dtbrfs_obj->aref, 
                                  dtbrfs_obj->lda, dtbrfs_obj->xref,
                                  dtbrfs_obj->ldb );

    dtbrfs_obj->info = LAPACKE_dtbtrs( dtbrfs_obj->matrix_layout, dtbrfs_obj->uplo,
                                  dtbrfs_obj->trans, dtbrfs_obj->diag, 
                                  dtbrfs_obj->n, dtbrfs_obj->kd, dtbrfs_obj->nrhs,
                                  (const double *)dtbrfs_obj->a, 
                                  dtbrfs_obj->lda, dtbrfs_obj->x,
                                  dtbrfs_obj->ldb );

    dtbrfs_obj->inforef = DTBRFS( dtbrfs_obj->matrix_layout, dtbrfs_obj->uplo,
                                  dtbrfs_obj->trans, dtbrfs_obj->diag, 
                                  dtbrfs_obj->n, dtbrfs_obj->kd, dtbrfs_obj->nrhs,
                                  (const double *)dtbrfs_obj->aref, 
                                  dtbrfs_obj->lda, 
								  dtbrfs_obj->bref, dtbrfs_obj->ldb,
								  dtbrfs_obj->xref, dtbrfs_obj->ldb,
								  dtbrfs_obj->ferrref, dtbrfs_obj->berrref
								  );

    dtbrfs_obj->info = LAPACKE_dtbrfs( dtbrfs_obj->matrix_layout, dtbrfs_obj->uplo,
                                  dtbrfs_obj->trans, dtbrfs_obj->diag, 
                                  dtbrfs_obj->n, dtbrfs_obj->kd, dtbrfs_obj->nrhs,
                                  (const double *)dtbrfs_obj->a, 
                                  dtbrfs_obj->lda, 
								  dtbrfs_obj->b, dtbrfs_obj->ldb,
								  dtbrfs_obj->x, dtbrfs_obj->ldb,
								  dtbrfs_obj->ferr, dtbrfs_obj->berr
								  );

    if( dtbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtbrfs is wrong\n", dtbrfs_obj->info );
    }
    if( dtbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtbrfs is wrong\n", 
        dtbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtbrfs_obj->diff_berr =  computeDiff_d( dtbrfs_obj->nrhs, 
                dtbrfs_obj->berr, dtbrfs_obj->berrref );
                
    dtbrfs_obj->diff_ferr =  computeDiff_d( dtbrfs_obj->nrhs, 
                dtbrfs_obj->ferr, dtbrfs_obj->ferrref );
}

TEST_F(dtbrfs_test, dtbrfs1) {
    EXPECT_NEAR(0.0, dtbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dtbrfs_test, dtbrfs2) {
    EXPECT_NEAR(0.0, dtbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dtbrfs_test, dtbrfs3) {
    EXPECT_NEAR(0.0, dtbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dtbrfs_test, dtbrfs4) {
    EXPECT_NEAR(0.0, dtbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin tbrfs_scomplex_parameters  class definition */
class tbrfs_scomplex_parameters{
   public:
      int b_bufsize;
      lapack_complex_float diff; // capture difference between ref o/p & libflame lapacke o/p.
      void *hModule, *dModule;
      float diff_berr, diff_ferr, diff_xerr;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd; // The number of superdiagonals or subdiagonals in A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
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
      tbrfs_scomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i,
								lapack_int kd_i);
              
      ~tbrfs_scomplex_parameters (); 
};  /* end of tbrfs_scomplex_parameters  class definition */


/* Constructor tbrfs_scomplex_parameters definition */
tbrfs_scomplex_parameters:: tbrfs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
					char trans_i, char diag_i, lapack_int n_i, lapack_int nrhs_i,
								lapack_int kd_i)
{
    matrix_layout = matrix_layout_i;
    n = n_i;
	kd = kd_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = kd+1;
    nrhs = nrhs_i;
    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
		lda = n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n ctbrfs : matrix_layout_i: %d  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d kd: %d \n",  matrix_layout, n, uplo, trans, diag, lda, ldb, nrhs, kd);
#endif

    b_bufsize = n*nrhs;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*lda));
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
       tbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_float)));

    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    
   } /* end of Constructor  */

tbrfs_scomplex_parameters:: ~tbrfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbrfs_scomplex_parameters object: destructor invoked. \n");
#endif
   tbrfs_free();
}

//  Test fixture class definition
class ctbrfs_test  : public  ::testing::Test {
public:
   tbrfs_scomplex_parameters  *ctbrfs_obj;
   void SetUp();  
   void TearDown () { delete ctbrfs_obj; }
};


void ctbrfs_test::SetUp(){

    /* LAPACKE CTBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ctbtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int kd, lapack_int nrhs, const lapack_complex_float *a, 
                                    lapack_int lda, lapack_complex_float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ctbtrs CTBTRS;

    typedef int (*Fptr_NL_LAPACKE_ctbrfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int kd,
			lapack_int nrhs, const lapack_complex_float* ab, lapack_int ldab,
			const lapack_complex_float* b, lapack_int ldb, const lapack_complex_float* x,
			lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_ctbrfs CTBRFS;

    ctbrfs_obj = new tbrfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd );


    ctbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CTBTRS = (Fptr_NL_LAPACKE_ctbtrs)dlsym(ctbrfs_obj->hModule, "LAPACKE_ctbtrs");
    ASSERT_TRUE(CTBTRS != NULL) << "failed to get the Netlib LAPACKE_ctbtrs symbol";

    CTBRFS = (Fptr_NL_LAPACKE_ctbrfs)dlsym(ctbrfs_obj->hModule, "LAPACKE_ctbrfs");
    ASSERT_TRUE(CTBRFS != NULL) << "failed to get the Netlib LAPACKE_ctbrfs symbol";

    /* Generate i/ps to tbrfs API from 'tbrfs' API */
    ctbrfs_obj->inforef = CTBTRS( ctbrfs_obj->matrix_layout, ctbrfs_obj->uplo,
                                  ctbrfs_obj->trans, ctbrfs_obj->diag, 
                                  ctbrfs_obj->n, ctbrfs_obj->kd, ctbrfs_obj->nrhs,
                                  (const lapack_complex_float *)ctbrfs_obj->aref, 
                                  ctbrfs_obj->lda, ctbrfs_obj->xref,
                                  ctbrfs_obj->ldb );

    ctbrfs_obj->info = LAPACKE_ctbtrs( ctbrfs_obj->matrix_layout, ctbrfs_obj->uplo,
                                  ctbrfs_obj->trans, ctbrfs_obj->diag, 
                                  ctbrfs_obj->n, ctbrfs_obj->kd, ctbrfs_obj->nrhs,
                                  (const lapack_complex_float *)ctbrfs_obj->a, 
                                  ctbrfs_obj->lda, ctbrfs_obj->x,
                                  ctbrfs_obj->ldb );
#if LAPACKE_TEST_VERBOSE
   printf(" \n ldx: %d  n: %d  lda: %d \n", ctbrfs_obj->ldb,  ctbrfs_obj->n, ctbrfs_obj->lda);
#endif

    ctbrfs_obj->inforef = CTBRFS( ctbrfs_obj->matrix_layout, ctbrfs_obj->uplo,
                                  ctbrfs_obj->trans, ctbrfs_obj->diag, 
                                  ctbrfs_obj->n, ctbrfs_obj->kd, ctbrfs_obj->nrhs,
                                  (const lapack_complex_float *)ctbrfs_obj->aref, 
                                  ctbrfs_obj->lda, 
								  ctbrfs_obj->bref, ctbrfs_obj->ldb,
								  ctbrfs_obj->xref, ctbrfs_obj->ldb,
								  ctbrfs_obj->ferrref, ctbrfs_obj->berrref
								  );

    ctbrfs_obj->info = LAPACKE_ctbrfs( ctbrfs_obj->matrix_layout, ctbrfs_obj->uplo,
                                  ctbrfs_obj->trans, ctbrfs_obj->diag, 
                                  ctbrfs_obj->n, ctbrfs_obj->kd, ctbrfs_obj->nrhs,
                                  (const lapack_complex_float *)ctbrfs_obj->a, 
                                  ctbrfs_obj->lda, 
								  ctbrfs_obj->b, ctbrfs_obj->ldb,
								  ctbrfs_obj->x, ctbrfs_obj->ldb,
								  ctbrfs_obj->ferr, ctbrfs_obj->berr
								  );

    if( ctbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctbrfs is wrong\n", ctbrfs_obj->info );
    }
    if( ctbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctbrfs is wrong\n", 
        ctbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctbrfs_obj->diff_berr =  computeDiff_s( ctbrfs_obj->nrhs, 
                ctbrfs_obj->berr, ctbrfs_obj->berrref );
                
    ctbrfs_obj->diff_ferr =  computeDiff_s( ctbrfs_obj->nrhs, 
                ctbrfs_obj->ferr, ctbrfs_obj->ferrref );
}

TEST_F(ctbrfs_test, ctbrfs1) {
    EXPECT_NEAR(0.0, ctbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ctbrfs_test, ctbrfs2) {
    EXPECT_NEAR(0.0, ctbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ctbrfs_test, ctbrfs3) {
    EXPECT_NEAR(0.0, ctbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ctbrfs_test, ctbrfs4) {
    EXPECT_NEAR(0.0, ctbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin tbrfs_dcomplex_parameters  class definition */
class tbrfs_dcomplex_parameters{
   public:
      int b_bufsize;
      lapack_complex_double diff; // capture difference between ref o/p & libflame lapacke o/p.
      void *hModule, *dModule;
      double diff_berr, diff_ferr, diff_xerr;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd; // The number of superdiagonals or subdiagonals in A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
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
      tbrfs_dcomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i,
								lapack_int kd_i);
              
      ~tbrfs_dcomplex_parameters (); 
};  /* end of tbrfs_dcomplex_parameters  class definition */


/* Constructor tbrfs_dcomplex_parameters definition */
tbrfs_dcomplex_parameters:: tbrfs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
					char trans_i, char diag_i, lapack_int n_i, lapack_int nrhs_i,
								lapack_int kd_i)
{
    matrix_layout = matrix_layout_i;
    n = n_i;
	kd = kd_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = kd+1;
    nrhs = nrhs_i;
    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
		lda = n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n ztbrfs : matrix_layout_i: %d  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d kd: %d \n",  matrix_layout, n, uplo, trans, diag, lda, ldb, nrhs, kd);
#endif

    b_bufsize = n*nrhs;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*lda));
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
       tbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_double)));

    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    
   } /* end of Constructor  */

tbrfs_dcomplex_parameters:: ~tbrfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbrfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   tbrfs_free();
}

//  Test fixture class definition
class ztbrfs_test  : public  ::testing::Test {
public:
   tbrfs_dcomplex_parameters  *ztbrfs_obj;
   void SetUp();  
   void TearDown () { delete ztbrfs_obj; }
};


void ztbrfs_test::SetUp(){

    /* LAPACKE ZTBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ztbtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int kd, lapack_int nrhs, const lapack_complex_double *a, 
                                    lapack_int lda, lapack_complex_double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ztbtrs ZTBTRS;

    typedef int (*Fptr_NL_LAPACKE_ztbrfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int kd,
			lapack_int nrhs, const lapack_complex_double* ab, lapack_int ldab,
			const lapack_complex_double* b, lapack_int ldb, const lapack_complex_double* x,
			lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_ztbrfs ZTBRFS;

    ztbrfs_obj = new tbrfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd );


    ztbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZTBTRS = (Fptr_NL_LAPACKE_ztbtrs)dlsym(ztbrfs_obj->hModule, "LAPACKE_ztbtrs");
    ASSERT_TRUE(ZTBTRS != NULL) << "failed to get the Netlib LAPACKE_ztbtrs symbol";

    ZTBRFS = (Fptr_NL_LAPACKE_ztbrfs)dlsym(ztbrfs_obj->hModule, "LAPACKE_ztbrfs");
    ASSERT_TRUE(ZTBRFS != NULL) << "failed to get the Netlib LAPACKE_ztbrfs symbol";

    /* Generate i/ps to tbrfs API from 'tbrfs' API */
    ztbrfs_obj->inforef = ZTBTRS( ztbrfs_obj->matrix_layout, ztbrfs_obj->uplo,
                                  ztbrfs_obj->trans, ztbrfs_obj->diag, 
                                  ztbrfs_obj->n, ztbrfs_obj->kd, ztbrfs_obj->nrhs,
                                  (const lapack_complex_double *)ztbrfs_obj->aref, 
                                  ztbrfs_obj->lda, ztbrfs_obj->xref,
                                  ztbrfs_obj->ldb );

    ztbrfs_obj->info = LAPACKE_ztbtrs( ztbrfs_obj->matrix_layout, ztbrfs_obj->uplo,
                                  ztbrfs_obj->trans, ztbrfs_obj->diag, 
                                  ztbrfs_obj->n, ztbrfs_obj->kd, ztbrfs_obj->nrhs,
                                  (const lapack_complex_double *)ztbrfs_obj->a, 
                                  ztbrfs_obj->lda, ztbrfs_obj->x,
                                  ztbrfs_obj->ldb );

    ztbrfs_obj->inforef = ZTBRFS( ztbrfs_obj->matrix_layout, ztbrfs_obj->uplo,
                                  ztbrfs_obj->trans, ztbrfs_obj->diag, 
                                  ztbrfs_obj->n, ztbrfs_obj->kd, ztbrfs_obj->nrhs,
                                  (const lapack_complex_double *)ztbrfs_obj->aref, 
                                  ztbrfs_obj->lda, 
								  ztbrfs_obj->bref, ztbrfs_obj->ldb,
								  ztbrfs_obj->xref, ztbrfs_obj->ldb,
								  ztbrfs_obj->ferrref, ztbrfs_obj->berrref
								  );

    ztbrfs_obj->info = LAPACKE_ztbrfs( ztbrfs_obj->matrix_layout, ztbrfs_obj->uplo,
                                  ztbrfs_obj->trans, ztbrfs_obj->diag, 
                                  ztbrfs_obj->n, ztbrfs_obj->kd, ztbrfs_obj->nrhs,
                                  (const lapack_complex_double *)ztbrfs_obj->a, 
                                  ztbrfs_obj->lda, 
								  ztbrfs_obj->b, ztbrfs_obj->ldb,
								  ztbrfs_obj->x, ztbrfs_obj->ldb,
								  ztbrfs_obj->ferr, ztbrfs_obj->berr
								  );

    if( ztbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztbrfs is wrong\n", ztbrfs_obj->info );
    }
    if( ztbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztbrfs is wrong\n", 
        ztbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztbrfs_obj->diff_berr =  computeDiff_d( ztbrfs_obj->nrhs, 
                ztbrfs_obj->berr, ztbrfs_obj->berrref );
                
    ztbrfs_obj->diff_ferr =  computeDiff_d( ztbrfs_obj->nrhs, 
                ztbrfs_obj->ferr, ztbrfs_obj->ferrref );
}

TEST_F(ztbrfs_test, ztbrfs1) {
    EXPECT_NEAR(0.0, ztbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ztbrfs_test, ztbrfs2) {
    EXPECT_NEAR(0.0, ztbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ztbrfs_test, ztbrfs3) {
    EXPECT_NEAR(0.0, ztbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ztbrfs_test, ztbrfs4) {
    EXPECT_NEAR(0.0, ztbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}
