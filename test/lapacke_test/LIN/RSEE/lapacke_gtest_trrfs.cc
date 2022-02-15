#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define trrfs_free() \
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

/* Begin trrfs_float_parameters  class definition */
class trrfs_float_parameters{
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
      trrfs_float_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~trrfs_float_parameters (); 
};  /* end of trrfs_float_parameters  class definition */


/* Constructor trrfs_float_parameters definition */
trrfs_float_parameters:: trrfs_float_parameters ( int matrix_layout_i, 
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
   printf(" \n trrfs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
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
       trrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(float)));
    memcpy(xref, b, ( b_bufsize*sizeof(float)));

    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Constructor  */

trrfs_float_parameters:: ~trrfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trrfs_float_parameters object: destructor invoked. \n");
#endif
   trrfs_free();
}

//  Test fixture class definition
class strrfs_test  : public  ::testing::Test {
public:
   trrfs_float_parameters  *strrfs_obj;
   void SetUp();  
   void TearDown () { delete strrfs_obj; }
};


void strrfs_test::SetUp(){

    /* LAPACKE STRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_strtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const float *a, 
                                    lapack_int lda, float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_strtrs STRTRS;

    typedef int (*Fptr_NL_LAPACKE_strrfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int nrhs, const float* a,
			lapack_int lda, const float* b, lapack_int ldb, const float* x,
			lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_strrfs STRRFS;

    strrfs_obj = new trrfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    strrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    strrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(strrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(strrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    STRTRS = (Fptr_NL_LAPACKE_strtrs)dlsym(strrfs_obj->hModule, "LAPACKE_strtrs");
    ASSERT_TRUE(STRTRS != NULL) << "failed to get the Netlib LAPACKE_strtrs symbol";

    STRRFS = (Fptr_NL_LAPACKE_strrfs)dlsym(strrfs_obj->hModule, "LAPACKE_strrfs");
    ASSERT_TRUE(STRRFS != NULL) << "failed to get the Netlib LAPACKE_strrfs symbol";

    /* Generate i/ps to trrfs API from 'trtrs' API */
    strrfs_obj->inforef = STRTRS( strrfs_obj->matrix_layout, strrfs_obj->uplo,
                                  strrfs_obj->trans, strrfs_obj->diag, 
                                  strrfs_obj->n, strrfs_obj->nrhs,
                                  (const float *)strrfs_obj->aref, 
                                  strrfs_obj->lda, strrfs_obj->xref,
                                  strrfs_obj->ldb );

    strrfs_obj->info = LAPACKE_strtrs( strrfs_obj->matrix_layout, strrfs_obj->uplo,
                                  strrfs_obj->trans, strrfs_obj->diag, 
                                  strrfs_obj->n, strrfs_obj->nrhs,
                                  (const float *)strrfs_obj->a, 
                                  strrfs_obj->lda, strrfs_obj->x,
                                  strrfs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    strrfs_obj->inforef = STRRFS( strrfs_obj->matrix_layout, strrfs_obj->uplo,
                                  strrfs_obj->trans, strrfs_obj->diag, 
                                  strrfs_obj->n, strrfs_obj->nrhs,
                                  (const float *)strrfs_obj->aref, 
                                  strrfs_obj->lda,
                                  strrfs_obj->bref, strrfs_obj->ldb,
								  strrfs_obj->xref, strrfs_obj->ldx,
                                  strrfs_obj->ferrref,
                                  strrfs_obj->berrref
								  );

    strrfs_obj->info = LAPACKE_strrfs( strrfs_obj->matrix_layout, strrfs_obj->uplo,
                                  strrfs_obj->trans, strrfs_obj->diag, 
                                  strrfs_obj->n, strrfs_obj->nrhs,
                                  (const float *)strrfs_obj->a, 
                                  strrfs_obj->lda,
                                  strrfs_obj->b, strrfs_obj->ldb,
								  strrfs_obj->x, strrfs_obj->ldx,
                                  strrfs_obj->ferr,
                                  strrfs_obj->berr
								  );
								  
    if( strrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_strrfs is wrong\n", strrfs_obj->info );
    }
    if( strrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_strrfs is wrong\n", 
        strrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    strrfs_obj->diff_berr =  computeDiff_s( strrfs_obj->nrhs, 
                strrfs_obj->berr, strrfs_obj->berrref );
                
    strrfs_obj->diff_ferr =  computeDiff_s( strrfs_obj->nrhs, 
                strrfs_obj->ferr, strrfs_obj->ferrref );

}

TEST_F(strrfs_test, strrfs1) {
    EXPECT_NEAR(0.0, strrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, strrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(strrfs_test, strrfs2) {
    EXPECT_NEAR(0.0, strrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, strrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(strrfs_test, strrfs3) {
    EXPECT_NEAR(0.0, strrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, strrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(strrfs_test, strrfs4) {
    EXPECT_NEAR(0.0, strrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, strrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin trrfs_double_parameters  class definition */
class trrfs_double_parameters{
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
      trrfs_double_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~trrfs_double_parameters (); 
};  /* end of trrfs_double_parameters  class definition */


/* Condtructor trrfs_double_parameters definition */
trrfs_double_parameters:: trrfs_double_parameters ( int matrix_layout_i, 
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
   printf(" \n trrfs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
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
       trrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(double)));
    memcpy(xref, b, ( b_bufsize*sizeof(double)));

    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Condtructor  */

trrfs_double_parameters:: ~trrfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trrfs_double_parameters object: dedtructor invoked. \n");
#endif
   trrfs_free();
}

//  Test fixture class definition
class dtrrfs_test  : public  ::testing::Test {
public:
   trrfs_double_parameters  *dtrrfs_obj;
   void SetUp();  
   void TearDown () { delete dtrrfs_obj; }
};


void dtrrfs_test::SetUp(){

    /* LAPACKE DTRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dtrtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const double *a, 
                                    lapack_int lda, double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_dtrtrs DTRTRS;

    typedef int (*Fptr_NL_LAPACKE_dtrrfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int nrhs, const double* a,
			lapack_int lda, const double* b, lapack_int ldb, const double* x,
			lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_dtrrfs DTRRFS;

    dtrrfs_obj = new trrfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    dtrrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtrrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtrrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtrrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DTRTRS = (Fptr_NL_LAPACKE_dtrtrs)dlsym(dtrrfs_obj->hModule, "LAPACKE_dtrtrs");
    ASSERT_TRUE(DTRTRS != NULL) << "failed to get the Netlib LAPACKE_dtrtrs symbol";

    DTRRFS = (Fptr_NL_LAPACKE_dtrrfs)dlsym(dtrrfs_obj->hModule, "LAPACKE_dtrrfs");
    ASSERT_TRUE(DTRRFS != NULL) << "failed to get the Netlib LAPACKE_dtrrfs symbol";

    /* Generate i/ps to trrfs API from 'trtrs' API */
    dtrrfs_obj->inforef = DTRTRS( dtrrfs_obj->matrix_layout, dtrrfs_obj->uplo,
                                  dtrrfs_obj->trans, dtrrfs_obj->diag, 
                                  dtrrfs_obj->n, dtrrfs_obj->nrhs,
                                  (const double *)dtrrfs_obj->aref, 
                                  dtrrfs_obj->lda, dtrrfs_obj->xref,
                                  dtrrfs_obj->ldb );

    dtrrfs_obj->info = LAPACKE_dtrtrs( dtrrfs_obj->matrix_layout, dtrrfs_obj->uplo,
                                  dtrrfs_obj->trans, dtrrfs_obj->diag, 
                                  dtrrfs_obj->n, dtrrfs_obj->nrhs,
                                  (const double *)dtrrfs_obj->a, 
                                  dtrrfs_obj->lda, dtrrfs_obj->x,
                                  dtrrfs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    dtrrfs_obj->inforef = DTRRFS( dtrrfs_obj->matrix_layout, dtrrfs_obj->uplo,
                                  dtrrfs_obj->trans, dtrrfs_obj->diag, 
                                  dtrrfs_obj->n, dtrrfs_obj->nrhs,
                                  (const double *)dtrrfs_obj->aref, 
                                  dtrrfs_obj->lda,
                                  dtrrfs_obj->bref, dtrrfs_obj->ldb,
								  dtrrfs_obj->xref, dtrrfs_obj->ldx,
                                  dtrrfs_obj->ferrref,
                                  dtrrfs_obj->berrref
								  );

    dtrrfs_obj->info = LAPACKE_dtrrfs( dtrrfs_obj->matrix_layout, dtrrfs_obj->uplo,
                                  dtrrfs_obj->trans, dtrrfs_obj->diag, 
                                  dtrrfs_obj->n, dtrrfs_obj->nrhs,
                                  (const double *)dtrrfs_obj->a, 
                                  dtrrfs_obj->lda,
                                  dtrrfs_obj->b, dtrrfs_obj->ldb,
								  dtrrfs_obj->x, dtrrfs_obj->ldx,
                                  dtrrfs_obj->ferr,
                                  dtrrfs_obj->berr
								  );
								  
    if( dtrrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtrrfs is wrong\n", dtrrfs_obj->info );
    }
    if( dtrrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtrrfs is wrong\n", 
        dtrrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtrrfs_obj->diff_berr =  computeDiff_d( dtrrfs_obj->nrhs, 
                dtrrfs_obj->berr, dtrrfs_obj->berrref );
                
    dtrrfs_obj->diff_ferr =  computeDiff_d( dtrrfs_obj->nrhs, 
                dtrrfs_obj->ferr, dtrrfs_obj->ferrref );

}

TEST_F(dtrrfs_test, dtrrfs1) {
    EXPECT_NEAR(0.0, dtrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dtrrfs_test, dtrrfs2) {
    EXPECT_NEAR(0.0, dtrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dtrrfs_test, dtrrfs3) {
    EXPECT_NEAR(0.0, dtrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dtrrfs_test, dtrrfs4) {
    EXPECT_NEAR(0.0, dtrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dtrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin trrfs_scomplex_parameters  class definition */
class trrfs_scomplex_parameters{
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
      trrfs_scomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~trrfs_scomplex_parameters (); 
};  /* end of trrfs_scomplex_parameters  class definition */


/* Conctructor trrfs_scomplex_parameters definition */
trrfs_scomplex_parameters:: trrfs_scomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n trrfs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
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
       trrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_float)));

    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Conctructor  */

trrfs_scomplex_parameters:: ~trrfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trrfs_scomplex_parameters object: dectructor invoked. \n");
#endif
   trrfs_free();
}

//  Test fixture class definition
class ctrrfs_test  : public  ::testing::Test {
public:
   trrfs_scomplex_parameters  *ctrrfs_obj;
   void SetUp();  
   void TearDown () { delete ctrrfs_obj; }
};


void ctrrfs_test::SetUp(){

    /* LAPACKE CTRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ctrtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const lapack_complex_float *a, 
                                    lapack_int lda, lapack_complex_float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ctrtrs CTRTRS;

    typedef int (*Fptr_NL_LAPACKE_ctrrfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int nrhs, const lapack_complex_float* a,
			lapack_int lda, const lapack_complex_float* b, lapack_int ldb, const lapack_complex_float* x,
			lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_ctrrfs CTRRFS;

    ctrrfs_obj = new trrfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    ctrrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctrrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctrrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctrrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CTRTRS = (Fptr_NL_LAPACKE_ctrtrs)dlsym(ctrrfs_obj->hModule, "LAPACKE_ctrtrs");
    ASSERT_TRUE(CTRTRS != NULL) << "failed to get the Netlib LAPACKE_ctrtrs symbol";

    CTRRFS = (Fptr_NL_LAPACKE_ctrrfs)dlsym(ctrrfs_obj->hModule, "LAPACKE_ctrrfs");
    ASSERT_TRUE(CTRRFS != NULL) << "failed to get the Netlib LAPACKE_ctrrfs symbol";

    /* Generate i/ps to trrfs API from 'trtrs' API */
    ctrrfs_obj->inforef = CTRTRS( ctrrfs_obj->matrix_layout, ctrrfs_obj->uplo,
                                  ctrrfs_obj->trans, ctrrfs_obj->diag, 
                                  ctrrfs_obj->n, ctrrfs_obj->nrhs,
                                  (const lapack_complex_float *)ctrrfs_obj->aref, 
                                  ctrrfs_obj->lda, ctrrfs_obj->xref,
                                  ctrrfs_obj->ldb );

    ctrrfs_obj->info = LAPACKE_ctrtrs( ctrrfs_obj->matrix_layout, ctrrfs_obj->uplo,
                                  ctrrfs_obj->trans, ctrrfs_obj->diag, 
                                  ctrrfs_obj->n, ctrrfs_obj->nrhs,
                                  (const lapack_complex_float *)ctrrfs_obj->a, 
                                  ctrrfs_obj->lda, ctrrfs_obj->x,
                                  ctrrfs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    ctrrfs_obj->inforef = CTRRFS( ctrrfs_obj->matrix_layout, ctrrfs_obj->uplo,
                                  ctrrfs_obj->trans, ctrrfs_obj->diag, 
                                  ctrrfs_obj->n, ctrrfs_obj->nrhs,
                                  (const lapack_complex_float *)ctrrfs_obj->aref, 
                                  ctrrfs_obj->lda,
                                  ctrrfs_obj->bref, ctrrfs_obj->ldb,
								  ctrrfs_obj->xref, ctrrfs_obj->ldx,
                                  ctrrfs_obj->ferrref,
                                  ctrrfs_obj->berrref
								  );

    ctrrfs_obj->info = LAPACKE_ctrrfs( ctrrfs_obj->matrix_layout, ctrrfs_obj->uplo,
                                  ctrrfs_obj->trans, ctrrfs_obj->diag, 
                                  ctrrfs_obj->n, ctrrfs_obj->nrhs,
                                  (const lapack_complex_float *)ctrrfs_obj->a, 
                                  ctrrfs_obj->lda,
                                  ctrrfs_obj->b, ctrrfs_obj->ldb,
								  ctrrfs_obj->x, ctrrfs_obj->ldx,
                                  ctrrfs_obj->ferr,
                                  ctrrfs_obj->berr
								  );
								  
    if( ctrrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctrrfs is wrong\n", ctrrfs_obj->info );
    }
    if( ctrrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctrrfs is wrong\n", 
        ctrrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctrrfs_obj->diff_berr =  computeDiff_s( ctrrfs_obj->nrhs, 
                ctrrfs_obj->berr, ctrrfs_obj->berrref );
                
    ctrrfs_obj->diff_ferr =  computeDiff_s( ctrrfs_obj->nrhs, 
                ctrrfs_obj->ferr, ctrrfs_obj->ferrref );

}

TEST_F(ctrrfs_test, ctrrfs1) {
    EXPECT_NEAR(0.0, ctrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ctrrfs_test, ctrrfs2) {
    EXPECT_NEAR(0.0, ctrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ctrrfs_test, ctrrfs3) {
    EXPECT_NEAR(0.0, ctrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ctrrfs_test, ctrrfs4) {
    EXPECT_NEAR(0.0, ctrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ctrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin trrfs_dcomplex_parameters  class definition */
class trrfs_dcomplex_parameters{
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
      trrfs_dcomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int nrhs_i);
              
      ~trrfs_dcomplex_parameters (); 
};  /* end of trrfs_dcomplex_parameters  class definition */


/* Conztructor trrfs_dcomplex_parameters definition */
trrfs_dcomplex_parameters:: trrfs_dcomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n trrfs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
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
       trrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, ( b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, ( b_bufsize*sizeof(lapack_complex_double)));

    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

   } /* end of Conztructor  */

trrfs_dcomplex_parameters:: ~trrfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trrfs_dcomplex_parameters object: deztructor invoked. \n");
#endif
   trrfs_free();
}

//  Test fixture class definition
class ztrrfs_test  : public  ::testing::Test {
public:
   trrfs_dcomplex_parameters  *ztrrfs_obj;
   void SetUp();  
   void TearDown () { delete ztrrfs_obj; }
};


void ztrrfs_test::SetUp(){

    /* LAPACKE CTRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ztrtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const lapack_complex_double *a, 
                                    lapack_int lda, lapack_complex_double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ztrtrs CTRTRS;

    typedef int (*Fptr_NL_LAPACKE_ztrrfs) (int matrix_layout, char uplo,
			char trans, char diag, lapack_int n, lapack_int nrhs, const lapack_complex_double* a,
			lapack_int lda, const lapack_complex_double* b, lapack_int ldb, const lapack_complex_double* x,
			lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_ztrrfs CTRRFS;

    ztrrfs_obj = new trrfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    ztrrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztrrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztrrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztrrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CTRTRS = (Fptr_NL_LAPACKE_ztrtrs)dlsym(ztrrfs_obj->hModule, "LAPACKE_ztrtrs");
    ASSERT_TRUE(CTRTRS != NULL) << "failed to get the Netlib LAPACKE_ztrtrs symbol";

    CTRRFS = (Fptr_NL_LAPACKE_ztrrfs)dlsym(ztrrfs_obj->hModule, "LAPACKE_ztrrfs");
    ASSERT_TRUE(CTRRFS != NULL) << "failed to get the Netlib LAPACKE_ztrrfs symbol";

    /* Generate i/ps to trrfs API from 'trtrs' API */
    ztrrfs_obj->inforef = CTRTRS( ztrrfs_obj->matrix_layout, ztrrfs_obj->uplo,
                                  ztrrfs_obj->trans, ztrrfs_obj->diag, 
                                  ztrrfs_obj->n, ztrrfs_obj->nrhs,
                                  (const lapack_complex_double *)ztrrfs_obj->aref, 
                                  ztrrfs_obj->lda, ztrrfs_obj->xref,
                                  ztrrfs_obj->ldb );

    ztrrfs_obj->info = LAPACKE_ztrtrs( ztrrfs_obj->matrix_layout, ztrrfs_obj->uplo,
                                  ztrrfs_obj->trans, ztrrfs_obj->diag, 
                                  ztrrfs_obj->n, ztrrfs_obj->nrhs,
                                  (const lapack_complex_double *)ztrrfs_obj->a, 
                                  ztrrfs_obj->lda, ztrrfs_obj->x,
                                  ztrrfs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    ztrrfs_obj->inforef = CTRRFS( ztrrfs_obj->matrix_layout, ztrrfs_obj->uplo,
                                  ztrrfs_obj->trans, ztrrfs_obj->diag, 
                                  ztrrfs_obj->n, ztrrfs_obj->nrhs,
                                  (const lapack_complex_double *)ztrrfs_obj->aref, 
                                  ztrrfs_obj->lda,
                                  ztrrfs_obj->bref, ztrrfs_obj->ldb,
								  ztrrfs_obj->xref, ztrrfs_obj->ldx,
                                  ztrrfs_obj->ferrref,
                                  ztrrfs_obj->berrref
								  );

    ztrrfs_obj->info = LAPACKE_ztrrfs( ztrrfs_obj->matrix_layout, ztrrfs_obj->uplo,
                                  ztrrfs_obj->trans, ztrrfs_obj->diag, 
                                  ztrrfs_obj->n, ztrrfs_obj->nrhs,
                                  (const lapack_complex_double *)ztrrfs_obj->a, 
                                  ztrrfs_obj->lda,
                                  ztrrfs_obj->b, ztrrfs_obj->ldb,
								  ztrrfs_obj->x, ztrrfs_obj->ldx,
                                  ztrrfs_obj->ferr,
                                  ztrrfs_obj->berr
								  );
								  
    if( ztrrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztrrfs is wrong\n", ztrrfs_obj->info );
    }
    if( ztrrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztrrfs is wrong\n", 
        ztrrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztrrfs_obj->diff_berr =  computeDiff_d( ztrrfs_obj->nrhs, 
                ztrrfs_obj->berr, ztrrfs_obj->berrref );
                
    ztrrfs_obj->diff_ferr =  computeDiff_d( ztrrfs_obj->nrhs, 
                ztrrfs_obj->ferr, ztrrfs_obj->ferrref );

}

TEST_F(ztrrfs_test, ztrrfs1) {
    EXPECT_NEAR(0.0, ztrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ztrrfs_test, ztrrfs2) {
    EXPECT_NEAR(0.0, ztrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ztrrfs_test, ztrrfs3) {
    EXPECT_NEAR(0.0, ztrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ztrrfs_test, ztrrfs4) {
    EXPECT_NEAR(0.0, ztrrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ztrrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

