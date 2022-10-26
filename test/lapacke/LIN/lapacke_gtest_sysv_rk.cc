#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define sysv_rk_free() \
       if (ipiv != NULL) free (ipiv); \
       if (bref != NULL) free (bref); \
       if (e != NULL)    free (e   ); \
       if (eref != NULL) free (eref); \
       if (b != NULL)    free (b   ); \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (ipivref != NULL)free (ipivref); \
       if( hModule != NULL) dlclose(hModule)  
//; \
//       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin sysv_rk_double_parameters  class definition */
class sysv_rk_double_parameters{
   public:
      int b_bufsize;
	  double diff; // capture difference between Netlib, libflame o/p
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array 'a' contains the matrix A
      void *hModule;//, *dModule;
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Output parameters */
      /*  superdiagonal (or subdiagonal) elements of the symmetric block diagonal 
           matrix D with 1-by-1 or 2-by-2 diagonal blocks.  */
      double *e,*eref;
      lapack_int *ipiv, *ipivref; // The pivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      sysv_rk_double_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_rk_double_parameters ();
};  /* end of sysv_rk_double_parameters  class definition */


/* Constructor sysv_rk_double_parameters definition */
sysv_rk_double_parameters:: sysv_rk_double_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
 //   dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sysv_rk Double:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sysv_rk_double_parameters object: malloc error.";
       sysv_rk_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_rk_double_parameters:: ~sysv_rk_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_rk_double_parameters object: destructor invoked. \n");
#endif
   sysv_rk_free();
}


//  Test fixture class definition
class dsysv_rk_test  : public  ::testing::Test {
public:
   sysv_rk_double_parameters  *dsysv_rk_obj;
   void SetUp();  
   void TearDown () { delete dsysv_rk_obj; }
};


void dsysv_rk_test::SetUp(){

    /* LAPACKE DSYSV_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_dsysv_rk) ( int matrix_layout, char uplo,
                                 lapack_int n, lapack_int nrhs,  double * a,
                             lapack_int lda,  double *e,  lapack_int * ipiv,
                                            double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dsysv_rk DSYSV_RK;

    dsysv_rk_obj = new  sysv_rk_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    //dsysv_rk_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsysv_rk_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    //ASSERT_TRUE(dsysv_rk_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsysv_rk_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYSV_RK = (Fptr_NL_LAPACKE_dsysv_rk)dlsym(dsysv_rk_obj->hModule, "LAPACKE_dsysv_rk");
    ASSERT_TRUE(DSYSV_RK != NULL) << "failed to get the Netlib LAPACKE_dsysv_rk symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    dsysv_rk_obj->inforef = DSYSV_RK( dsysv_rk_obj->matrix_layout,
                                  dsysv_rk_obj->uplo, dsysv_rk_obj->n,
                                  dsysv_rk_obj->nrhs,dsysv_rk_obj->aref,
                                  dsysv_rk_obj->lda,
                               dsysv_rk_obj->eref, dsysv_rk_obj->ipivref,
                                dsysv_rk_obj->bref, dsysv_rk_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dsysv_rk_obj->info = LAPACKE_dsysv_rk( dsysv_rk_obj->matrix_layout,
                dsysv_rk_obj->uplo, dsysv_rk_obj->n, dsysv_rk_obj->nrhs,
                                  dsysv_rk_obj->a, dsysv_rk_obj->lda, 
								  dsysv_rk_obj->e, dsysv_rk_obj->ipiv,
                                 dsysv_rk_obj->b, dsysv_rk_obj->ldb );


    if( dsysv_rk_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsysv_rk is wrong\n",
                    dsysv_rk_obj->info );
    }
    if( dsysv_rk_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsysv_rk is wrong\n",
        dsysv_rk_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dsysv_rk_obj->diff =  computeDiff_d( dsysv_rk_obj->n,
                           dsysv_rk_obj->e, dsysv_rk_obj->eref );
}

TEST_F(dsysv_rk_test, dsysv_rk1) {
    EXPECT_NEAR(0.0, dsysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsysv_rk_test, dsysv_rk2) {
    EXPECT_NEAR(0.0, dsysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsysv_rk_test, dsysv_rk3) {
    EXPECT_NEAR(0.0, dsysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsysv_rk_test, dsysv_rk4) {
    EXPECT_NEAR(0.0, dsysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sysv_rk_float_parameters  class definition */
class sysv_rk_float_parameters{

   public:
      int b_bufsize;
	  float diff; // captures difference between netlib and libflame o/ps
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array 'a' contains the matrix A
      void *hModule, *dModule;
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Output parameters */
      lapack_int *ipiv, *ipivref; // The pivot indices
      /*  superdiagonal (or subdiagonal) elements of the symmetric block diagonal 
           matrix D with 1-by-1 or 2-by-2 diagonal blocks.  */
      float *e,*eref;

      /* Return Values */
      lapack_int info, inforef;

   public:
      sysv_rk_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~sysv_rk_float_parameters ();
};  /* end of sysv_rk_float_parameters  class definition */


/* Constructor sysv_rk_float_parameters definition */
sysv_rk_float_parameters:: sysv_rk_float_parameters ( int matrix_layout_i,
                       char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sysv_rk float:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif
    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (lda*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sysv_rk_double_parameters object: malloc error.";
       sysv_rk_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

sysv_rk_float_parameters:: ~sysv_rk_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_rk_float_parameters object: destructor invoked. \n");
#endif
   sysv_rk_free();
}


//  Test fixture class definition
class ssysv_rk_test  : public  ::testing::Test {
public:
   sysv_rk_float_parameters  *ssysv_rk_obj;
   void SetUp();  
   void TearDown () { delete ssysv_rk_obj; }
};


void ssysv_rk_test::SetUp(){

    /* LAPACKE SSYSV_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_ssysv_rk) ( int matrix_layout, char uplo,
                                 lapack_int n, lapack_int nrhs,  float * a,
                               lapack_int lda, float *e, lapack_int * ipiv,
                                              float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_ssysv_rk SSYSV_RK;

    ssysv_rk_obj = new  sysv_rk_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    ssysv_rk_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssysv_rk_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssysv_rk_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssysv_rk_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYSV_RK = (Fptr_NL_LAPACKE_ssysv_rk)dlsym(ssysv_rk_obj->hModule, "LAPACKE_ssysv_rk");
    ASSERT_TRUE(SSYSV_RK != NULL) << "failed to get the Netlib LAPACKE_ssysv_rk symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    ssysv_rk_obj->inforef = SSYSV_RK( ssysv_rk_obj->matrix_layout,
                                  ssysv_rk_obj->uplo, ssysv_rk_obj->n,
                               ssysv_rk_obj->nrhs, ssysv_rk_obj->aref,
                                 ssysv_rk_obj->lda,ssysv_rk_obj->eref,
                                                ssysv_rk_obj->ipivref,
                                  ssysv_rk_obj->bref, ssysv_rk_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssysv_rk_obj->info = LAPACKE_ssysv_rk( ssysv_rk_obj->matrix_layout,
                ssysv_rk_obj->uplo, ssysv_rk_obj->n, ssysv_rk_obj->nrhs,
                                     ssysv_rk_obj->a, ssysv_rk_obj->lda,
                                    ssysv_rk_obj->e, ssysv_rk_obj->ipiv,
                                   ssysv_rk_obj->b, ssysv_rk_obj->ldb );
    if( ssysv_rk_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssysv_rk is wrong\n",
                    ssysv_rk_obj->info );
    }
    if( ssysv_rk_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssysv_rk is wrong\n",
        ssysv_rk_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ssysv_rk_obj->diff =  computeDiff_s( ssysv_rk_obj->b_bufsize,
                           ssysv_rk_obj->b, ssysv_rk_obj->bref );
}

TEST_F(ssysv_rk_test, ssysv_rk1) {
    EXPECT_NEAR(0.0, ssysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_rk_test, ssysv_rk2) {
    EXPECT_NEAR(0.0, ssysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_rk_test, ssysv_rk3) {
    EXPECT_NEAR(0.0, ssysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_rk_test, ssysv_rk4) {
    EXPECT_NEAR(0.0, ssysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sysv_rk_scomplex_parameters  class definition */
class sysv_rk_scomplex_parameters{
   public:
	  float diff; // captures difference between netlib and libflame o/ps
      int b_bufsize;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array 'a' contains the matrix A
      void *hModule, *dModule;
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Output parameters */
      /*  superdiagonal (or subdiagonal) elements of the symmetric block diagonal 
           matrix D with 1-by-1 or 2-by-2 diagonal blocks.  */
      lapack_complex_float *e,*eref;
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      sysv_rk_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_rk_scomplex_parameters ();
};  /* end of sysv_rk_scomplex_parameters  class definition */


/* Constructor sysv_rk_scomplex_parameters definition */
sysv_rk_scomplex_parameters:: sysv_rk_scomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sysv_rk scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sysv_rk_scomplex_parameters object: malloc error.";
       sysv_rk_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_rk_scomplex_parameters:: ~sysv_rk_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_rk_scomplex_parameters object: destructor invoked. \n");
#endif
   sysv_rk_free();
}


//  Test fixture class definition
class csysv_rk_test  : public  ::testing::Test {
public:
   sysv_rk_scomplex_parameters  *csysv_rk_obj;
   void SetUp();  
   void TearDown () { delete csysv_rk_obj; }
};


void csysv_rk_test::SetUp(){

    /* LAPACKE CSYSV_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_csysv_rk) ( int matrix_layout, char uplo,
                                             lapack_int n, lapack_int nrhs,
                                  lapack_complex_float * a, lapack_int lda, 
                               lapack_complex_float * e, lapack_int * ipiv,
                               lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_csysv_rk CSYSV_RK;

    csysv_rk_obj = new  sysv_rk_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    csysv_rk_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csysv_rk_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csysv_rk_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csysv_rk_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSYSV_RK = (Fptr_NL_LAPACKE_csysv_rk)dlsym(csysv_rk_obj->hModule, "LAPACKE_csysv_rk");
    ASSERT_TRUE(CSYSV_RK != NULL) << "failed to get the Netlib LAPACKE_csysv_rk symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    csysv_rk_obj->inforef = CSYSV_RK( csysv_rk_obj->matrix_layout,
                                  csysv_rk_obj->uplo, csysv_rk_obj->n,
                                  csysv_rk_obj->nrhs, csysv_rk_obj->aref,
                                  csysv_rk_obj->lda, csysv_rk_obj->eref,
								  csysv_rk_obj->ipivref,
                                  csysv_rk_obj->bref, csysv_rk_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    csysv_rk_obj->info = LAPACKE_csysv_rk( csysv_rk_obj->matrix_layout,
                csysv_rk_obj->uplo, csysv_rk_obj->n, csysv_rk_obj->nrhs,
                                  csysv_rk_obj->a, csysv_rk_obj->lda, 
								  csysv_rk_obj->e, csysv_rk_obj->ipiv,
                                 csysv_rk_obj->b, csysv_rk_obj->ldb );


    if( csysv_rk_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csysv_rk is wrong\n",
                    csysv_rk_obj->info );
    }
    if( csysv_rk_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csysv_rk is wrong\n",
        csysv_rk_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    csysv_rk_obj->diff =  computeDiff_c( csysv_rk_obj->b_bufsize,
                           csysv_rk_obj->b, csysv_rk_obj->bref );
}

TEST_F(csysv_rk_test, csysv_rk1) {
    EXPECT_NEAR(0.0, csysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csysv_rk_test, csysv_rk2) {
    EXPECT_NEAR(0.0, csysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csysv_rk_test, csysv_rk3) {
    EXPECT_NEAR(0.0, csysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csysv_rk_test, csysv_rk4) {
    EXPECT_NEAR(0.0, csysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sysv_rk_dcomplex_parameters  class definition */
class sysv_rk_dcomplex_parameters{
   public:
	  double diff; // captures difference between netlib and libflame o/ps
      int b_bufsize;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array 'a' contains the matrix A
      void *hModule, *dModule;
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Output parameters */
      /*  superdiagonal (or subdiagonal) elements of the symmetric block diagonal 
           matrix D with 1-by-1 or 2-by-2 diagonal blocks.  */
      lapack_complex_double *e,*eref;
      lapack_int *ipiv,*ipivref; // The ipivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      sysv_rk_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_rk_dcomplex_parameters ();
};  /* end of sysv_rk_dcomplex_parameters  class definition */


/* Constructor sysv_rk_dcomplex_parameters definition */
sysv_rk_dcomplex_parameters:: sysv_rk_dcomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sysv_rk DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sysv_rk_dcomplex_parameters object: malloc error.";
       sysv_rk_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_rk_dcomplex_parameters:: ~sysv_rk_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_rk_dcomplex_parameters object: destructor invoked. \n");
#endif
   sysv_rk_free();
}


//  Test fixture class definition
class zsysv_rk_test  : public  ::testing::Test {
public:
   sysv_rk_dcomplex_parameters  *zsysv_rk_obj;
   void SetUp();  
   void TearDown () { delete zsysv_rk_obj; }
};


void zsysv_rk_test::SetUp(){

    /* LAPACKE ZSYSV_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_zsysv_rk) ( int matrix_layout, char uplo,
                                             lapack_int n, lapack_int nrhs,
								  lapack_complex_double *a, lapack_int lda,  
							   lapack_complex_double *e, lapack_int * ipiv,
								lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zsysv_rk ZSYSV_RK;

    zsysv_rk_obj = new  sysv_rk_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zsysv_rk_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsysv_rk_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsysv_rk_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsysv_rk_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYSV_RK = (Fptr_NL_LAPACKE_zsysv_rk)dlsym(zsysv_rk_obj->hModule, "LAPACKE_zsysv_rk");
    ASSERT_TRUE(ZSYSV_RK != NULL) << "failed to get the Netlib LAPACKE_zsysv_rk symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    zsysv_rk_obj->inforef = ZSYSV_RK( zsysv_rk_obj->matrix_layout,
                                  zsysv_rk_obj->uplo, zsysv_rk_obj->n,
                                  zsysv_rk_obj->nrhs, zsysv_rk_obj->aref,
                                  zsysv_rk_obj->lda, zsysv_rk_obj->eref,
								  zsysv_rk_obj->ipivref,
                                  zsysv_rk_obj->bref, zsysv_rk_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zsysv_rk_obj->info = LAPACKE_zsysv_rk( zsysv_rk_obj->matrix_layout,
                zsysv_rk_obj->uplo, zsysv_rk_obj->n, zsysv_rk_obj->nrhs,
                                  zsysv_rk_obj->a, zsysv_rk_obj->lda,
                                zsysv_rk_obj->e,zsysv_rk_obj->ipiv,
                                 zsysv_rk_obj->b, zsysv_rk_obj->ldb );


    if( zsysv_rk_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsysv_rk is wrong\n",
                    zsysv_rk_obj->info );
    }
    if( zsysv_rk_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsysv_rk is wrong\n",
        zsysv_rk_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zsysv_rk_obj->diff =  computeDiff_z( zsysv_rk_obj->b_bufsize,
                           zsysv_rk_obj->b, zsysv_rk_obj->bref );
}

TEST_F(zsysv_rk_test, zsysv_rk1) {
    EXPECT_NEAR(0.0, zsysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsysv_rk_test, zsysv_rk2) {
    EXPECT_NEAR(0.0, zsysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsysv_rk_test, zsysv_rk3) {
    EXPECT_NEAR(0.0, zsysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsysv_rk_test, zsysv_rk4) {
    EXPECT_NEAR(0.0, zsysv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

