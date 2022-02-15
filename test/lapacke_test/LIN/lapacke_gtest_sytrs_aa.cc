#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define LAPACKE_TEST_VERBOSE (1)
#define sytrs_aa_free() \
       if (ipiv != NULL) free (ipiv); \
       if (bref != NULL) free (bref); \
       if (b != NULL)    free (b   ); \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (ipivref != NULL)free (ipivref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

class sytrs_aa_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

      float threshold;
   public:
      sytrs_aa_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_aa_scomplex_parameters ();
};  /* end of sytrs_aa_scomplex_parameters  class definition */


/* Constructor sytrs_aa_scomplex_parameters definition */
sytrs_aa_scomplex_parameters:: sytrs_aa_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_aa sComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
		ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
		ldb = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

     b_bufsize = ldb*nrhs;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_aa_scomplex_parameters object: malloc error.";
       sytrs_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, n, n,uplo);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_int_buffer_pair_with_constant( ipiv, ipivref, n, 0);

   } /* end of Constructor  */

sytrs_aa_scomplex_parameters:: ~sytrs_aa_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_aa_scomplex_parameters object: destructor invoked. \n");
#endif
   sytrs_aa_free();
}


//  Test fixture class definition
class csytrs_aa_test  : public  ::testing::Test {
public:
   sytrs_aa_scomplex_parameters  *csytrs_aa_obj;
   void SetUp();  
   void TearDown () { delete csytrs_aa_obj; }
};


void csytrs_aa_test::SetUp(){

    /* LAPACKE ZSYTRS_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrs_aa) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_float * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_csytrs_aa ZSYTRS_AA;

     /* LAPACKE ZSYTRF_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf_aa) ( int matrix_layout,char uplo ,lapack_int n,
                                    lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_csytrf_aa ZSYTRF_AA;


    csytrs_aa_obj = new  sytrs_aa_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    csytrs_aa_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
    idx = Circular_Increment_Index(idx);


    csytrs_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csytrs_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csytrs_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csytrs_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYTRS_AA = (Fptr_NL_LAPACKE_csytrs_aa)dlsym(csytrs_aa_obj->hModule, "LAPACKE_csytrs_aa");
    ASSERT_TRUE(ZSYTRS_AA != NULL) << "failed to get the Netlib LAPACKE_csytrs_aa symbol";

    ZSYTRF_AA = (Fptr_NL_LAPACKE_csytrf_aa)dlsym(csytrs_aa_obj->hModule,"LAPACKE_csytrf_aa");
    ASSERT_TRUE(ZSYTRF_AA != NULL) << "failed to get the Netlib LAPACKE_csytrf_aa symbol";

    /* Pre condition: need to call sytrf_aa - before calling sytrs_aa function */

    /* Compute the Netlib-Lapacke's reference o/p */
    csytrs_aa_obj->inforef = ZSYTRF_AA( csytrs_aa_obj->matrix_layout,
                                    csytrs_aa_obj->uplo, csytrs_aa_obj->n,
                                     csytrs_aa_obj->aref,
                               csytrs_aa_obj->lda, csytrs_aa_obj->ipivref);

    csytrs_aa_obj->inforef = ZSYTRS_AA( csytrs_aa_obj->matrix_layout,
                                  csytrs_aa_obj->uplo, csytrs_aa_obj->n,
                                  csytrs_aa_obj->nrhs,
                                  (const lapack_complex_float *)csytrs_aa_obj->aref,
                                  csytrs_aa_obj->lda, csytrs_aa_obj->ipivref,
                                  csytrs_aa_obj->bref, csytrs_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    csytrs_aa_obj->info = LAPACKE_csytrf_aa( csytrs_aa_obj->matrix_layout,
                                 csytrs_aa_obj->uplo, csytrs_aa_obj->n,
                                     csytrs_aa_obj->a,
                               csytrs_aa_obj->lda, csytrs_aa_obj->ipiv);

    csytrs_aa_obj->info = LAPACKE_csytrs_aa( csytrs_aa_obj->matrix_layout,
                csytrs_aa_obj->uplo, csytrs_aa_obj->n, csytrs_aa_obj->nrhs,
                                  (const lapack_complex_float *)csytrs_aa_obj->a,
                               csytrs_aa_obj->lda, csytrs_aa_obj->ipiv,
                                 csytrs_aa_obj->b, csytrs_aa_obj->ldb );


    if( csytrs_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csytrs_aa is wrong\n",
                    csytrs_aa_obj->info );
    }
    if( csytrs_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrs_aa is wrong\n",
        csytrs_aa_obj->inforef );
    }
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_aa_obj->b_bufsize,
                           csytrs_aa_obj->b, csytrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, csytrs_aa_obj->threshold);

}

TEST_F(csytrs_aa_test, csytrs_aa1) {}
TEST_F(csytrs_aa_test, csytrs_aa2) {}
TEST_F(csytrs_aa_test, csytrs_aa3) {}
TEST_F(csytrs_aa_test, csytrs_aa4) {}


/* Begin sytrs_aa_double_parameters  class definition */
class sytrs_aa_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;
	  
	  float threshold;

   public:
      sytrs_aa_double_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_aa_double_parameters ();
};  /* end of sytrs_aa_double_parameters  class definition */


/* Constructor sytrs_aa_double_parameters definition */
sytrs_aa_double_parameters:: sytrs_aa_double_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_aa Double:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
		ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
		ldb = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

     b_bufsize = ldb*nrhs;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_aa_double_parameters object: malloc error.";
       sytrs_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref, n, lda,uplo);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_int_buffer_pair_with_constant( ipiv, ipivref, n, 0);

   } /* end of Constructor  */

sytrs_aa_double_parameters:: ~sytrs_aa_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_aa_double_parameters object: destructor invoked. \n");
#endif
//   sytrs_aa_free();
       if( hModule != NULL) dlclose(hModule); 
       if(dModule != NULL) dlclose(dModule);

		if (ipiv != NULL) free (ipiv); 
       if (bref != NULL) free (bref); 
       if (b != NULL)    free (b   ); 
       if (a != NULL)    free (a   ); 
       if (aref != NULL) free (aref); 
       if (ipivref != NULL)free (ipivref); 

}


//  Test fixture class definition
class dsytrs_aa_test  : public  ::testing::Test {
public:
   sytrs_aa_double_parameters  *dsytrs_aa_obj;
   void SetUp();  
   void TearDown () { delete dsytrs_aa_obj; }
};


void dsytrs_aa_test::SetUp(){

    /* LAPACKE DSYTRS_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrs_aa) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs, const double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                                            double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dsytrs_aa DSYTRS_AA;

     /* LAPACKE DSYTRF_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf_aa) ( int matrix_layout, char uplo,
             lapack_int n, double* a, lapack_int lda, lapack_int* ipiv );

    Fptr_NL_LAPACKE_dsytrf_aa DSYTRF_AA;

    dsytrs_aa_obj = new  sytrs_aa_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    dsytrs_aa_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
    idx = Circular_Increment_Index(idx);

    dsytrs_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsytrs_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsytrs_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsytrs_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYTRS_AA = (Fptr_NL_LAPACKE_dsytrs_aa)dlsym(dsytrs_aa_obj->hModule, "LAPACKE_dsytrs_aa");
    ASSERT_TRUE(DSYTRS_AA != NULL) << "failed to get the Netlib LAPACKE_dsytrs_aa symbol";

    DSYTRF_AA = (Fptr_NL_LAPACKE_dsytrf_aa)dlsym(dsytrs_aa_obj->hModule,"LAPACKE_dsytrf_aa");
    ASSERT_TRUE(DSYTRF_AA != NULL) << "failed to get the Netlib LAPACKE_dsytrf_aa symbol";

    /* Pre condition: need to call sytrf_aa - before calling sytrs_aa function */

    /* Compute the Netlib-Lapacke's reference o/p */
    dsytrs_aa_obj->inforef = DSYTRF_AA( dsytrs_aa_obj->matrix_layout,
                            dsytrs_aa_obj->uplo, dsytrs_aa_obj->n,
                                     dsytrs_aa_obj->aref,
                      dsytrs_aa_obj->lda, dsytrs_aa_obj->ipivref);

    dsytrs_aa_obj->inforef = DSYTRS_AA( dsytrs_aa_obj->matrix_layout,
                                  dsytrs_aa_obj->uplo, dsytrs_aa_obj->n,
                                  dsytrs_aa_obj->nrhs,
                                  (const double *)dsytrs_aa_obj->aref,
                              dsytrs_aa_obj->lda, dsytrs_aa_obj->ipivref,
                                dsytrs_aa_obj->bref, dsytrs_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dsytrs_aa_obj->info = LAPACKE_dsytrf_aa( dsytrs_aa_obj->matrix_layout,
                                 dsytrs_aa_obj->uplo, dsytrs_aa_obj->n,
                                     dsytrs_aa_obj->a,
                               dsytrs_aa_obj->lda, dsytrs_aa_obj->ipiv);

    dsytrs_aa_obj->info = LAPACKE_dsytrs_aa( dsytrs_aa_obj->matrix_layout,
                dsytrs_aa_obj->uplo, dsytrs_aa_obj->n, dsytrs_aa_obj->nrhs,
                                  (const double *)dsytrs_aa_obj->a,
                               dsytrs_aa_obj->lda, dsytrs_aa_obj->ipiv,
                                 dsytrs_aa_obj->b, dsytrs_aa_obj->ldb );


    if( dsytrs_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsytrs_aa is wrong\n",
                    dsytrs_aa_obj->info );
    }
    if( dsytrs_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrs_aa is wrong\n",
        dsytrs_aa_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    double diff;
    diff =  computeDiff_d( dsytrs_aa_obj->b_bufsize,
                           dsytrs_aa_obj->b, dsytrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, dsytrs_aa_obj->threshold);

}

TEST_F(dsytrs_aa_test, dsytrs_aa1) {}
TEST_F(dsytrs_aa_test, dsytrs_aa2) {}
TEST_F(dsytrs_aa_test, dsytrs_aa3) {}
TEST_F(dsytrs_aa_test, dsytrs_aa4) {}

/* Begin sytrs_aa_float_parameters  class definition */
class sytrs_aa_float_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;
	  
	  float threshold;

   public:
      sytrs_aa_float_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_aa_float_parameters ();
};  /* end of sytrs_aa_float_parameters  class definition */


/* Constructor sytrs_aa_float_parameters definition */
sytrs_aa_float_parameters:: sytrs_aa_float_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_aa Double:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
		ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
		ldb = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

     b_bufsize = ldb*nrhs;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_aa_float_parameters object: malloc error.";
       sytrs_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, n, lda,uplo);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_int_buffer_pair_with_constant( ipiv, ipivref, n, 0);

   } /* end of Constructor  */

sytrs_aa_float_parameters:: ~sytrs_aa_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_aa_float_parameters object: destructor invoked. \n");
#endif
//   sytrs_aa_free();
       if( hModule != NULL) dlclose(hModule); 
       if(dModule != NULL) dlclose(dModule);

		if (ipiv != NULL) free (ipiv); 
       if (bref != NULL) free (bref); 
       if (b != NULL)    free (b   ); 
       if (a != NULL)    free (a   ); 
       if (aref != NULL) free (aref); 
       if (ipivref != NULL)free (ipivref); 

}


//  Test fixture class definition
class ssytrs_aa_test  : public  ::testing::Test {
public:
   sytrs_aa_float_parameters  *ssytrs_aa_obj;
   void SetUp();  
   void TearDown () { delete ssytrs_aa_obj; }
};


void ssytrs_aa_test::SetUp(){

    /* LAPACKE DSYTRS_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrs_aa) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs, const float * a,
                                  lapack_int lda, const lapack_int * ipiv,
                                            float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_ssytrs_aa DSYTRS_AA;

     /* LAPACKE DSYTRF_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf_aa) ( int matrix_layout, char uplo,
             lapack_int n, float* a, lapack_int lda, lapack_int* ipiv );

    Fptr_NL_LAPACKE_ssytrf_aa DSYTRF_AA;

    ssytrs_aa_obj = new  sytrs_aa_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    ssytrs_aa_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
    idx = Circular_Increment_Index(idx);

    ssytrs_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssytrs_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssytrs_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssytrs_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYTRS_AA = (Fptr_NL_LAPACKE_ssytrs_aa)dlsym(ssytrs_aa_obj->hModule, "LAPACKE_ssytrs_aa");
    ASSERT_TRUE(DSYTRS_AA != NULL) << "failed to get the Netlib LAPACKE_ssytrs_aa symbol";

    DSYTRF_AA = (Fptr_NL_LAPACKE_ssytrf_aa)dlsym(ssytrs_aa_obj->hModule,"LAPACKE_ssytrf_aa");
    ASSERT_TRUE(DSYTRF_AA != NULL) << "failed to get the Netlib LAPACKE_ssytrf_aa symbol";

    /* Pre condition: need to call sytrf_aa - before calling sytrs_aa function */

    /* Compute the Netlib-Lapacke's reference o/p */
    ssytrs_aa_obj->inforef = DSYTRF_AA( ssytrs_aa_obj->matrix_layout,
                            ssytrs_aa_obj->uplo, ssytrs_aa_obj->n,
                                     ssytrs_aa_obj->aref,
                      ssytrs_aa_obj->lda, ssytrs_aa_obj->ipivref);

    ssytrs_aa_obj->inforef = DSYTRS_AA( ssytrs_aa_obj->matrix_layout,
                                  ssytrs_aa_obj->uplo, ssytrs_aa_obj->n,
                                  ssytrs_aa_obj->nrhs,
                                  (const float *)ssytrs_aa_obj->aref,
                              ssytrs_aa_obj->lda, ssytrs_aa_obj->ipivref,
                                ssytrs_aa_obj->bref, ssytrs_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    ssytrs_aa_obj->info = LAPACKE_ssytrf_aa( ssytrs_aa_obj->matrix_layout,
                                 ssytrs_aa_obj->uplo, ssytrs_aa_obj->n,
                                     ssytrs_aa_obj->a,
                               ssytrs_aa_obj->lda, ssytrs_aa_obj->ipiv);

    ssytrs_aa_obj->info = LAPACKE_ssytrs_aa( ssytrs_aa_obj->matrix_layout,
                ssytrs_aa_obj->uplo, ssytrs_aa_obj->n, ssytrs_aa_obj->nrhs,
                                  (const float *)ssytrs_aa_obj->a,
                               ssytrs_aa_obj->lda, ssytrs_aa_obj->ipiv,
                                 ssytrs_aa_obj->b, ssytrs_aa_obj->ldb );


    if( ssytrs_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssytrs_aa is wrong\n",
                    ssytrs_aa_obj->info );
    }
    if( ssytrs_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrs_aa is wrong\n",
        ssytrs_aa_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    float diff;
    diff =  computeDiff_s( ssytrs_aa_obj->b_bufsize,
                           ssytrs_aa_obj->b, ssytrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, ssytrs_aa_obj->threshold);

}

TEST_F(ssytrs_aa_test, ssytrs_aa1) {}
TEST_F(ssytrs_aa_test, ssytrs_aa2) {}
TEST_F(ssytrs_aa_test, ssytrs_aa3) {}
TEST_F(ssytrs_aa_test, ssytrs_aa4) {}


/* Begin sytrs_aa_dcomplex_parameters  class definition */
class sytrs_aa_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

      float threshold;
   public:
      sytrs_aa_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_aa_dcomplex_parameters ();
};  /* end of sytrs_aa_dcomplex_parameters  class definition */


/* Constructor sytrs_aa_dcomplex_parameters definition */
sytrs_aa_dcomplex_parameters:: sytrs_aa_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_aa DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
		ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
		ldb = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

     b_bufsize = ldb*nrhs;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_aa_dcomplex_parameters object: malloc error.";
       sytrs_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, n, lda,uplo);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_int_buffer_pair_with_constant( ipiv, ipivref, n, 0);

   } /* end of Constructor  */

sytrs_aa_dcomplex_parameters:: ~sytrs_aa_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_aa_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytrs_aa_free();
}


//  Test fixture class definition
class zsytrs_aa_test  : public  ::testing::Test {
public:
   sytrs_aa_dcomplex_parameters  *zsytrs_aa_obj;
   void SetUp();  
   void TearDown () { delete zsytrs_aa_obj; }
};


void zsytrs_aa_test::SetUp(){

    /* LAPACKE ZSYTRS_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrs_aa) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zsytrs_aa ZSYTRS_AA;

     /* LAPACKE ZSYTRF_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf_aa) ( int matrix_layout,char uplo ,lapack_int n,
                                    lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zsytrf_aa ZSYTRF_AA;


    zsytrs_aa_obj = new  sytrs_aa_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    zsytrs_aa_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
    idx = Circular_Increment_Index(idx);


    zsytrs_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsytrs_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsytrs_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsytrs_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYTRS_AA = (Fptr_NL_LAPACKE_zsytrs_aa)dlsym(zsytrs_aa_obj->hModule, "LAPACKE_zsytrs_aa");
    ASSERT_TRUE(ZSYTRS_AA != NULL) << "failed to get the Netlib LAPACKE_zsytrs_aa symbol";

    ZSYTRF_AA = (Fptr_NL_LAPACKE_zsytrf_aa)dlsym(zsytrs_aa_obj->hModule,"LAPACKE_zsytrf_aa");
    ASSERT_TRUE(ZSYTRF_AA != NULL) << "failed to get the Netlib LAPACKE_zsytrf_aa symbol";

    /* Pre condition: need to call sytrf_aa - before calling sytrs_aa function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zsytrs_aa_obj->inforef = ZSYTRF_AA( zsytrs_aa_obj->matrix_layout,
                                    zsytrs_aa_obj->uplo, zsytrs_aa_obj->n,
                                     zsytrs_aa_obj->aref,
                               zsytrs_aa_obj->lda, zsytrs_aa_obj->ipivref);

    zsytrs_aa_obj->inforef = ZSYTRS_AA( zsytrs_aa_obj->matrix_layout,
                                  zsytrs_aa_obj->uplo, zsytrs_aa_obj->n,
                                  zsytrs_aa_obj->nrhs,
                                  (const lapack_complex_double *)zsytrs_aa_obj->aref,
                                  zsytrs_aa_obj->lda, zsytrs_aa_obj->ipivref,
                                  zsytrs_aa_obj->bref, zsytrs_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zsytrs_aa_obj->info = LAPACKE_zsytrf_aa( zsytrs_aa_obj->matrix_layout,
                                 zsytrs_aa_obj->uplo, zsytrs_aa_obj->n,
                                     zsytrs_aa_obj->a,
                               zsytrs_aa_obj->lda, zsytrs_aa_obj->ipiv);

    zsytrs_aa_obj->info = LAPACKE_zsytrs_aa( zsytrs_aa_obj->matrix_layout,
                zsytrs_aa_obj->uplo, zsytrs_aa_obj->n, zsytrs_aa_obj->nrhs,
                                  (const lapack_complex_double *)zsytrs_aa_obj->a,
                               zsytrs_aa_obj->lda, zsytrs_aa_obj->ipiv,
                                 zsytrs_aa_obj->b, zsytrs_aa_obj->ldb );


    if( zsytrs_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsytrs_aa is wrong\n",
                    zsytrs_aa_obj->info );
    }
    if( zsytrs_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrs_aa is wrong\n",
        zsytrs_aa_obj->inforef );
    }
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_aa_obj->b_bufsize,
                           zsytrs_aa_obj->b, zsytrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, zsytrs_aa_obj->threshold);

}

TEST_F(zsytrs_aa_test, zsytrs_aa1) {}
TEST_F(zsytrs_aa_test, zsytrs_aa2) {}
TEST_F(zsytrs_aa_test, zsytrs_aa3) {}
TEST_F(zsytrs_aa_test, zsytrs_aa4) {}

