#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

// #define LAPACKE_TEST_VERBOSE (1)
#define gbsv_free() \
       if (ipiv != NULL) free (ipiv); \
       if (ipivref != NULL)free (ipivref); \
       if (b != NULL)    free (b   ); \
       if (bref != NULL) free (bref); \
       if (ab != NULL)    free (ab  ); \
       if (abref != NULL) free (abref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin gbsv_double_parameters  class definition */
class gbsv_double_parameters{
   public:
      int b_bufsize;
      double diff; // to capture the netlib, libflame o/ps

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kl;// The number of subdiagonals within the band of A
      lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double *ab, *abref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gbsv_double_parameters ( int matrix_layout_i, 
                                lapack_int n_i, lapack_int ldab_i,
                                 lapack_int kl_i, lapack_int ku_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gbsv_double_parameters ();
};  /* end of gbsv_double_parameters  class definition */


/* Constructor gbsv_double_parameters definition */
gbsv_double_parameters:: gbsv_double_parameters ( int matrix_layout_i,
                         lapack_int n_i, lapack_int ldab_i,
                                       lapack_int kl_i, lapack_int ku_i,
                                 lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    kl = kl_i;
    ku = ku_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbsv Double:  n: %d,  lda: %d ldb: %d nrhs: %d \n",
             n, ldab, ldb, nrhs);
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
    lapacke_gtest_alloc_double_buffer_pair( &ab, &abref, (n*ldab));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "gbsv_double_parameters object: malloc error.";
       gbsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

gbsv_double_parameters:: ~gbsv_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbsv_double_parameters object: destructor invoked. \n");
#endif
   gbsv_free();
}


//  Test fixture class definition
class dgbsv_test  : public  ::testing::Test {
public:
   gbsv_double_parameters  *dgbsv_obj;
   void SetUp();  
   void TearDown () { delete dgbsv_obj; }
};


void dgbsv_test::SetUp(){

    /* LAPACKE DGBSV prototype */
    typedef int (*Fptr_NL_LAPACKE_dgbsv) ( int matrix_layout, 
                               lapack_int n, lapack_int kl, lapack_int ku,
                                        lapack_int nrhs,  double * ab,
                                  lapack_int ldab,  lapack_int * ipiv,
                                            double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dgbsv DGBSV;

    dgbsv_obj = new  gbsv_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kl,
                                         lin_solver_paramslist[idx].ku,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    dgbsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgbsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgbsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgbsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGBSV = (Fptr_NL_LAPACKE_dgbsv)dlsym(dgbsv_obj->hModule, "LAPACKE_dgbsv");
    ASSERT_TRUE(DGBSV != NULL) << "failed to get the Netlib LAPACKE_dgbsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */

    dgbsv_obj->inforef = DGBSV( dgbsv_obj->matrix_layout,
                                   dgbsv_obj->n, dgbsv_obj->kl, dgbsv_obj->ku,
                                  dgbsv_obj->nrhs, dgbsv_obj->abref,
                                  dgbsv_obj->ldab, dgbsv_obj->ipivref,
                                  dgbsv_obj->bref, dgbsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */

    dgbsv_obj->info = LAPACKE_dgbsv( dgbsv_obj->matrix_layout,
                                 dgbsv_obj->n, dgbsv_obj->kl,
                                 dgbsv_obj->ku, dgbsv_obj->nrhs,
                                  ( double *)dgbsv_obj->ab,
                               dgbsv_obj->ldab, dgbsv_obj->ipiv,
                                 dgbsv_obj->b, dgbsv_obj->ldb );


    if( dgbsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgbsv is wrong\n",
                    dgbsv_obj->info );
    }
    if( dgbsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgbsv is wrong\n",
        dgbsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dgbsv_obj->diff =  computeDiff_d( dgbsv_obj->b_bufsize,
                           dgbsv_obj->b, dgbsv_obj->bref );
}

TEST_F(dgbsv_test, dgbsv1) {
    EXPECT_NEAR(0.0, dgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbsv_test, dgbsv2) {
    EXPECT_NEAR(0.0, dgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbsv_test, dgbsv3) {
    EXPECT_NEAR(0.0, dgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbsv_test, dgbsv4) {
    EXPECT_NEAR(0.0, dgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin gbsv_float_parameters  class definition */
class gbsv_float_parameters{

   public:
      int b_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kl;// The number of subdiagonals within the band of A
      lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float *ab, *abref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gbsv_float_parameters (int matrix_layout_i, 
                               lapack_int n_i, lapack_int ldab_i,
                                lapack_int kl_i, lapack_int ku_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~gbsv_float_parameters ();
};  /* end of gbsv_float_parameters  class definition */


/* Constructor gbsv_float_parameters definition */
gbsv_float_parameters:: gbsv_float_parameters ( int matrix_layout_i,
                                    lapack_int n_i, lapack_int ldab_i,
                                     lapack_int kl_i, lapack_int ku_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    kl = kl_i;
    ku = ku_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbsv float:  n: %d,  lda: %d ldb: %d nrhs: %d \n",
             n, ldab, ldb, nrhs);
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
    lapacke_gtest_alloc_float_buffer_pair( &ab, &abref, (ldab*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "gbsv_double_parameters object: malloc error.";
       gbsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( ab, abref, ldab*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

gbsv_float_parameters:: ~gbsv_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbsv_float_parameters object: destructor invoked. \n");
#endif
   gbsv_free();
}


//  Test fixture class definition
class sgbsv_test  : public  ::testing::Test {
public:
   gbsv_float_parameters  *sgbsv_obj;
   void SetUp();  
   void TearDown () { delete sgbsv_obj; }
};


void sgbsv_test::SetUp(){

    /* LAPACKE SGBSV prototype */
    typedef int (*Fptr_NL_LAPACKE_sgbsv) ( int matrix_layout, 
                               lapack_int n, lapack_int kl, lapack_int ku, 
                                         lapack_int nrhs,  float *ab,
                                 lapack_int ldab,  lapack_int * ipiv,
                                            float *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_sgbsv SGBSV;

    sgbsv_obj = new  gbsv_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kl,
                                         lin_solver_paramslist[idx].ku,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    sgbsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgbsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgbsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgbsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGBSV = (Fptr_NL_LAPACKE_sgbsv)dlsym(sgbsv_obj->hModule, "LAPACKE_sgbsv");
    ASSERT_TRUE(SGBSV != NULL) << "failed to get the Netlib LAPACKE_sgbsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    sgbsv_obj->inforef = SGBSV( sgbsv_obj->matrix_layout,
                                   sgbsv_obj->n, sgbsv_obj->kl,
                                   sgbsv_obj->ku, sgbsv_obj->nrhs,
                                  ( float *)sgbsv_obj->abref,
                                  sgbsv_obj->ldab, sgbsv_obj->ipivref,
                                  sgbsv_obj->bref, sgbsv_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API  */
    sgbsv_obj->info = LAPACKE_sgbsv( sgbsv_obj->matrix_layout,
                                 sgbsv_obj->n, sgbsv_obj->kl,
                                sgbsv_obj->ku,  sgbsv_obj->nrhs,
                                 ( float *)sgbsv_obj->ab,
                               sgbsv_obj->ldab, sgbsv_obj->ipiv,
                                 sgbsv_obj->b, sgbsv_obj->ldb );
    if( sgbsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgbsv is wrong\n",
                    sgbsv_obj->info );
    }
    if( sgbsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgbsv is wrong\n",
        sgbsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    sgbsv_obj->diff =  computeDiff_s( sgbsv_obj->b_bufsize,
                           sgbsv_obj->b, sgbsv_obj->bref );
}

TEST_F(sgbsv_test, sgbsv1) {
    EXPECT_NEAR(0.0, sgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbsv_test, sgbsv2) {
    EXPECT_NEAR(0.0, sgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbsv_test, sgbsv3) {
    EXPECT_NEAR(0.0, sgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbsv_test, sgbsv4) {
    EXPECT_NEAR(0.0, sgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin gbsv_scomplex_parameters  class definition */
class gbsv_scomplex_parameters{

   public:
      int b_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kl;// The number of subdiagonals within the band of A
      lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *ab, *abref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gbsv_scomplex_parameters ( int matrix_layout_i, 
                                lapack_int n_i, lapack_int ldab_i,
                                 lapack_int kl_i, lapack_int ku_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gbsv_scomplex_parameters ();
};  /* end of gbsv_scomplex_parameters  class definition */


/* Constructor gbsv_scomplex_parameters definition */
gbsv_scomplex_parameters:: gbsv_scomplex_parameters ( int matrix_layout_i,
                                        lapack_int n_i, lapack_int ldab_i,
                                       lapack_int kl_i, lapack_int ku_i,
                                    lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    kl = kl_i;
    ku = ku_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbsv scomplex:  n: %d,  lda: %d ldb: %d nrhs: %d \n",
             n, ldab, ldb, nrhs);
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
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &ab, &abref, (n*ldab));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "gbsv_scomplex_parameters object: malloc error.";
       gbsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

gbsv_scomplex_parameters:: ~gbsv_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbsv_scomplex_parameters object: destructor invoked. \n");
#endif
   gbsv_free();
}


//  Test fixture class definition
class cgbsv_test  : public  ::testing::Test {
public:
   gbsv_scomplex_parameters  *cgbsv_obj;
   void SetUp();  
   void TearDown () { delete cgbsv_obj; }
};


void cgbsv_test::SetUp(){

    /* LAPACKE CGBSV prototype */
    typedef int (*Fptr_NL_LAPACKE_cgbsv) ( int matrix_layout, 
                                lapack_int n, lapack_int kl, lapack_int ku, 
                                lapack_int nrhs,  lapack_complex_float * ab,
                                lapack_int ldab,  lapack_int * ipiv,
                                lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_cgbsv CGBSV;


    cgbsv_obj = new  gbsv_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kl,
                                         lin_solver_paramslist[idx].ku,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    cgbsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgbsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgbsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgbsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGBSV = (Fptr_NL_LAPACKE_cgbsv)dlsym(cgbsv_obj->hModule, "LAPACKE_cgbsv");
    ASSERT_TRUE(CGBSV != NULL) << "failed to get the Netlib LAPACKE_cgbsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    cgbsv_obj->inforef = CGBSV( cgbsv_obj->matrix_layout,
                               cgbsv_obj->n, cgbsv_obj->kl,
                                cgbsv_obj->ku, cgbsv_obj->nrhs,
                            ( lapack_complex_float *)cgbsv_obj->abref,
                                  cgbsv_obj->ldab, cgbsv_obj->ipivref,
                                  cgbsv_obj->bref, cgbsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    cgbsv_obj->info = LAPACKE_cgbsv( cgbsv_obj->matrix_layout,
                                    cgbsv_obj->n, cgbsv_obj->kl,
                                    cgbsv_obj->ku,cgbsv_obj->nrhs,
                            ( lapack_complex_float *)cgbsv_obj->ab,
                                  cgbsv_obj->ldab, cgbsv_obj->ipiv,
                                    cgbsv_obj->b, cgbsv_obj->ldb );

    if( cgbsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgbsv is wrong\n",
                    cgbsv_obj->info );
    }
    if( cgbsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgbsv is wrong\n",
        cgbsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cgbsv_obj->diff =  computeDiff_c( cgbsv_obj->b_bufsize,
                           cgbsv_obj->b, cgbsv_obj->bref );
}

TEST_F(cgbsv_test, cgbsv1) {
    EXPECT_NEAR(0.0, cgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbsv_test, cgbsv2) {
    EXPECT_NEAR(0.0, cgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbsv_test, cgbsv3) {
    EXPECT_NEAR(0.0, cgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbsv_test, cgbsv4) {
    EXPECT_NEAR(0.0, cgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin gbsv_dcomplex_parameters  class definition */
class gbsv_dcomplex_parameters{
 
 public:
      double diff; // to capture the netlib, libflame o/ps
      int b_bufsize;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kl;// The number of subdiagonals within the band of A
      lapack_int ku; // The number of superdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *ab, *abref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gbsv_dcomplex_parameters ( int matrix_layout_i, 
                                lapack_int n_i, lapack_int ldab_i,
                                 lapack_int kl_i, lapack_int ku_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gbsv_dcomplex_parameters ();
};  /* end of gbsv_dcomplex_parameters  class definition */


/* Constructor gbsv_dcomplex_parameters definition */
gbsv_dcomplex_parameters:: gbsv_dcomplex_parameters ( int matrix_layout_i,
                                        lapack_int n_i, lapack_int ldab_i,
                                        lapack_int kl_i, lapack_int ku_i,
                                    lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    kl = kl_i;
    ku = ku_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gbsv DComplex:  n: %d,  lda: %d ldb: %d nrhs: %d \n",
             n, ldab, ldb, nrhs);
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
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &ab, &abref, (n*ldab));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "gbsv_dcomplex_parameters object: malloc error.";
       gbsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

gbsv_dcomplex_parameters:: ~gbsv_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbsv_dcomplex_parameters object: destructor invoked. \n");
#endif
   gbsv_free();
}


//  Test fixture class definition
class zgbsv_test  : public  ::testing::Test {
public:
   gbsv_dcomplex_parameters  *zgbsv_obj;
   void SetUp();  
   void TearDown () { delete zgbsv_obj; }
};


void zgbsv_test::SetUp(){

    /* LAPACKE ZGBSV prototype */
    typedef int (*Fptr_NL_LAPACKE_zgbsv) ( int matrix_layout, 
                          lapack_int n, lapack_int kl, lapack_int ku,
                        lapack_int nrhs,  lapack_complex_double * ab,
                                 lapack_int ldab,  lapack_int * ipiv,
                          lapack_complex_double * b, lapack_int ldb );

    Fptr_NL_LAPACKE_zgbsv ZGBSV;

    zgbsv_obj = new  gbsv_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kl,
                                         lin_solver_paramslist[idx].ku,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zgbsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgbsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgbsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgbsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGBSV = (Fptr_NL_LAPACKE_zgbsv)dlsym(zgbsv_obj->hModule, "LAPACKE_zgbsv");
    ASSERT_TRUE(ZGBSV != NULL) << "failed to get the Netlib LAPACKE_zgbsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zgbsv_obj->inforef = ZGBSV( zgbsv_obj->matrix_layout,
                                   zgbsv_obj->n,
                                  zgbsv_obj->kl, zgbsv_obj->ku,
                                  zgbsv_obj->nrhs,
                       ( lapack_complex_double *)zgbsv_obj->abref,
                                  zgbsv_obj->ldab, zgbsv_obj->ipivref,
                                  zgbsv_obj->bref, zgbsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zgbsv_obj->info = LAPACKE_zgbsv( zgbsv_obj->matrix_layout,
                                    zgbsv_obj->n, zgbsv_obj->kl,
                                    zgbsv_obj->ku,zgbsv_obj->nrhs,
                                  ( lapack_complex_double *)zgbsv_obj->ab,
                                    zgbsv_obj->ldab, zgbsv_obj->ipiv,
                                    zgbsv_obj->b, zgbsv_obj->ldb );

    if( zgbsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgbsv is wrong\n",
                    zgbsv_obj->info );
    }
    if( zgbsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgbsv is wrong\n",
        zgbsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zgbsv_obj->diff =  computeDiff_z( zgbsv_obj->b_bufsize,
                           zgbsv_obj->b, zgbsv_obj->bref );
}

TEST_F(zgbsv_test, zgbsv1) {
    EXPECT_NEAR(0.0, zgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbsv_test, zgbsv2) {
    EXPECT_NEAR(0.0, zgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbsv_test, zgbsv3) {
    EXPECT_NEAR(0.0, zgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbsv_test, zgbsv4) {
    EXPECT_NEAR(0.0, zgbsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
