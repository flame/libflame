#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define gtsv_free() \
       if (b != NULL)    free (b   ); \
       if (bref != NULL) free (bref); \
       if (d != NULL)    free (d   ); \
       if (dref != NULL) free (dref); \
       if (du != NULL)    free (du   ); \
       if (duref != NULL) free (duref); \
       if (dl != NULL)    free (dl   ); \
       if (dlref != NULL) free (dlref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin gtsv_double_parameters  class definition */
class gtsv_double_parameters{
   public:
      int b_bufsize;
      double diff; // to capture the netlib, libflame o/ps
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double * dl, *dlref; // (n - 1) multipliers that define the matrix L
      double * d, *dref; //   n diagonal elements of the upper triangular matrix U
      double * du, *duref; //  (n - 1) elements of the first superdiagonal of U.
      void *hModule, *dModule;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gtsv_double_parameters ( int matrix_layout_i, 
              lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gtsv_double_parameters ();
};  /* end of gtsv_double_parameters  class definition */


/* Constructor gtsv_double_parameters definition */
gtsv_double_parameters:: gtsv_double_parameters ( int matrix_layout_i, 
      lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gtsv Double:  n: %d, ldb: %d nrhs: %d \n",
             n, ldb, nrhs);
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
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_double_buffer_pair( &dl, &dlref, (n-1));
    lapacke_gtest_alloc_double_buffer_pair( &du, &duref, (n-1));

    if( (b==NULL) || (bref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (d==NULL) || (dref==NULL) ||  \
        (du==NULL) || (duref==NULL)  ){
       gtsv_free();
       EXPECT_FALSE( true) << "gtsv_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_double_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_double_buffer_pair_rand( dl, dlref, (n-1));
    lapacke_gtest_init_double_buffer_pair_rand( du, duref, (n-1));
   

   } /* end of Constructor  */

gtsv_double_parameters:: ~gtsv_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtsv_double_parameters object: destructor invoked. \n");
#endif
   gtsv_free();
}


//  Test fixture class definition
class dgtsv_test  : public  ::testing::Test {
public:
   gtsv_double_parameters  *dgtsv_obj;
   void SetUp();  
   void TearDown () { delete dgtsv_obj; }
};


void dgtsv_test::SetUp(){

    /* LAPACKE DGTSV prototype */
    typedef int (*Fptr_NL_LAPACKE_dgtsv) ( int matrix_layout,  
    lapack_int n, lapack_int nrhs, const double * dl, const double * d, 
    const double * du, double * b, lapack_int ldb );

    Fptr_NL_LAPACKE_dgtsv DGTSV;

    dgtsv_obj = new  gtsv_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dgtsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgtsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgtsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgtsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGTSV = (Fptr_NL_LAPACKE_dgtsv)dlsym(dgtsv_obj->hModule, "LAPACKE_dgtsv");
    ASSERT_TRUE(DGTSV != NULL) << "failed to get the Netlib LAPACKE_dgtsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    dgtsv_obj->inforef = DGTSV( dgtsv_obj->matrix_layout,
                                   dgtsv_obj->n,
                                  dgtsv_obj->nrhs,
                                  dgtsv_obj->dlref,
                                  dgtsv_obj->dref,
                                  dgtsv_obj->duref,
                                  dgtsv_obj->bref, dgtsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */

    dgtsv_obj->info = LAPACKE_dgtsv( dgtsv_obj->matrix_layout,
                 dgtsv_obj->n, dgtsv_obj->nrhs,
                                  dgtsv_obj->dl,
                                  dgtsv_obj->d,
                                  dgtsv_obj->du,
                                 dgtsv_obj->b, dgtsv_obj->ldb );

    if( dgtsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
                   LAPACKE_dgtsv is wrong\n", dgtsv_obj->info );
    }
    if( dgtsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgtsv \
                            is wrong\n",  dgtsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dgtsv_obj->diff =  computeDiff_d( dgtsv_obj->b_bufsize,
                           dgtsv_obj->b, dgtsv_obj->bref );
}

TEST_F(dgtsv_test, dgtsv1) {
    EXPECT_NEAR(0.0, dgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgtsv_test, dgtsv2) {
    EXPECT_NEAR(0.0, dgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgtsv_test, dgtsv3) {
    EXPECT_NEAR(0.0, dgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgtsv_test, dgtsv4) {
    EXPECT_NEAR(0.0, dgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}



/* Begin gtsv_float_parameters  class definition */
class gtsv_float_parameters{
   public:
      int b_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float * dl, *dlref; // (n - 1) multipliers that define the matrix L
      float * d, *dref; //   n diagonal elements of the upper triangular matrix U
      float * du, *duref; //  (n - 1) elements of the first superdiagonal of U.
      void *hModule, *dModule;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gtsv_float_parameters ( int matrix_layout_i, 
              lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gtsv_float_parameters ();
};  /* end of gtsv_float_parameters  class definition */


/* Constructor gtsv_float_parameters definition */
gtsv_float_parameters:: gtsv_float_parameters ( int matrix_layout_i, 
      lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gtsv Double:  n: %d, ldb: %d nrhs: %d \n",
             n, ldb, nrhs);
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
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_float_buffer_pair( &dl, &dlref, (n-1));
    lapacke_gtest_alloc_float_buffer_pair( &du, &duref, (n-1));

    if( (b==NULL) || (bref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (d==NULL) || (dref==NULL) ||  \
        (du==NULL) || (duref==NULL) ){
       gtsv_free();
       EXPECT_FALSE( true) << "gtsv_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_float_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_float_buffer_pair_rand( dl, dlref, (n-1));
    lapacke_gtest_init_float_buffer_pair_rand( du, duref, (n-1));
   

   } /* end of Constructor  */

gtsv_float_parameters:: ~gtsv_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtsv_float_parameters object: destructor invoked. \n");
#endif
   gtsv_free();
}


//  Test fixture class definition
class sgtsv_test  : public  ::testing::Test {
public:
   gtsv_float_parameters  *sgtsv_obj;
   void SetUp();  
   void TearDown () { delete sgtsv_obj; }
};


void sgtsv_test::SetUp(){

    /* LAPACKE SGTSV prototype */
    typedef int (*Fptr_NL_LAPACKE_sgtsv) ( int matrix_layout,  
    lapack_int n, lapack_int nrhs, const float * dl, const float * d, 
    const float * du, float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_sgtsv SGTSV;

    sgtsv_obj = new  gtsv_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    sgtsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgtsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgtsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgtsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGTSV = (Fptr_NL_LAPACKE_sgtsv)dlsym(sgtsv_obj->hModule, "LAPACKE_sgtsv");
    ASSERT_TRUE(SGTSV != NULL) << "failed to get the Netlib LAPACKE_sgtsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    sgtsv_obj->inforef = SGTSV( sgtsv_obj->matrix_layout,
                                   sgtsv_obj->n,
                                  sgtsv_obj->nrhs,
                                  sgtsv_obj->dlref,
                                  sgtsv_obj->dref,
                                  sgtsv_obj->duref,
                                  sgtsv_obj->bref, sgtsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    sgtsv_obj->info = LAPACKE_sgtsv( sgtsv_obj->matrix_layout,
                 sgtsv_obj->n, sgtsv_obj->nrhs,
                                  sgtsv_obj->dl,
                                  sgtsv_obj->d,
                                  sgtsv_obj->du,
                                 sgtsv_obj->b, sgtsv_obj->ldb );


    if( sgtsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgtsv \
        is wrong\n", sgtsv_obj->info );
    }
    if( sgtsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgtsv is wrong\n",
        sgtsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    sgtsv_obj->diff =  computeDiff_s( sgtsv_obj->b_bufsize,
                           sgtsv_obj->b, sgtsv_obj->bref );
}

TEST_F(sgtsv_test, sgtsv1) {
    EXPECT_NEAR(0.0, sgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgtsv_test, sgtsv2) {
    EXPECT_NEAR(0.0, sgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgtsv_test, sgtsv3) {
    EXPECT_NEAR(0.0, sgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgtsv_test, sgtsv4) {
    EXPECT_NEAR(0.0, sgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gtsv_scomplex_parameters  class definition */
class gtsv_scomplex_parameters{
   public:
      void *hModule, *dModule;
      int b_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *dl, *dlref; // (n - 1) multipliers that define the matrix L
      lapack_complex_float *d, *dref; //   n diagonal elements of the upper triangular matrix U
      lapack_complex_float *du, *duref; //  (n - 1) elements of the first superdiagonal of U.

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gtsv_scomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gtsv_scomplex_parameters ();
};  /* end of gtsv_scomplex_parameters  class definition */


/* Constructor gtsv_scomplex_parameters definition */
gtsv_scomplex_parameters:: gtsv_scomplex_parameters ( int matrix_layout_i, 
      lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gtsv Double:  n: %d, ldb: %d nrhs: %d \n",
             n, ldb, nrhs);
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
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &dl, &dlref, (n-1));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &du, &duref, (n-1));

    if( (b==NULL) || (bref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (d==NULL) || (dref==NULL) ||  \
        (du==NULL) || (duref==NULL) ){
       gtsv_free();
       EXPECT_FALSE( true) << "gtsv_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_scomplex_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( dl, dlref, (n-1));
    lapacke_gtest_init_scomplex_buffer_pair_rand( du, duref, (n-1));

   } /* end of Constructor  */

gtsv_scomplex_parameters:: ~gtsv_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtsv_scomplex_parameters object: destructor invoked. \n");
#endif
   gtsv_free();
}


//  Test fixture class definition
class cgtsv_test  : public  ::testing::Test {
public:
   gtsv_scomplex_parameters  *cgtsv_obj;
   void SetUp();  
   void TearDown () { delete cgtsv_obj; }
};


void cgtsv_test::SetUp(){

    /* LAPACKE CGTSV prototype */
    typedef int (*Fptr_NL_LAPACKE_cgtsv) ( int matrix_layout,  
    lapack_int n, lapack_int nrhs, const lapack_complex_float * dl, 
    const lapack_complex_float * d, 
    const lapack_complex_float * du, 
    lapack_complex_float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_cgtsv CGTSV;

     /* LAPACKE CGTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cgttrf) ( lapack_int n,
    lapack_complex_float * dl, lapack_complex_float * d,
    lapack_complex_float * du, lapack_complex_float * du2,
    lapack_int * ipiv);

    Fptr_NL_LAPACKE_cgttrf CGTTRF;

    cgtsv_obj = new  gtsv_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    cgtsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgtsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgtsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgtsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGTSV = (Fptr_NL_LAPACKE_cgtsv)dlsym(cgtsv_obj->hModule, "LAPACKE_cgtsv");
    ASSERT_TRUE(CGTSV != NULL) << "failed to get the Netlib LAPACKE_cgtsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    cgtsv_obj->inforef = CGTSV( cgtsv_obj->matrix_layout,
                                   cgtsv_obj->n,
                                  cgtsv_obj->nrhs,
                                  cgtsv_obj->dlref,
                                  cgtsv_obj->dref,
                                  cgtsv_obj->duref,
                                  cgtsv_obj->bref, cgtsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    cgtsv_obj->info = LAPACKE_cgtsv( cgtsv_obj->matrix_layout,
                 cgtsv_obj->n, cgtsv_obj->nrhs,
                                  cgtsv_obj->dl,
                                  cgtsv_obj->d,
                                  cgtsv_obj->du,
                                 cgtsv_obj->b, cgtsv_obj->ldb );

    if( cgtsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgtsv is wrong\n",
                    cgtsv_obj->info );
    }
    if( cgtsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgtsv is wrong\n",
        cgtsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cgtsv_obj->diff =  computeDiff_c( cgtsv_obj->b_bufsize,
                           cgtsv_obj->b, cgtsv_obj->bref );
}

TEST_F(cgtsv_test, cgtsv1) {
    EXPECT_NEAR(0.0, cgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgtsv_test, cgtsv2) {
    EXPECT_NEAR(0.0, cgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgtsv_test, cgtsv3) {
    EXPECT_NEAR(0.0, cgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgtsv_test, cgtsv4) {
    EXPECT_NEAR(0.0, cgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gtsv_dcomplex_parameters  class definition */
class gtsv_dcomplex_parameters{
   public:
      void *hModule, *dModule;
      int b_bufsize;
      double diff; // to capture the netlib, libflame o/ps
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *dl, *dlref; // (n - 1) multipliers that define the matrix L
      lapack_complex_double *d, *dref; //   n diagonal elements of the upper triangular matrix U
      lapack_complex_double *du, *duref; //  (n - 1) elements of the first superdiagonal of U.

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      gtsv_dcomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~gtsv_dcomplex_parameters ();
};  /* end of gtsv_dcomplex_parameters  class definition */


/* Constructor gtsv_dcomplex_parameters definition */
gtsv_dcomplex_parameters:: gtsv_dcomplex_parameters ( int matrix_layout_i, 
      lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n gtsv Double:  n: %d, ldb: %d nrhs: %d \n",
             n, ldb, nrhs);
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
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &dl, &dlref, (n-1));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &du, &duref, (n-1));

    if( (b==NULL) || (bref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (d==NULL) || (dref==NULL) ||  \
        (du==NULL) || (duref==NULL) ){
       gtsv_free();
       EXPECT_FALSE( true) << "gtsv_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( dl, dlref, (n-1));
    lapacke_gtest_init_dcomplex_buffer_pair_rand( du, duref, (n-1));
   

   } /* end of Constructor  */

gtsv_dcomplex_parameters:: ~gtsv_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtsv_dcomplex_parameters object: destructor invoked. \n");
#endif
   gtsv_free();
}


//  Test fixture class definition
class zgtsv_test  : public  ::testing::Test {
public:
   gtsv_dcomplex_parameters  *zgtsv_obj;
   void SetUp();  
   void TearDown () { delete zgtsv_obj; }
};


void zgtsv_test::SetUp(){

    /* LAPACKE ZGTSV prototype */
    typedef int (*Fptr_NL_LAPACKE_zgtsv) ( int matrix_layout,  
    lapack_int n, lapack_int nrhs, const lapack_complex_double * dl,
    const lapack_complex_double * d, const lapack_complex_double * du,
    lapack_complex_double * b, lapack_int ldb );

    Fptr_NL_LAPACKE_zgtsv ZGTSV;

    zgtsv_obj = new  gtsv_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    zgtsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgtsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgtsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgtsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGTSV = (Fptr_NL_LAPACKE_zgtsv)dlsym(zgtsv_obj->hModule, "LAPACKE_zgtsv");
    ASSERT_TRUE(ZGTSV != NULL) << "failed to get the Netlib LAPACKE_zgtsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zgtsv_obj->inforef = ZGTSV( zgtsv_obj->matrix_layout,
                                   zgtsv_obj->n,
                                  zgtsv_obj->nrhs,
                                  zgtsv_obj->dlref,
                                  zgtsv_obj->dref,
                                  zgtsv_obj->duref,
                                  zgtsv_obj->bref, zgtsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zgtsv_obj->info = LAPACKE_zgtsv( zgtsv_obj->matrix_layout,
                 zgtsv_obj->n, zgtsv_obj->nrhs,
                                  zgtsv_obj->dl,
                                  zgtsv_obj->d,
                                  zgtsv_obj->du,
                                 zgtsv_obj->b, zgtsv_obj->ldb );

    if( zgtsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgtsv \
        is wrong\n", zgtsv_obj->info );
    }
    if( zgtsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgtsv is wrong\n",
        zgtsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zgtsv_obj->diff =  computeDiff_z( zgtsv_obj->b_bufsize,
                           zgtsv_obj->b, zgtsv_obj->bref );
}

TEST_F(zgtsv_test, zgtsv1) {
    EXPECT_NEAR(0.0, zgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgtsv_test, zgtsv2) {
    EXPECT_NEAR(0.0, zgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgtsv_test, zgtsv3) {
    EXPECT_NEAR(0.0, zgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgtsv_test, zgtsv4) {
    EXPECT_NEAR(0.0, zgtsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
