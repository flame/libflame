#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define geequ_free() \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (r != NULL)    free (r   ); \
       if (rref != NULL) free (rref); \
       if (c != NULL)    free (c   ); \
       if (cref != NULL) free (cref); \
       if( hModule != NULL) dlclose(hModule); \
       if( dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin geequ_float_parameters  class definition */
class geequ_float_parameters{
   public:
      int a_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int m; // rows in A
      lapack_int n; // Columns in A
      float *a, *aref; //Input matrix array
      lapack_int lda;  //  leading dimension of 'a'

      /* Output parameters */
      float * r, *rref; // row scale factors, array if size 'm'
      float * c, *cref; // colum scale factors, , array if size 'n'
      float rowcnd, rowcndref; // ratio of the smallest r[i] to the largest r[i]
      float colcnd, colcndref; // ratio of the smallest c[i] to the largest c[i]
      float amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      geequ_float_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i);
             
      ~geequ_float_parameters ();
};  /* end of geequ_float_parameters  class definition */


/* Constructor geequ_float_parameters definition */
geequ_float_parameters:: geequ_float_parameters ( int matrix_layout_i, 
      lapack_int m_i, lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    m = m_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

     lda = m; // as per API spec, lda≥ max(1, m).

#if LAPACKE_TEST_VERBOSE
   printf(" \n geequ float:  m: %d, n: %d lda: %d \n", m, n, lda);
#endif
    if(matrix_layout==LAPACK_COL_MAJOR){
        a_bufsize = lda*n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        a_bufsize = lda*m;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &r, &rref, m);
    lapacke_gtest_alloc_float_buffer_pair( &c, &cref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (r==NULL) || (rref==NULL) ||  \
        (c==NULL) || (cref==NULL) ){
       geequ_free();
       EXPECT_FALSE( true) << "geequ_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, a_bufsize);

   } /* end of Constructor  */

geequ_float_parameters:: ~geequ_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" geequ_float_parameters object: destructor invoked. \n");
#endif
   geequ_free();
}


//  Test fixture class definition
class sgeequ_test  : public  ::testing::Test {
public:
   geequ_float_parameters  *sgeequ_obj;
   void SetUp();  
   void TearDown () { delete sgeequ_obj; }
};


void sgeequ_test::SetUp(){

    /* LAPACKE SGEEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeequ) ( int matrix_layout, lapack_int m,
                               lapack_int n, const float *a, lapack_int lda,
                               float *r, float *c, float *rowcnd, 
                               float *colcnd, float *amax  );

    Fptr_NL_LAPACKE_sgeequ SGEEQU;

    sgeequ_obj = new  geequ_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    sgeequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGEEQU = (Fptr_NL_LAPACKE_sgeequ)dlsym(sgeequ_obj->hModule, "LAPACKE_sgeequ");
    ASSERT_TRUE(SGEEQU != NULL) << "failed to get the Netlib LAPACKE_sgeequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    sgeequ_obj->inforef = SGEEQU( sgeequ_obj->matrix_layout,
                                  sgeequ_obj->m,
                                  sgeequ_obj->n,
                                  sgeequ_obj->aref,
                                  sgeequ_obj->lda,
                                  sgeequ_obj->rref,
                                  sgeequ_obj->cref,
                                  &sgeequ_obj->rowcndref,
                                  &sgeequ_obj->colcndref,
                                  &sgeequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    sgeequ_obj->info = LAPACKE_sgeequ( sgeequ_obj->matrix_layout,
                                  sgeequ_obj->m,
                                  sgeequ_obj->n,
                                  sgeequ_obj->a,
                                  sgeequ_obj->lda,
                                  sgeequ_obj->r,
                                  sgeequ_obj->c,
                                  &sgeequ_obj->rowcnd,
                                  &sgeequ_obj->colcnd,
                                  &sgeequ_obj->amax);

    if( sgeequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgeequ \
        is wrong\n", sgeequ_obj->info );
    }
    if( sgeequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeequ is wrong\n",
        sgeequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    sgeequ_obj->diff =  computeDiff_s(  sgeequ_obj->m,
                                        sgeequ_obj->r,
                                        sgeequ_obj->rref );
    sgeequ_obj->diff +=  computeDiff_s(  sgeequ_obj->n,
                                        sgeequ_obj->c,
                                        sgeequ_obj->cref );
}

TEST_F(sgeequ_test, sgeequ1) {
    EXPECT_NEAR(0.0, sgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->rowcnd, sgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->colcnd, sgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->amax, sgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeequ_test, sgeequ2) {
    EXPECT_NEAR(0.0, sgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->rowcnd, sgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->colcnd, sgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->amax, sgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeequ_test, sgeequ3) {
    EXPECT_NEAR(0.0, sgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->rowcnd, sgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->colcnd, sgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->amax, sgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeequ_test, sgeequ4) {
    EXPECT_NEAR(0.0, sgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->rowcnd, sgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->colcnd, sgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequ_obj->amax, sgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin geequ_double_parameters  class definition */
class geequ_double_parameters{
   public:
      int a_bufsize;
      double diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int m; // rows in A
      lapack_int n; // Columns in A
      double *a, *aref; //Input matrix array
      lapack_int lda;  //  leading dimension of 'a'

      /* Output parameters */
      double * r, *rref; // row scale factors, array if size 'm'
      double * c, *cref; // colum scale factors, , array if size 'n'
      double rowcnd, rowcndref; // ratio of the smallest r[i] to the largest r[i]
      double colcnd, colcndref; // ratio of the smallest c[i] to the largest c[i]
      double amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      geequ_double_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i);
             
      ~geequ_double_parameters ();
};  /* end of geequ_double_parameters  class definition */


/* Constructor geequ_double_parameters definition */
geequ_double_parameters:: geequ_double_parameters ( int matrix_layout_i, 
      lapack_int m_i, lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    m = m_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

     lda = m; // as per API spec, lda≥ max(1, m).

#if LAPACKE_TEST_VERBOSE
   printf(" \n geequ double:  m: %d, n: %d lda: %d \n", m, n, lda);
#endif
    if(matrix_layout==LAPACK_COL_MAJOR){
        a_bufsize = lda*n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        a_bufsize = lda*m;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &r, &rref, m);
    lapacke_gtest_alloc_double_buffer_pair( &c, &cref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (r==NULL) || (rref==NULL) ||  \
        (c==NULL) || (cref==NULL) ){
       geequ_free();
       EXPECT_FALSE( true) << "geequ_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, a_bufsize);

   } /* end of Constructor  */

geequ_double_parameters:: ~geequ_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" geequ_double_parameters object: destructor invoked. \n");
#endif
   geequ_free();
}


//  Test fixture class definition
class dgeequ_test  : public  ::testing::Test {
public:
   geequ_double_parameters  *dgeequ_obj;
   void SetUp();  
   void TearDown () { delete dgeequ_obj; }
};


void dgeequ_test::SetUp(){

    /* LAPACKE DGEEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeequ) ( int matrix_layout, lapack_int m,
                               lapack_int n, const double *a, lapack_int lda,
                               double *r, double *c, double *rowcnd, 
                               double *colcnd, double *amax  );

    Fptr_NL_LAPACKE_dgeequ DGEEQU;

    dgeequ_obj = new  geequ_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    dgeequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGEEQU = (Fptr_NL_LAPACKE_dgeequ)dlsym(dgeequ_obj->hModule, "LAPACKE_dgeequ");
    ASSERT_TRUE(DGEEQU != NULL) << "failed to get the Netlib LAPACKE_dgeequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    dgeequ_obj->inforef = DGEEQU( dgeequ_obj->matrix_layout,
                                  dgeequ_obj->m,
                                  dgeequ_obj->n,
                                  dgeequ_obj->aref,
                                  dgeequ_obj->lda,
                                  dgeequ_obj->rref,
                                  dgeequ_obj->cref,
                                  &dgeequ_obj->rowcndref,
                                  &dgeequ_obj->colcndref,
                                  &dgeequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    dgeequ_obj->info = LAPACKE_dgeequ( dgeequ_obj->matrix_layout,
                                  dgeequ_obj->m,
                                  dgeequ_obj->n,
                                  dgeequ_obj->a,
                                  dgeequ_obj->lda,
                                  dgeequ_obj->r,
                                  dgeequ_obj->c,
                                  &dgeequ_obj->rowcnd,
                                  &dgeequ_obj->colcnd,
                                  &dgeequ_obj->amax);

    if( dgeequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgeequ \
        is wrong\n", dgeequ_obj->info );
    }
    if( dgeequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeequ is wrong\n",
        dgeequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dgeequ_obj->diff =  computeDiff_d(  dgeequ_obj->m,
                                        dgeequ_obj->r,
                                        dgeequ_obj->rref );
    dgeequ_obj->diff +=  computeDiff_d(  dgeequ_obj->n,
                                        dgeequ_obj->c,
                                        dgeequ_obj->cref );
}

TEST_F(dgeequ_test, dgeequ1) {
    EXPECT_NEAR(0.0, dgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->rowcnd, dgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->colcnd, dgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->amax, dgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeequ_test, dgeequ2) {
    EXPECT_NEAR(0.0, dgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->rowcnd, dgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->colcnd, dgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->amax, dgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeequ_test, dgeequ3) {
    EXPECT_NEAR(0.0, dgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->rowcnd, dgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->colcnd, dgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->amax, dgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeequ_test, dgeequ4) {
    EXPECT_NEAR(0.0, dgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->rowcnd, dgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->colcnd, dgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequ_obj->amax, dgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin geequ_scomplex_parameters  class definition */
class geequ_scomplex_parameters{
   public:
      int a_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int m; // rows in A
      lapack_int n; // Columns in A
      lapack_complex_float *a, *aref; //Input matrix array
      lapack_int lda;  //  leading dimension of 'a'

      /* Output parameters */
      float * r, *rref; // row scale factors, array if size 'm'
      float * c, *cref; // colum scale factors, , array if size 'n'
      float rowcnd, rowcndref; // ratio of the smallest r[i] to the largest r[i]
      float colcnd, colcndref; // ratio of the smallest c[i] to the largest c[i]
      float amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      geequ_scomplex_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i);
             
      ~geequ_scomplex_parameters ();
};  /* end of geequ_scomplex_parameters  class definition */


/* Constructor geequ_scomplex_parameters definition */
geequ_scomplex_parameters:: geequ_scomplex_parameters ( int matrix_layout_i, 
      lapack_int m_i, lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    m = m_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

     lda = m; // as per API spec, lda≥ max(1, m).

#if LAPACKE_TEST_VERBOSE
   printf(" \n geequ scomplex:  m: %d, n: %d lda: %d \n", m, n, lda);
#endif
    if(matrix_layout==LAPACK_COL_MAJOR){
        a_bufsize = lda*n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        a_bufsize = lda*m;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &r, &rref, m);
    lapacke_gtest_alloc_float_buffer_pair( &c, &cref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (r==NULL) || (rref==NULL) ||  \
        (c==NULL) || (cref==NULL) ){
       geequ_free();
       EXPECT_FALSE( true) << "geequ_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, a_bufsize);

   } /* end of Constructor  */

geequ_scomplex_parameters:: ~geequ_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" geequ_scomplex_parameters object: destructor invoked. \n");
#endif
   geequ_free();
}


//  Test fixture class definition
class cgeequ_test  : public  ::testing::Test {
public:
   geequ_scomplex_parameters  *cgeequ_obj;
   void SetUp();  
   void TearDown () { delete cgeequ_obj; }
};


void cgeequ_test::SetUp(){

    /* LAPACKE CGEEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeequ) ( int matrix_layout, lapack_int m,
                               lapack_int n, const lapack_complex_float *a, 
                               lapack_int lda, float *r, float *c, 
                               float *rowcnd, float *colcnd, float *amax  );

    Fptr_NL_LAPACKE_cgeequ CGEEQU;

    cgeequ_obj = new  geequ_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    cgeequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgeequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgeequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgeequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGEEQU = (Fptr_NL_LAPACKE_cgeequ)dlsym(cgeequ_obj->hModule, "LAPACKE_cgeequ");
    ASSERT_TRUE(CGEEQU != NULL) << "failed to get the Netlib LAPACKE_cgeequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    cgeequ_obj->inforef = CGEEQU( cgeequ_obj->matrix_layout,
                                  cgeequ_obj->m,
                                  cgeequ_obj->n,
                                  cgeequ_obj->aref,
                                  cgeequ_obj->lda,
                                  cgeequ_obj->rref,
                                  cgeequ_obj->cref,
                                  &cgeequ_obj->rowcndref,
                                  &cgeequ_obj->colcndref,
                                  &cgeequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    cgeequ_obj->info = LAPACKE_cgeequ( cgeequ_obj->matrix_layout,
                                  cgeequ_obj->m,
                                  cgeequ_obj->n,
                                  cgeequ_obj->a,
                                  cgeequ_obj->lda,
                                  cgeequ_obj->r,
                                  cgeequ_obj->c,
                                  &cgeequ_obj->rowcnd,
                                  &cgeequ_obj->colcnd,
                                  &cgeequ_obj->amax);

    if( cgeequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgeequ \
        is wrong\n", cgeequ_obj->info );
    }
    if( cgeequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeequ is wrong\n",
        cgeequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cgeequ_obj->diff =  computeDiff_s(  cgeequ_obj->m,
                                        cgeequ_obj->r,
                                        cgeequ_obj->rref );
    cgeequ_obj->diff +=  computeDiff_s(  cgeequ_obj->n,
                                        cgeequ_obj->c,
                                        cgeequ_obj->cref );
}

TEST_F(cgeequ_test, cgeequ1) {
    EXPECT_NEAR(0.0, cgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->rowcnd, cgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->colcnd, cgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->amax, cgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeequ_test, cgeequ2) {
    EXPECT_NEAR(0.0, cgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->rowcnd, cgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->colcnd, cgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->amax, cgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeequ_test, cgeequ3) {
    EXPECT_NEAR(0.0, cgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->rowcnd, cgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->colcnd, cgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->amax, cgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeequ_test, cgeequ4) {
    EXPECT_NEAR(0.0, cgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->rowcnd, cgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->colcnd, cgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequ_obj->amax, cgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}


/* Begin geequ_dcomplex_parameters  class definition */
class geequ_dcomplex_parameters{
   public:
      int a_bufsize;
      double diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int m; // rows in A
      lapack_int n; // Columns in A
      lapack_complex_double *a, *aref; //Input matrix array
      lapack_int lda;  //  leading dimension of 'a'

      /* Output parameters */
      double * r, *rref; // row scale factors, array if size 'm'
      double * c, *cref; // colum scale factors, , array if size 'n'
      double rowcnd, rowcndref; // ratio of the smallest r[i] to the largest r[i]
      double colcnd, colcndref; // ratio of the smallest c[i] to the largest c[i]
      double amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      geequ_dcomplex_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i);
             
      ~geequ_dcomplex_parameters ();
};  /* end of geequ_dcomplex_parameters  class definition */


/* Constructor geequ_dcomplex_parameters definition */
geequ_dcomplex_parameters:: geequ_dcomplex_parameters ( int matrix_layout_i, 
      lapack_int m_i, lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    m = m_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

     lda = m; // as per API spec, lda≥ max(1, m).

#if LAPACKE_TEST_VERBOSE
   printf(" \n geequ dcomplex:  m: %d, n: %d lda: %d \n", m, n, lda);
#endif
    if(matrix_layout==LAPACK_COL_MAJOR){
        a_bufsize = lda*n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        a_bufsize = lda*m;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &r, &rref, m);
    lapacke_gtest_alloc_double_buffer_pair( &c, &cref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (r==NULL) || (rref==NULL) ||  \
        (c==NULL) || (cref==NULL) ){
       geequ_free();
       EXPECT_FALSE( true) << "geequ_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, a_bufsize);

   } /* end of Constructor  */

geequ_dcomplex_parameters:: ~geequ_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" geequ_dcomplex_parameters object: destructor invoked. \n");
#endif
   geequ_free();
}


//  Test fixture class definition
class zgeequ_test  : public  ::testing::Test {
public:
   geequ_dcomplex_parameters  *zgeequ_obj;
   void SetUp();  
   void TearDown () { delete zgeequ_obj; }
};


void zgeequ_test::SetUp(){

    /* LAPACKE ZGEEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeequ) ( int matrix_layout, lapack_int m,
                               lapack_int n, const lapack_complex_double *a, 
                               lapack_int lda, double *r, double *c, 
                               double *rowcnd, double *colcnd, double *amax  );

    Fptr_NL_LAPACKE_zgeequ ZGEEQU;

    zgeequ_obj = new  geequ_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    zgeequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgeequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgeequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgeequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGEEQU = (Fptr_NL_LAPACKE_zgeequ)dlsym(zgeequ_obj->hModule, "LAPACKE_zgeequ");
    ASSERT_TRUE(ZGEEQU != NULL) << "failed to get the Netlib LAPACKE_zgeequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zgeequ_obj->inforef = ZGEEQU( zgeequ_obj->matrix_layout,
                                  zgeequ_obj->m,
                                  zgeequ_obj->n,
                                  zgeequ_obj->aref,
                                  zgeequ_obj->lda,
                                  zgeequ_obj->rref,
                                  zgeequ_obj->cref,
                                  &zgeequ_obj->rowcndref,
                                  &zgeequ_obj->colcndref,
                                  &zgeequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    zgeequ_obj->info = LAPACKE_zgeequ( zgeequ_obj->matrix_layout,
                                  zgeequ_obj->m,
                                  zgeequ_obj->n,
                                  zgeequ_obj->a,
                                  zgeequ_obj->lda,
                                  zgeequ_obj->r,
                                  zgeequ_obj->c,
                                  &zgeequ_obj->rowcnd,
                                  &zgeequ_obj->colcnd,
                                  &zgeequ_obj->amax);

    if( zgeequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgeequ \
        is wrong\n", zgeequ_obj->info );
    }
    if( zgeequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeequ is wrong\n",
        zgeequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zgeequ_obj->diff =  computeDiff_d(  zgeequ_obj->m,
                                        zgeequ_obj->r,
                                        zgeequ_obj->rref );
    zgeequ_obj->diff +=  computeDiff_d(  zgeequ_obj->n,
                                        zgeequ_obj->c,
                                        zgeequ_obj->cref );
}

TEST_F(zgeequ_test, zgeequ1) {
    EXPECT_NEAR(0.0, zgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->rowcnd, zgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->colcnd, zgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->amax, zgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeequ_test, zgeequ2) {
    EXPECT_NEAR(0.0, zgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->rowcnd, zgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->colcnd, zgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->amax, zgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeequ_test, zgeequ3) {
    EXPECT_NEAR(0.0, zgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->rowcnd, zgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->colcnd, zgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->amax, zgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeequ_test, zgeequ4) {
    EXPECT_NEAR(0.0, zgeequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->rowcnd, zgeequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->colcnd, zgeequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequ_obj->amax, zgeequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}
