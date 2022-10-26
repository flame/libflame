#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define gbequ_free() \
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

/* Begin gbequ_float_parameters  class definition */
class gbequ_float_parameters{
   public:
      int a_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int m; // rows in A
      lapack_int n; // Columns in A
      lapack_int kl; // number of subdiagonals
      lapack_int ku; // number of superdiagonals
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
      gbequ_float_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i,
              lapack_int kl_i, lapack_int ku_i);
             
      ~gbequ_float_parameters ();
};  /* end of gbequ_float_parameters  class definition */


/* Constructor gbequ_float_parameters definition */
gbequ_float_parameters:: gbequ_float_parameters ( int matrix_layout_i, 
      lapack_int m_i, lapack_int n_i, lapack_int kl_i, lapack_int ku_i) {

    matrix_layout = matrix_layout_i;
    m = m_i;
    n = n_i;
    kl = kl_i;
    ku = ku_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = kl+ku+2; // as per API spec, ldab≥kl+ku+1.

#if LAPACKE_TEST_VERBOSE
   printf(" \n gbequ float:  m: %d, n: %d lda: %d \n", m, n, lda);
#endif
    if(matrix_layout==LAPACK_COL_MAJOR){
        a_bufsize = lda*n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        lda = n;
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
       gbequ_free();
       EXPECT_FALSE( true) << "gbequ_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, a_bufsize);
    lapacke_gtest_init_float_buffer_pair_with_constant(r, rref, m, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(c, cref, n, 0.0);

   } /* end of Constructor  */

gbequ_float_parameters:: ~gbequ_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbequ_float_parameters object: destructor invoked. \n");
#endif
   gbequ_free();
}


//  Test fixture class definition
class sgbequ_test  : public  ::testing::Test {
public:
   gbequ_float_parameters  *sgbequ_obj;
   void SetUp();  
   void TearDown () { delete sgbequ_obj; }
};


void sgbequ_test::SetUp(){

    /* LAPACKE SGBEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_sgbequ) ( int matrix_layout, lapack_int m,
                               lapack_int n, lapack_int kl, lapack_int ku,
                               const float *a, lapack_int lda,
                               float *r, float *c, float *rowcnd, 
                               float *colcnd, float *amax  );

    Fptr_NL_LAPACKE_sgbequ SGBEQU;

    sgbequ_obj = new  gbequ_float_parameters(
                                lin_solver_paramslist[idx].matrix_layout,
                                            lin_solver_paramslist[idx].m,
                                            lin_solver_paramslist[idx].n,
                                            lin_solver_paramslist[idx].kl,
                                            lin_solver_paramslist[idx].ku
                                          );
    idx = Circular_Increment_Index(idx);

    sgbequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgbequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgbequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgbequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGBEQU = (Fptr_NL_LAPACKE_sgbequ)dlsym(sgbequ_obj->hModule, "LAPACKE_sgbequ");
    ASSERT_TRUE(SGBEQU != NULL) << "failed to get the Netlib LAPACKE_sgbequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    sgbequ_obj->inforef = SGBEQU( sgbequ_obj->matrix_layout,
                                  sgbequ_obj->m,
                                  sgbequ_obj->n,
                                  sgbequ_obj->kl,
                                  sgbequ_obj->ku,
                                  sgbequ_obj->aref,
                                  sgbequ_obj->lda,
                                  sgbequ_obj->rref,
                                  sgbequ_obj->cref,
                                  &sgbequ_obj->rowcndref,
                                  &sgbequ_obj->colcndref,
                                  &sgbequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    sgbequ_obj->info = LAPACKE_sgbequ( sgbequ_obj->matrix_layout,
                                  sgbequ_obj->m,
                                  sgbequ_obj->n,
                                  sgbequ_obj->kl,
                                  sgbequ_obj->ku,
                                  sgbequ_obj->a,
                                  sgbequ_obj->lda,
                                  sgbequ_obj->r,
                                  sgbequ_obj->c,
                                  &sgbequ_obj->rowcnd,
                                  &sgbequ_obj->colcnd,
                                  &sgbequ_obj->amax);

    if( sgbequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgbequ \
        is wrong\n", sgbequ_obj->info );
    }
    if( sgbequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgbequ is wrong\n",
        sgbequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    sgbequ_obj->diff =  computeDiff_s(  sgbequ_obj->m,
                                        sgbequ_obj->r,
                                        sgbequ_obj->rref );
    sgbequ_obj->diff +=  computeDiff_s(  sgbequ_obj->n,
                                        sgbequ_obj->c,
                                        sgbequ_obj->cref );
}

TEST_F(sgbequ_test, sgbequ1) {
    EXPECT_NEAR(0.0, sgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->rowcnd, sgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->colcnd, sgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->amax, sgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbequ_test, sgbequ2) {
    EXPECT_NEAR(0.0, sgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->rowcnd, sgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->colcnd, sgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->amax, sgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbequ_test, sgbequ3) {
    EXPECT_NEAR(0.0, sgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->rowcnd, sgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->colcnd, sgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->amax, sgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbequ_test, sgbequ4) {
    EXPECT_NEAR(0.0, sgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->rowcnd, sgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->colcnd, sgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequ_obj->amax, sgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gbequ_double_parameters  class definition */
class gbequ_double_parameters{
   public:
      int a_bufsize;
      double diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int m; // rows in A
      lapack_int n; // Columns in A
      lapack_int kl; // number of subdiagonals
      lapack_int ku; // number of superdiagonals
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
      gbequ_double_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i,
              lapack_int kl_i, lapack_int ku_i);
             
      ~gbequ_double_parameters ();
};  /* end of gbequ_double_parameters  class definition */


/* Constructor gbequ_double_parameters definition */
gbequ_double_parameters:: gbequ_double_parameters ( int matrix_layout_i, 
      lapack_int m_i, lapack_int n_i, lapack_int kl_i, lapack_int ku_i) {

    matrix_layout = matrix_layout_i;
    m = m_i;
    n = n_i;
    kl = kl_i;
    ku = ku_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = kl+ku+2; // as per API spec, ldab≥kl+ku+1.

#if LAPACKE_TEST_VERBOSE
   printf(" \n gbequ double:  m: %d, n: %d lda: %d \n", m, n, lda);
#endif
    if(matrix_layout==LAPACK_COL_MAJOR){
        a_bufsize = lda*n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        lda = n;
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
       gbequ_free();
       EXPECT_FALSE( true) << "gbequ_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, a_bufsize);
    lapacke_gtest_init_double_buffer_pair_with_constant(r, rref, m, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(c, cref, n, 0.0);

   } /* end of Constructor  */

gbequ_double_parameters:: ~gbequ_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbequ_double_parameters object: destructor invoked. \n");
#endif
   gbequ_free();
}


//  Test fixture class definition
class dgbequ_test  : public  ::testing::Test {
public:
   gbequ_double_parameters  *dgbequ_obj;
   void SetUp();  
   void TearDown () { delete dgbequ_obj; }
};


void dgbequ_test::SetUp(){

    /* LAPACKE DGBEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_dgbequ) ( int matrix_layout, lapack_int m,
                               lapack_int n, lapack_int kl, lapack_int ku,
                   const double *a, lapack_int lda,
                               double *r, double *c, double *rowcnd, 
                               double *colcnd, double *amax  );

    Fptr_NL_LAPACKE_dgbequ DGBEQU;

    dgbequ_obj = new  gbequ_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                            lin_solver_paramslist[idx].n,
                                            lin_solver_paramslist[idx].kl,
                                            lin_solver_paramslist[idx].ku
                                          );
    idx = Circular_Increment_Index(idx);

    dgbequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgbequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgbequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgbequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGBEQU = (Fptr_NL_LAPACKE_dgbequ)dlsym(dgbequ_obj->hModule, "LAPACKE_dgbequ");
    ASSERT_TRUE(DGBEQU != NULL) << "failed to get the Netlib LAPACKE_dgbequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    dgbequ_obj->inforef = DGBEQU( dgbequ_obj->matrix_layout,
                                  dgbequ_obj->m,
                                  dgbequ_obj->n,
                                  dgbequ_obj->kl,
                                  dgbequ_obj->ku,
                                  dgbequ_obj->aref,
                                  dgbequ_obj->lda,
                                  dgbequ_obj->rref,
                                  dgbequ_obj->cref,
                                  &dgbequ_obj->rowcndref,
                                  &dgbequ_obj->colcndref,
                                  &dgbequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    dgbequ_obj->info = LAPACKE_dgbequ( dgbequ_obj->matrix_layout,
                                  dgbequ_obj->m,
                                  dgbequ_obj->n,
                                  dgbequ_obj->kl,
                                  dgbequ_obj->ku,
                                  dgbequ_obj->a,
                                  dgbequ_obj->lda,
                                  dgbequ_obj->r,
                                  dgbequ_obj->c,
                                  &dgbequ_obj->rowcnd,
                                  &dgbequ_obj->colcnd,
                                  &dgbequ_obj->amax);

    if( dgbequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgbequ \
        is wrong\n", dgbequ_obj->info );
    }
    if( dgbequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgbequ is wrong\n",
        dgbequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dgbequ_obj->diff =  computeDiff_d(  dgbequ_obj->m,
                                        dgbequ_obj->r,
                                        dgbequ_obj->rref );
    dgbequ_obj->diff +=  computeDiff_d(  dgbequ_obj->n,
                                        dgbequ_obj->c,
                                        dgbequ_obj->cref );
}

TEST_F(dgbequ_test, dgbequ1) {
    EXPECT_NEAR(0.0, dgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->rowcnd, dgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->colcnd, dgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->amax, dgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbequ_test, dgbequ2) {
    EXPECT_NEAR(0.0, dgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->rowcnd, dgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->colcnd, dgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->amax, dgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbequ_test, dgbequ3) {
    EXPECT_NEAR(0.0, dgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->rowcnd, dgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->colcnd, dgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->amax, dgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbequ_test, dgbequ4) {
    EXPECT_NEAR(0.0, dgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->rowcnd, dgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->colcnd, dgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequ_obj->amax, dgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gbequ_scomplex_parameters  class definition */
class gbequ_scomplex_parameters{
   public:
      int a_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int m; // rows in A
      lapack_int n; // Columns in A
      lapack_int kl; // number of subdiagonals
      lapack_int ku; // number of superdiagonals
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
      gbequ_scomplex_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i,
              lapack_int kl_i, lapack_int ku_i);
             
      ~gbequ_scomplex_parameters ();
};  /* end of gbequ_scomplex_parameters  class definition */


/* Constructor gbequ_scomplex_parameters definition */
gbequ_scomplex_parameters:: gbequ_scomplex_parameters ( int matrix_layout_i, 
      lapack_int m_i, lapack_int n_i, lapack_int kl_i, lapack_int ku_i) {

    matrix_layout = matrix_layout_i;
    m = m_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = kl+ku+2; // as per API spec, ldab≥kl+ku+1.

#if LAPACKE_TEST_VERBOSE
   printf(" \n gbequ scomplex:  m: %d, n: %d lda: %d \n", m, n, lda);
#endif
    if(matrix_layout==LAPACK_COL_MAJOR){
        a_bufsize = lda*n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        lda = n;
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
       gbequ_free();
       EXPECT_FALSE( true) << "gbequ_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, a_bufsize);
    lapacke_gtest_init_float_buffer_pair_with_constant(r, rref, m, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(c, cref, n, 0.0);

   } /* end of Constructor  */

gbequ_scomplex_parameters:: ~gbequ_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbequ_scomplex_parameters object: destructor invoked. \n");
#endif
   gbequ_free();
}


//  Test fixture class definition
class cgbequ_test  : public  ::testing::Test {
public:
   gbequ_scomplex_parameters  *cgbequ_obj;
   void SetUp();  
   void TearDown () { delete cgbequ_obj; }
};


void cgbequ_test::SetUp(){

    /* LAPACKE CGBEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_cgbequ) ( int matrix_layout, lapack_int m,
                               lapack_int n, lapack_int kl, lapack_int ku,
                   const lapack_complex_float *a, 
                               lapack_int lda, float *r, float *c, 
                               float *rowcnd, float *colcnd, float *amax  );

    Fptr_NL_LAPACKE_cgbequ CGBEQU;

    cgbequ_obj = new  gbequ_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                            lin_solver_paramslist[idx].n,
                                            lin_solver_paramslist[idx].kl,
                                            lin_solver_paramslist[idx].ku
                                          );
    idx = Circular_Increment_Index(idx);

    cgbequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgbequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgbequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgbequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGBEQU = (Fptr_NL_LAPACKE_cgbequ)dlsym(cgbequ_obj->hModule, "LAPACKE_cgbequ");
    ASSERT_TRUE(CGBEQU != NULL) << "failed to get the Netlib LAPACKE_cgbequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    cgbequ_obj->inforef = CGBEQU( cgbequ_obj->matrix_layout,
                                  cgbequ_obj->m,
                                  cgbequ_obj->n,
                                  cgbequ_obj->kl,
                                  cgbequ_obj->ku,
                                  cgbequ_obj->aref,
                                  cgbequ_obj->lda,
                                  cgbequ_obj->rref,
                                  cgbequ_obj->cref,
                                  &cgbequ_obj->rowcndref,
                                  &cgbequ_obj->colcndref,
                                  &cgbequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    cgbequ_obj->info = LAPACKE_cgbequ( cgbequ_obj->matrix_layout,
                                  cgbequ_obj->m,
                                  cgbequ_obj->n,
                                  cgbequ_obj->kl,
                                  cgbequ_obj->ku,
                                  cgbequ_obj->a,
                                  cgbequ_obj->lda,
                                  cgbequ_obj->r,
                                  cgbequ_obj->c,
                                  &cgbequ_obj->rowcnd,
                                  &cgbequ_obj->colcnd,
                                  &cgbequ_obj->amax);

    if( cgbequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgbequ \
        is wrong\n", cgbequ_obj->info );
    }
    if( cgbequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgbequ is wrong\n",
        cgbequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cgbequ_obj->diff =  computeDiff_s(  cgbequ_obj->m,
                                        cgbequ_obj->r,
                                        cgbequ_obj->rref );
    cgbequ_obj->diff +=  computeDiff_s(  cgbequ_obj->n,
                                        cgbequ_obj->c,
                                        cgbequ_obj->cref );
}

TEST_F(cgbequ_test, cgbequ1) {
    EXPECT_NEAR(0.0, cgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->rowcnd, cgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->colcnd, cgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->amax, cgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbequ_test, cgbequ2) {
    EXPECT_NEAR(0.0, cgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->rowcnd, cgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->colcnd, cgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->amax, cgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbequ_test, cgbequ3) {
    EXPECT_NEAR(0.0, cgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->rowcnd, cgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->colcnd, cgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->amax, cgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbequ_test, cgbequ4) {
    EXPECT_NEAR(0.0, cgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->rowcnd, cgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->colcnd, cgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequ_obj->amax, cgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}


/* Begin gbequ_dcomplex_parameters  class definition */
class gbequ_dcomplex_parameters{
   public:
      int a_bufsize;
      double diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int m; // rows in A
      lapack_int n; // Columns in A
      lapack_int kl; // number of subdiagonals
      lapack_int ku; // number of superdiagonals
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
      gbequ_dcomplex_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i,
              lapack_int kl_i, lapack_int ku_i);
             
      ~gbequ_dcomplex_parameters ();
};  /* end of gbequ_dcomplex_parameters  class definition */


/* Constructor gbequ_dcomplex_parameters definition */
gbequ_dcomplex_parameters:: gbequ_dcomplex_parameters ( int matrix_layout_i, 
      lapack_int m_i, lapack_int n_i, lapack_int kl_i, lapack_int ku_i) {

    matrix_layout = matrix_layout_i;
    m = m_i;
    n = n_i;
    kl = kl_i;
    ku = ku_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = kl+ku+2; // as per API spec, ldab≥kl+ku+1.

#if LAPACKE_TEST_VERBOSE
   printf(" \n gbequ double:  m: %d, n: %d lda: %d \n", m, n, lda);
#endif
    if(matrix_layout==LAPACK_COL_MAJOR){
        a_bufsize = lda*n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        lda = n;
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
       gbequ_free();
       EXPECT_FALSE( true) << "gbequ_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, a_bufsize);
    lapacke_gtest_init_double_buffer_pair_with_constant(r, rref, m, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(c, cref, n, 0.0);

   } /* end of Constructor  */

gbequ_dcomplex_parameters:: ~gbequ_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbequ_dcomplex_parameters object: destructor invoked. \n");
#endif
   gbequ_free();
}


//  Test fixture class definition
class zgbequ_test  : public  ::testing::Test {
public:
   gbequ_dcomplex_parameters  *zgbequ_obj;
   void SetUp();  
   void TearDown () { delete zgbequ_obj; }
};


void zgbequ_test::SetUp(){

    /* LAPACKE ZGBEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_zgbequ) ( int matrix_layout, lapack_int m,
                               lapack_int n, lapack_int kl, lapack_int ku,
                   const lapack_complex_double *a, 
                               lapack_int lda, double *r, double *c, 
                               double *rowcnd, double *colcnd, double *amax  );

    Fptr_NL_LAPACKE_zgbequ ZGBEQU;

    zgbequ_obj = new  gbequ_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].kl,
                                         lin_solver_paramslist[idx].ku
                                          );
    idx = Circular_Increment_Index(idx);

    zgbequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgbequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgbequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgbequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGBEQU = (Fptr_NL_LAPACKE_zgbequ)dlsym(zgbequ_obj->hModule, "LAPACKE_zgbequ");
    ASSERT_TRUE(ZGBEQU != NULL) << "failed to get the Netlib LAPACKE_zgbequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zgbequ_obj->inforef = ZGBEQU( zgbequ_obj->matrix_layout,
                                  zgbequ_obj->m,
                                  zgbequ_obj->n,
                                  zgbequ_obj->kl,
                                  zgbequ_obj->ku,
                                  zgbequ_obj->aref,
                                  zgbequ_obj->lda,
                                  zgbequ_obj->rref,
                                  zgbequ_obj->cref,
                                  &zgbequ_obj->rowcndref,
                                  &zgbequ_obj->colcndref,
                                  &zgbequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    zgbequ_obj->info = LAPACKE_zgbequ( zgbequ_obj->matrix_layout,
                                  zgbequ_obj->m,
                                  zgbequ_obj->n,
                                  zgbequ_obj->kl,
                                  zgbequ_obj->ku,
                                  zgbequ_obj->a,
                                  zgbequ_obj->lda,
                                  zgbequ_obj->r,
                                  zgbequ_obj->c,
                                  &zgbequ_obj->rowcnd,
                                  &zgbequ_obj->colcnd,
                                  &zgbequ_obj->amax);

    if( zgbequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgbequ \
        is wrong\n", zgbequ_obj->info );
    }
    if( zgbequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgbequ is wrong\n",
        zgbequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zgbequ_obj->diff =  computeDiff_d(  zgbequ_obj->m,
                                        zgbequ_obj->r,
                                        zgbequ_obj->rref );
    zgbequ_obj->diff +=  computeDiff_d(  zgbequ_obj->n,
                                        zgbequ_obj->c,
                                        zgbequ_obj->cref );
}

TEST_F(zgbequ_test, zgbequ1) {
    EXPECT_NEAR(0.0, zgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->rowcnd, zgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->colcnd, zgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->amax, zgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbequ_test, zgbequ2) {
    EXPECT_NEAR(0.0, zgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->rowcnd, zgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->colcnd, zgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->amax, zgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbequ_test, zgbequ3) {
    EXPECT_NEAR(0.0, zgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->rowcnd, zgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->colcnd, zgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->amax, zgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbequ_test, zgbequ4) {
    EXPECT_NEAR(0.0, zgbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->rowcnd, zgbequ_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->colcnd, zgbequ_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequ_obj->amax, zgbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}
