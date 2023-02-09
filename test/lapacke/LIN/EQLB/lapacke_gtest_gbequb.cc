#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define gbequb_free() \
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

/* Begin gbequb_float_parameters  class definition */
class gbequb_float_parameters{
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
      gbequb_float_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i,
              lapack_int kl_i, lapack_int ku_i);
             
      ~gbequb_float_parameters ();
};  /* end of gbequb_float_parameters  class definition */


/* Constructor gbequb_float_parameters definition */
gbequb_float_parameters:: gbequb_float_parameters ( int matrix_layout_i, 
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

    lda = m; // as per API spec, lda≥ max(1, m).

#if LAPACKE_TEST_VERBOSE
   printf(" \n gbequb float:  m: %d, n: %d lda: %d \n", m, n, lda);
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
       gbequb_free();
       EXPECT_FALSE( true) << "gbequb_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, a_bufsize);
    lapacke_gtest_init_float_buffer_pair_with_constant(r, rref, m, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(c, cref, n, 0.0);

   } /* end of Constructor  */

gbequb_float_parameters:: ~gbequb_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbequb_float_parameters object: destructor invoked. \n");
#endif
   gbequb_free();
}


//  Test fixture class definition
class sgbequb_test  : public  ::testing::Test {
public:
   gbequb_float_parameters  *sgbequb_obj;
   void SetUp();  
   void TearDown () { delete sgbequb_obj; }
};


void sgbequb_test::SetUp(){

    /* LAPACKE SGBEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_sgbequb) ( int matrix_layout, lapack_int m,
                               lapack_int n, lapack_int kl, lapack_int ku,
                               const float *a, lapack_int lda,
                               float *r, float *c, float *rowcnd, 
                               float *colcnd, float *amax  );

    Fptr_NL_LAPACKE_sgbequb SGBEQUB;

    sgbequb_obj = new  gbequb_float_parameters(
                                lin_solver_paramslist[idx].matrix_layout,
                                            lin_solver_paramslist[idx].m,
                                            lin_solver_paramslist[idx].n,
                                            lin_solver_paramslist[idx].kl,
                                            lin_solver_paramslist[idx].ku
                                          );
    idx = Circular_Increment_Index(idx);

    sgbequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgbequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgbequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgbequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGBEQUB = (Fptr_NL_LAPACKE_sgbequb)dlsym(sgbequb_obj->hModule, "LAPACKE_sgbequb");
    ASSERT_TRUE(SGBEQUB != NULL) << "failed to get the Netlib LAPACKE_sgbequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    sgbequb_obj->inforef = SGBEQUB( sgbequb_obj->matrix_layout,
                                  sgbequb_obj->m,
                                  sgbequb_obj->n,
                                  sgbequb_obj->kl,
                                  sgbequb_obj->ku,
                                  sgbequb_obj->aref,
                                  sgbequb_obj->lda,
                                  sgbequb_obj->rref,
                                  sgbequb_obj->cref,
                                  &sgbequb_obj->rowcndref,
                                  &sgbequb_obj->colcndref,
                                  &sgbequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    sgbequb_obj->info = LAPACKE_sgbequb( sgbequb_obj->matrix_layout,
                                  sgbequb_obj->m,
                                  sgbequb_obj->n,
                                  sgbequb_obj->kl,
                                  sgbequb_obj->ku,
                                  sgbequb_obj->a,
                                  sgbequb_obj->lda,
                                  sgbequb_obj->r,
                                  sgbequb_obj->c,
                                  &sgbequb_obj->rowcnd,
                                  &sgbequb_obj->colcnd,
                                  &sgbequb_obj->amax);

    if( sgbequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgbequb \
        is wrong\n", sgbequb_obj->info );
    }
    if( sgbequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgbequb is wrong\n",
        sgbequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    sgbequb_obj->diff =  computeDiff_s(  sgbequb_obj->m,
                                        sgbequb_obj->r,
                                        sgbequb_obj->rref );
    sgbequb_obj->diff +=  computeDiff_s(  sgbequb_obj->n,
                                        sgbequb_obj->c,
                                        sgbequb_obj->cref );
}

TEST_F(sgbequb_test, sgbequb1) {
    EXPECT_NEAR(0.0, sgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->rowcnd, sgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->colcnd, sgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->amax, sgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbequb_test, sgbequb2) {
    EXPECT_NEAR(0.0, sgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->rowcnd, sgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->colcnd, sgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->amax, sgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbequb_test, sgbequb3) {
    EXPECT_NEAR(0.0, sgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->rowcnd, sgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->colcnd, sgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->amax, sgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgbequb_test, sgbequb4) {
    EXPECT_NEAR(0.0, sgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->rowcnd, sgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->colcnd, sgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgbequb_obj->amax, sgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gbequb_double_parameters  class definition */
class gbequb_double_parameters{
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
      gbequb_double_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i,
              lapack_int kl_i, lapack_int ku_i);
             
      ~gbequb_double_parameters ();
};  /* end of gbequb_double_parameters  class definition */


/* Constructor gbequb_double_parameters definition */
gbequb_double_parameters:: gbequb_double_parameters ( int matrix_layout_i, 
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

    lda = m; // as per API spec, lda≥ max(1, m).

#if LAPACKE_TEST_VERBOSE
   printf(" \n gbequb double:  m: %d, n: %d lda: %d \n", m, n, lda);
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
       gbequb_free();
       EXPECT_FALSE( true) << "gbequb_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, a_bufsize);
    lapacke_gtest_init_double_buffer_pair_with_constant(r, rref, m, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(c, cref, n, 0.0);

   } /* end of Constructor  */

gbequb_double_parameters:: ~gbequb_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbequb_double_parameters object: destructor invoked. \n");
#endif
   gbequb_free();
}


//  Test fixture class definition
class dgbequb_test  : public  ::testing::Test {
public:
   gbequb_double_parameters  *dgbequb_obj;
   void SetUp();  
   void TearDown () { delete dgbequb_obj; }
};


void dgbequb_test::SetUp(){

    /* LAPACKE DGBEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_dgbequb) ( int matrix_layout, lapack_int m,
                               lapack_int n, lapack_int kl, lapack_int ku,
                   const double *a, lapack_int lda,
                               double *r, double *c, double *rowcnd, 
                               double *colcnd, double *amax  );

    Fptr_NL_LAPACKE_dgbequb DGBEQUB;

    dgbequb_obj = new  gbequb_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                            lin_solver_paramslist[idx].n,
                                            lin_solver_paramslist[idx].kl,
                                            lin_solver_paramslist[idx].ku
                                          );
    idx = Circular_Increment_Index(idx);

    dgbequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgbequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgbequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgbequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGBEQUB = (Fptr_NL_LAPACKE_dgbequb)dlsym(dgbequb_obj->hModule, "LAPACKE_dgbequb");
    ASSERT_TRUE(DGBEQUB != NULL) << "failed to get the Netlib LAPACKE_dgbequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    dgbequb_obj->inforef = DGBEQUB( dgbequb_obj->matrix_layout,
                                  dgbequb_obj->m,
                                  dgbequb_obj->n,
                                  dgbequb_obj->kl,
                                  dgbequb_obj->ku,
                                  dgbequb_obj->aref,
                                  dgbequb_obj->lda,
                                  dgbequb_obj->rref,
                                  dgbequb_obj->cref,
                                  &dgbequb_obj->rowcndref,
                                  &dgbequb_obj->colcndref,
                                  &dgbequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    dgbequb_obj->info = LAPACKE_dgbequb( dgbequb_obj->matrix_layout,
                                  dgbequb_obj->m,
                                  dgbequb_obj->n,
                                  dgbequb_obj->kl,
                                  dgbequb_obj->ku,
                                  dgbequb_obj->a,
                                  dgbequb_obj->lda,
                                  dgbequb_obj->r,
                                  dgbequb_obj->c,
                                  &dgbequb_obj->rowcnd,
                                  &dgbequb_obj->colcnd,
                                  &dgbequb_obj->amax);

    if( dgbequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgbequb \
        is wrong\n", dgbequb_obj->info );
    }
    if( dgbequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgbequb is wrong\n",
        dgbequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dgbequb_obj->diff =  computeDiff_d(  dgbequb_obj->m,
                                        dgbequb_obj->r,
                                        dgbequb_obj->rref );
    dgbequb_obj->diff +=  computeDiff_d(  dgbequb_obj->n,
                                        dgbequb_obj->c,
                                        dgbequb_obj->cref );
}

TEST_F(dgbequb_test, dgbequb1) {
    EXPECT_NEAR(0.0, dgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->rowcnd, dgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->colcnd, dgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->amax, dgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbequb_test, dgbequb2) {
    EXPECT_NEAR(0.0, dgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->rowcnd, dgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->colcnd, dgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->amax, dgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbequb_test, dgbequb3) {
    EXPECT_NEAR(0.0, dgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->rowcnd, dgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->colcnd, dgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->amax, dgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgbequb_test, dgbequb4) {
    EXPECT_NEAR(0.0, dgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->rowcnd, dgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->colcnd, dgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgbequb_obj->amax, dgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin gbequb_scomplex_parameters  class definition */
class gbequb_scomplex_parameters{
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
      gbequb_scomplex_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i,
              lapack_int kl_i, lapack_int ku_i);
             
      ~gbequb_scomplex_parameters ();
};  /* end of gbequb_scomplex_parameters  class definition */


/* Constructor gbequb_scomplex_parameters definition */
gbequb_scomplex_parameters:: gbequb_scomplex_parameters ( int matrix_layout_i, 
      lapack_int m_i, lapack_int n_i, lapack_int kl_i, lapack_int ku_i) {

    matrix_layout = matrix_layout_i;
    m = m_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = m; // as per API spec, lda≥ max(1, m).

#if LAPACKE_TEST_VERBOSE
   printf(" \n gbequb scomplex:  m: %d, n: %d lda: %d \n", m, n, lda);
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
       gbequb_free();
       EXPECT_FALSE( true) << "gbequb_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, a_bufsize);
    lapacke_gtest_init_float_buffer_pair_with_constant(r, rref, m, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(c, cref, n, 0.0);

   } /* end of Constructor  */

gbequb_scomplex_parameters:: ~gbequb_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbequb_scomplex_parameters object: destructor invoked. \n");
#endif
   gbequb_free();
}


//  Test fixture class definition
class cgbequb_test  : public  ::testing::Test {
public:
   gbequb_scomplex_parameters  *cgbequb_obj;
   void SetUp();  
   void TearDown () { delete cgbequb_obj; }
};


void cgbequb_test::SetUp(){

    /* LAPACKE CGBEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_cgbequb) ( int matrix_layout, lapack_int m,
                               lapack_int n, lapack_int kl, lapack_int ku,
                   const lapack_complex_float *a, 
                               lapack_int lda, float *r, float *c, 
                               float *rowcnd, float *colcnd, float *amax  );

    Fptr_NL_LAPACKE_cgbequb CGBEQUB;

    cgbequb_obj = new  gbequb_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                            lin_solver_paramslist[idx].n,
                                            lin_solver_paramslist[idx].kl,
                                            lin_solver_paramslist[idx].ku
                                          );
    idx = Circular_Increment_Index(idx);

    cgbequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgbequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgbequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgbequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGBEQUB = (Fptr_NL_LAPACKE_cgbequb)dlsym(cgbequb_obj->hModule, "LAPACKE_cgbequb");
    ASSERT_TRUE(CGBEQUB != NULL) << "failed to get the Netlib LAPACKE_cgbequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    cgbequb_obj->inforef = CGBEQUB( cgbequb_obj->matrix_layout,
                                  cgbequb_obj->m,
                                  cgbequb_obj->n,
                                  cgbequb_obj->kl,
                                  cgbequb_obj->ku,
                                  cgbequb_obj->aref,
                                  cgbequb_obj->lda,
                                  cgbequb_obj->rref,
                                  cgbequb_obj->cref,
                                  &cgbequb_obj->rowcndref,
                                  &cgbequb_obj->colcndref,
                                  &cgbequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    cgbequb_obj->info = LAPACKE_cgbequb( cgbequb_obj->matrix_layout,
                                  cgbequb_obj->m,
                                  cgbequb_obj->n,
                                  cgbequb_obj->kl,
                                  cgbequb_obj->ku,
                                  cgbequb_obj->a,
                                  cgbequb_obj->lda,
                                  cgbequb_obj->r,
                                  cgbequb_obj->c,
                                  &cgbequb_obj->rowcnd,
                                  &cgbequb_obj->colcnd,
                                  &cgbequb_obj->amax);

    if( cgbequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgbequb \
        is wrong\n", cgbequb_obj->info );
    }
    if( cgbequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgbequb is wrong\n",
        cgbequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cgbequb_obj->diff =  computeDiff_s(  cgbequb_obj->m,
                                        cgbequb_obj->r,
                                        cgbequb_obj->rref );
    cgbequb_obj->diff +=  computeDiff_s(  cgbequb_obj->n,
                                        cgbequb_obj->c,
                                        cgbequb_obj->cref );
}

TEST_F(cgbequb_test, cgbequb1) {
    EXPECT_NEAR(0.0, cgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->rowcnd, cgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->colcnd, cgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->amax, cgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbequb_test, cgbequb2) {
    EXPECT_NEAR(0.0, cgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->rowcnd, cgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->colcnd, cgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->amax, cgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbequb_test, cgbequb3) {
    EXPECT_NEAR(0.0, cgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->rowcnd, cgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->colcnd, cgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->amax, cgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgbequb_test, cgbequb4) {
    EXPECT_NEAR(0.0, cgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->rowcnd, cgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->colcnd, cgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgbequb_obj->amax, cgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}


/* Begin gbequb_dcomplex_parameters  class definition */
class gbequb_dcomplex_parameters{
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
      gbequb_dcomplex_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i,
              lapack_int kl_i, lapack_int ku_i);
             
      ~gbequb_dcomplex_parameters ();
};  /* end of gbequb_dcomplex_parameters  class definition */


/* Constructor gbequb_dcomplex_parameters definition */
gbequb_dcomplex_parameters:: gbequb_dcomplex_parameters ( int matrix_layout_i, 
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

    lda = m; // as per API spec, lda≥ max(1, m).

#if LAPACKE_TEST_VERBOSE
   printf(" \n gbequb double:  m: %d, n: %d lda: %d \n", m, n, lda);
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
       gbequb_free();
       EXPECT_FALSE( true) << "gbequb_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, a_bufsize);
    lapacke_gtest_init_double_buffer_pair_with_constant(r, rref, m, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(c, cref, n, 0.0);

   } /* end of Constructor  */

gbequb_dcomplex_parameters:: ~gbequb_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbequb_dcomplex_parameters object: destructor invoked. \n");
#endif
   gbequb_free();
}


//  Test fixture class definition
class zgbequb_test  : public  ::testing::Test {
public:
   gbequb_dcomplex_parameters  *zgbequb_obj;
   void SetUp();  
   void TearDown () { delete zgbequb_obj; }
};


void zgbequb_test::SetUp(){

    /* LAPACKE ZGBEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_zgbequb) ( int matrix_layout, lapack_int m,
                               lapack_int n, lapack_int kl, lapack_int ku,
                   const lapack_complex_double *a, 
                               lapack_int lda, double *r, double *c, 
                               double *rowcnd, double *colcnd, double *amax  );

    Fptr_NL_LAPACKE_zgbequb ZGBEQUB;

    zgbequb_obj = new  gbequb_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].kl,
                                         lin_solver_paramslist[idx].ku
                                          );
    idx = Circular_Increment_Index(idx);

    zgbequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgbequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgbequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgbequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGBEQUB = (Fptr_NL_LAPACKE_zgbequb)dlsym(zgbequb_obj->hModule, "LAPACKE_zgbequb");
    ASSERT_TRUE(ZGBEQUB != NULL) << "failed to get the Netlib LAPACKE_zgbequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zgbequb_obj->inforef = ZGBEQUB( zgbequb_obj->matrix_layout,
                                  zgbequb_obj->m,
                                  zgbequb_obj->n,
                                  zgbequb_obj->kl,
                                  zgbequb_obj->ku,
                                  zgbequb_obj->aref,
                                  zgbequb_obj->lda,
                                  zgbequb_obj->rref,
                                  zgbequb_obj->cref,
                                  &zgbequb_obj->rowcndref,
                                  &zgbequb_obj->colcndref,
                                  &zgbequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    zgbequb_obj->info = LAPACKE_zgbequb( zgbequb_obj->matrix_layout,
                                  zgbequb_obj->m,
                                  zgbequb_obj->n,
                                  zgbequb_obj->kl,
                                  zgbequb_obj->ku,
                                  zgbequb_obj->a,
                                  zgbequb_obj->lda,
                                  zgbequb_obj->r,
                                  zgbequb_obj->c,
                                  &zgbequb_obj->rowcnd,
                                  &zgbequb_obj->colcnd,
                                  &zgbequb_obj->amax);

    if( zgbequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgbequb \
        is wrong\n", zgbequb_obj->info );
    }
    if( zgbequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgbequb is wrong\n",
        zgbequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zgbequb_obj->diff =  computeDiff_d(  zgbequb_obj->m,
                                        zgbequb_obj->r,
                                        zgbequb_obj->rref );
    zgbequb_obj->diff +=  computeDiff_d(  zgbequb_obj->n,
                                        zgbequb_obj->c,
                                        zgbequb_obj->cref );
}

TEST_F(zgbequb_test, zgbequb1) {
    EXPECT_NEAR(0.0, zgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->rowcnd, zgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->colcnd, zgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->amax, zgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbequb_test, zgbequb2) {
    EXPECT_NEAR(0.0, zgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->rowcnd, zgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->colcnd, zgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->amax, zgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbequb_test, zgbequb3) {
    EXPECT_NEAR(0.0, zgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->rowcnd, zgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->colcnd, zgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->amax, zgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgbequb_test, zgbequb4) {
    EXPECT_NEAR(0.0, zgbequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->rowcnd, zgbequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->colcnd, zgbequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgbequb_obj->amax, zgbequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}
