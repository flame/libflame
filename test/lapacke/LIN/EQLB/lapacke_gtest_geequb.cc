#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define geequb_free() \
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

/* Begin geequb_float_parameters  class definition */
class geequb_float_parameters{
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
      geequb_float_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i);
             
      ~geequb_float_parameters ();
};  /* end of geequb_float_parameters  class definition */


/* Constructor geequb_float_parameters definition */
geequb_float_parameters:: geequb_float_parameters ( int matrix_layout_i, 
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
   printf(" \n geequb float:  m: %d, n: %d lda: %d \n", m, n, lda);
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
       geequb_free();
       EXPECT_FALSE( true) << "geequb_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, a_bufsize);

   } /* end of Constructor  */

geequb_float_parameters:: ~geequb_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" geequb_float_parameters object: destructor invoked. \n");
#endif
   geequb_free();
}


//  Test fixture class definition
class sgeequb_test  : public  ::testing::Test {
public:
   geequb_float_parameters  *sgeequb_obj;
   void SetUp();  
   void TearDown () { delete sgeequb_obj; }
};


void sgeequb_test::SetUp(){

    /* LAPACKE SGEEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeequb) ( int matrix_layout, lapack_int m,
                               lapack_int n, const float *a, lapack_int lda,
                               float *r, float *c, float *rowcnd, 
                               float *colcnd, float *amax  );

    Fptr_NL_LAPACKE_sgeequb SGEEQUB;

    sgeequb_obj = new  geequb_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    sgeequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGEEQUB = (Fptr_NL_LAPACKE_sgeequb)dlsym(sgeequb_obj->hModule, "LAPACKE_sgeequb");
    ASSERT_TRUE(SGEEQUB != NULL) << "failed to get the Netlib LAPACKE_sgeequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    sgeequb_obj->inforef = SGEEQUB( sgeequb_obj->matrix_layout,
                                  sgeequb_obj->m,
                                  sgeequb_obj->n,
                                  sgeequb_obj->aref,
                                  sgeequb_obj->lda,
                                  sgeequb_obj->rref,
                                  sgeequb_obj->cref,
                                  &sgeequb_obj->rowcndref,
                                  &sgeequb_obj->colcndref,
                                  &sgeequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    sgeequb_obj->info = LAPACKE_sgeequb( sgeequb_obj->matrix_layout,
                                  sgeequb_obj->m,
                                  sgeequb_obj->n,
                                  sgeequb_obj->a,
                                  sgeequb_obj->lda,
                                  sgeequb_obj->r,
                                  sgeequb_obj->c,
                                  &sgeequb_obj->rowcnd,
                                  &sgeequb_obj->colcnd,
                                  &sgeequb_obj->amax);

    if( sgeequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgeequb \
        is wrong\n", sgeequb_obj->info );
    }
    if( sgeequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeequb is wrong\n",
        sgeequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    sgeequb_obj->diff =  computeDiff_s(  sgeequb_obj->m,
                                        sgeequb_obj->r,
                                        sgeequb_obj->rref );
    sgeequb_obj->diff +=  computeDiff_s(  sgeequb_obj->n,
                                        sgeequb_obj->c,
                                        sgeequb_obj->cref );
}

TEST_F(sgeequb_test, sgeequb1) {
    EXPECT_NEAR(0.0, sgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->rowcnd, sgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->colcnd, sgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->amax, sgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeequb_test, sgeequb2) {
    EXPECT_NEAR(0.0, sgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->rowcnd, sgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->colcnd, sgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->amax, sgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeequb_test, sgeequb3) {
    EXPECT_NEAR(0.0, sgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->rowcnd, sgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->colcnd, sgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->amax, sgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeequb_test, sgeequb4) {
    EXPECT_NEAR(0.0, sgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->rowcnd, sgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->colcnd, sgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sgeequb_obj->amax, sgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin geequb_double_parameters  class definition */
class geequb_double_parameters{
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
      geequb_double_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i);
             
      ~geequb_double_parameters ();
};  /* end of geequb_double_parameters  class definition */


/* Constructor geequb_double_parameters definition */
geequb_double_parameters:: geequb_double_parameters ( int matrix_layout_i, 
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
   printf(" \n geequb double:  m: %d, n: %d lda: %d \n", m, n, lda);
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
       geequb_free();
       EXPECT_FALSE( true) << "geequb_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, a_bufsize);

   } /* end of Constructor  */

geequb_double_parameters:: ~geequb_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" geequb_double_parameters object: destructor invoked. \n");
#endif
   geequb_free();
}


//  Test fixture class definition
class dgeequb_test  : public  ::testing::Test {
public:
   geequb_double_parameters  *dgeequb_obj;
   void SetUp();  
   void TearDown () { delete dgeequb_obj; }
};


void dgeequb_test::SetUp(){

    /* LAPACKE DGEEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeequb) ( int matrix_layout, lapack_int m,
                               lapack_int n, const double *a, lapack_int lda,
                               double *r, double *c, double *rowcnd, 
                               double *colcnd, double *amax  );

    Fptr_NL_LAPACKE_dgeequb DGEEQUB;

    dgeequb_obj = new  geequb_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    dgeequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGEEQUB = (Fptr_NL_LAPACKE_dgeequb)dlsym(dgeequb_obj->hModule, "LAPACKE_dgeequb");
    ASSERT_TRUE(DGEEQUB != NULL) << "failed to get the Netlib LAPACKE_dgeequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    dgeequb_obj->inforef = DGEEQUB( dgeequb_obj->matrix_layout,
                                  dgeequb_obj->m,
                                  dgeequb_obj->n,
                                  dgeequb_obj->aref,
                                  dgeequb_obj->lda,
                                  dgeequb_obj->rref,
                                  dgeequb_obj->cref,
                                  &dgeequb_obj->rowcndref,
                                  &dgeequb_obj->colcndref,
                                  &dgeequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    dgeequb_obj->info = LAPACKE_dgeequb( dgeequb_obj->matrix_layout,
                                  dgeequb_obj->m,
                                  dgeequb_obj->n,
                                  dgeequb_obj->a,
                                  dgeequb_obj->lda,
                                  dgeequb_obj->r,
                                  dgeequb_obj->c,
                                  &dgeequb_obj->rowcnd,
                                  &dgeequb_obj->colcnd,
                                  &dgeequb_obj->amax);

    if( dgeequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dgeequb \
        is wrong\n", dgeequb_obj->info );
    }
    if( dgeequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeequb is wrong\n",
        dgeequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dgeequb_obj->diff =  computeDiff_d(  dgeequb_obj->m,
                                        dgeequb_obj->r,
                                        dgeequb_obj->rref );
    dgeequb_obj->diff +=  computeDiff_d(  dgeequb_obj->n,
                                        dgeequb_obj->c,
                                        dgeequb_obj->cref );
}

TEST_F(dgeequb_test, dgeequb1) {
    EXPECT_NEAR(0.0, dgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->rowcnd, dgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->colcnd, dgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->amax, dgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeequb_test, dgeequb2) {
    EXPECT_NEAR(0.0, dgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->rowcnd, dgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->colcnd, dgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->amax, dgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeequb_test, dgeequb3) {
    EXPECT_NEAR(0.0, dgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->rowcnd, dgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->colcnd, dgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->amax, dgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeequb_test, dgeequb4) {
    EXPECT_NEAR(0.0, dgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->rowcnd, dgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->colcnd, dgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dgeequb_obj->amax, dgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin geequb_scomplex_parameters  class definition */
class geequb_scomplex_parameters{
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
      geequb_scomplex_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i);
             
      ~geequb_scomplex_parameters ();
};  /* end of geequb_scomplex_parameters  class definition */


/* Constructor geequb_scomplex_parameters definition */
geequb_scomplex_parameters:: geequb_scomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n geequb scomplex:  m: %d, n: %d lda: %d \n", m, n, lda);
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
       geequb_free();
       EXPECT_FALSE( true) << "geequb_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, a_bufsize);

   } /* end of Constructor  */

geequb_scomplex_parameters:: ~geequb_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" geequb_scomplex_parameters object: destructor invoked. \n");
#endif
   geequb_free();
}


//  Test fixture class definition
class cgeequb_test  : public  ::testing::Test {
public:
   geequb_scomplex_parameters  *cgeequb_obj;
   void SetUp();  
   void TearDown () { delete cgeequb_obj; }
};


void cgeequb_test::SetUp(){

    /* LAPACKE CGEEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeequb) ( int matrix_layout, lapack_int m,
                               lapack_int n, const lapack_complex_float *a, 
                               lapack_int lda, float *r, float *c, 
                               float *rowcnd, float *colcnd, float *amax  );

    Fptr_NL_LAPACKE_cgeequb CGEEQUB;

    cgeequb_obj = new  geequb_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    cgeequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgeequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgeequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgeequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGEEQUB = (Fptr_NL_LAPACKE_cgeequb)dlsym(cgeequb_obj->hModule, "LAPACKE_cgeequb");
    ASSERT_TRUE(CGEEQUB != NULL) << "failed to get the Netlib LAPACKE_cgeequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    cgeequb_obj->inforef = CGEEQUB( cgeequb_obj->matrix_layout,
                                  cgeequb_obj->m,
                                  cgeequb_obj->n,
                                  cgeequb_obj->aref,
                                  cgeequb_obj->lda,
                                  cgeequb_obj->rref,
                                  cgeequb_obj->cref,
                                  &cgeequb_obj->rowcndref,
                                  &cgeequb_obj->colcndref,
                                  &cgeequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    cgeequb_obj->info = LAPACKE_cgeequb( cgeequb_obj->matrix_layout,
                                  cgeequb_obj->m,
                                  cgeequb_obj->n,
                                  cgeequb_obj->a,
                                  cgeequb_obj->lda,
                                  cgeequb_obj->r,
                                  cgeequb_obj->c,
                                  &cgeequb_obj->rowcnd,
                                  &cgeequb_obj->colcnd,
                                  &cgeequb_obj->amax);

    if( cgeequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cgeequb \
        is wrong\n", cgeequb_obj->info );
    }
    if( cgeequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeequb is wrong\n",
        cgeequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cgeequb_obj->diff =  computeDiff_s(  cgeequb_obj->m,
                                        cgeequb_obj->r,
                                        cgeequb_obj->rref );
    cgeequb_obj->diff +=  computeDiff_s(  cgeequb_obj->n,
                                        cgeequb_obj->c,
                                        cgeequb_obj->cref );
}

TEST_F(cgeequb_test, cgeequb1) {
    EXPECT_NEAR(0.0, cgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->rowcnd, cgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->colcnd, cgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->amax, cgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeequb_test, cgeequb2) {
    EXPECT_NEAR(0.0, cgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->rowcnd, cgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->colcnd, cgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->amax, cgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeequb_test, cgeequb3) {
    EXPECT_NEAR(0.0, cgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->rowcnd, cgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->colcnd, cgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->amax, cgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeequb_test, cgeequb4) {
    EXPECT_NEAR(0.0, cgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->rowcnd, cgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->colcnd, cgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cgeequb_obj->amax, cgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}


/* Begin geequb_dcomplex_parameters  class definition */
class geequb_dcomplex_parameters{
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
      geequb_dcomplex_parameters ( int matrix_layout_i, 
              lapack_int m_i, lapack_int n_i);
             
      ~geequb_dcomplex_parameters ();
};  /* end of geequb_dcomplex_parameters  class definition */


/* Constructor geequb_dcomplex_parameters definition */
geequb_dcomplex_parameters:: geequb_dcomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n geequb dcomplex:  m: %d, n: %d lda: %d \n", m, n, lda);
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
       geequb_free();
       EXPECT_FALSE( true) << "geequb_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, a_bufsize);

   } /* end of Constructor  */

geequb_dcomplex_parameters:: ~geequb_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" geequb_dcomplex_parameters object: destructor invoked. \n");
#endif
   geequb_free();
}


//  Test fixture class definition
class zgeequb_test  : public  ::testing::Test {
public:
   geequb_dcomplex_parameters  *zgeequb_obj;
   void SetUp();  
   void TearDown () { delete zgeequb_obj; }
};


void zgeequb_test::SetUp(){

    /* LAPACKE ZGEEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeequb) ( int matrix_layout, lapack_int m,
                               lapack_int n, const lapack_complex_double *a, 
                               lapack_int lda, double *r, double *c, 
                               double *rowcnd, double *colcnd, double *amax  );

    Fptr_NL_LAPACKE_zgeequb ZGEEQUB;

    zgeequb_obj = new  geequb_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].m,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    zgeequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgeequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgeequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgeequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGEEQUB = (Fptr_NL_LAPACKE_zgeequb)dlsym(zgeequb_obj->hModule, "LAPACKE_zgeequb");
    ASSERT_TRUE(ZGEEQUB != NULL) << "failed to get the Netlib LAPACKE_zgeequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zgeequb_obj->inforef = ZGEEQUB( zgeequb_obj->matrix_layout,
                                  zgeequb_obj->m,
                                  zgeequb_obj->n,
                                  zgeequb_obj->aref,
                                  zgeequb_obj->lda,
                                  zgeequb_obj->rref,
                                  zgeequb_obj->cref,
                                  &zgeequb_obj->rowcndref,
                                  &zgeequb_obj->colcndref,
                                  &zgeequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    zgeequb_obj->info = LAPACKE_zgeequb( zgeequb_obj->matrix_layout,
                                  zgeequb_obj->m,
                                  zgeequb_obj->n,
                                  zgeequb_obj->a,
                                  zgeequb_obj->lda,
                                  zgeequb_obj->r,
                                  zgeequb_obj->c,
                                  &zgeequb_obj->rowcnd,
                                  &zgeequb_obj->colcnd,
                                  &zgeequb_obj->amax);

    if( zgeequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zgeequb \
        is wrong\n", zgeequb_obj->info );
    }
    if( zgeequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeequb is wrong\n",
        zgeequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zgeequb_obj->diff =  computeDiff_d(  zgeequb_obj->m,
                                        zgeequb_obj->r,
                                        zgeequb_obj->rref );
    zgeequb_obj->diff +=  computeDiff_d(  zgeequb_obj->n,
                                        zgeequb_obj->c,
                                        zgeequb_obj->cref );
}

TEST_F(zgeequb_test, zgeequb1) {
    EXPECT_NEAR(0.0, zgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->rowcnd, zgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->colcnd, zgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->amax, zgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeequb_test, zgeequb2) {
    EXPECT_NEAR(0.0, zgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->rowcnd, zgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->colcnd, zgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->amax, zgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeequb_test, zgeequb3) {
    EXPECT_NEAR(0.0, zgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->rowcnd, zgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->colcnd, zgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->amax, zgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeequb_test, zgeequb4) {
    EXPECT_NEAR(0.0, zgeequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->rowcnd, zgeequb_obj->rowcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->colcnd, zgeequb_obj->colcndref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zgeequb_obj->amax, zgeequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}
