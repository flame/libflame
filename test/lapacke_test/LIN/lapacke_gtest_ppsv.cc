#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define ppsv_free() \
       if (b != NULL)    free (b   ); \
       if (bref != NULL) free (bref); \
       if (a != NULL)    free (a  ); \
       if (aref != NULL) free (aref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ppsv_double_parameters  class definition */
class ppsv_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array ab contains the matrix A
      int b_bufsize;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ppsv_double_parameters ( int matrix_layout_i, char uplo_i, 
                                    lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~ppsv_double_parameters (); 
};  /* end of ppsv_double_parameters  class definition */


/* Constructor ppsv_double_parameters definition */
ppsv_double_parameters:: ppsv_double_parameters ( int matrix_layout_i, 
                              char uplo_i, lapack_int n_i,
                 lapack_int nrhs_i, lapack_int ldb_i) {
    int buf_siz;                                
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
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
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*(n + 1)/2) );
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "ppsv_double_parameters object: malloc error. Exiting ";
       ppsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

ppsv_double_parameters:: ~ppsv_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppsv_double_parameters object: destructor invoked. \n");
#endif
//   ppsv_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class dppsv_test  : public  ::testing::Test {
public:
   ppsv_double_parameters  *dppsv_obj;
   void SetUp();  
   void TearDown () { delete dppsv_obj; }
};

void dppsv_test::SetUp(){

    /* LAPACKE DPPSV prototype */
    typedef int (*Fptr_NL_LAPACKE_dppsv) ( int matrix_layout, 
         char uplo, lapack_int n, lapack_int nrhs,
               double * a, double * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_dppsv DPPSV;

    void *hModule, *dModule;
    dppsv_obj = new  ppsv_double_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    DPPSV = (Fptr_NL_LAPACKE_dppsv)dlsym(hModule, "LAPACKE_dppsv");
    ASSERT_TRUE(DPPSV != NULL) << "failed to get the Netlib LAPACKE_dppsv symbol";
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dppsv_obj->inforef = DPPSV( dppsv_obj->matrix_layout, 
          dppsv_obj->uplo, dppsv_obj->n,
         dppsv_obj->nrhs, ( double *)dppsv_obj->aref, 
                          dppsv_obj->bref, dppsv_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    dppsv_obj->info = LAPACKE_dppsv( dppsv_obj->matrix_layout, 
               dppsv_obj->uplo, dppsv_obj->n,
                 dppsv_obj->nrhs, ( double *)dppsv_obj->a, 
                                 dppsv_obj->b, dppsv_obj->ldb );

    if( dppsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dppsv is wrong\n", 
                    dppsv_obj->info );
    }
    if( dppsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dppsv is wrong\n", 
        dppsv_obj->inforef );
    }
    
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(dppsv_test, dppsv1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dppsv_obj->b_bufsize,
                           dppsv_obj->b, dppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dppsv_test, dppsv2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dppsv_obj->b_bufsize,
                           dppsv_obj->b, dppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dppsv_test, dppsv3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dppsv_obj->b_bufsize,
                           dppsv_obj->b, dppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dppsv_test, dppsv4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dppsv_obj->b_bufsize,
                           dppsv_obj->b, dppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin ppsv_float_parameters  class definition */
class ppsv_float_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array ab contains the matrix A
      int b_bufsize;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ppsv_float_parameters ( int matrix_layout_i, char uplo_i, 
                                    lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~ppsv_float_parameters (); 
};  /* end of ppsv_float_parameters  class definition */


/* Constructor ppsv_float_parameters definition */
ppsv_float_parameters:: ppsv_float_parameters ( int matrix_layout_i, 
                              char uplo_i, lapack_int n_i,
                 lapack_int nrhs_i, lapack_int ldb_i) {
    int buf_siz;                                
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
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
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*(n + 1)/2) );
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "ppsv_float_parameters object: malloc error. Exiting ";
       ppsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

ppsv_float_parameters:: ~ppsv_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppsv_float_parameters object: destructor invoked. \n");
#endif
//   ppsv_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class sppsv_test  : public  ::testing::Test {
public:
   ppsv_float_parameters  *sppsv_obj;
   void SetUp();  
   void TearDown () { delete sppsv_obj; }
};

void sppsv_test::SetUp(){

    /* LAPACKE SPPSV prototype */
    typedef int (*Fptr_NL_LAPACKE_sppsv) ( int matrix_layout, 
         char uplo, lapack_int n, lapack_int nrhs,
               float * a, float * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_sppsv SPPSV;
    
    void *hModule, *dModule;
    sppsv_obj = new  ppsv_float_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    SPPSV = (Fptr_NL_LAPACKE_sppsv)dlsym(hModule, "LAPACKE_sppsv");
    ASSERT_TRUE(SPPSV != NULL) << "failed to get the Netlib LAPACKE_sppsv symbol";
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    sppsv_obj->inforef = SPPSV( sppsv_obj->matrix_layout, 
          sppsv_obj->uplo, sppsv_obj->n,
         sppsv_obj->nrhs, ( float *)sppsv_obj->aref, 
                          sppsv_obj->bref, sppsv_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    sppsv_obj->info = LAPACKE_sppsv( sppsv_obj->matrix_layout, 
               sppsv_obj->uplo, sppsv_obj->n,
                 sppsv_obj->nrhs, ( float *)sppsv_obj->a, 
                                 sppsv_obj->b, sppsv_obj->ldb );

    if( sppsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sppsv is wrong\n", 
                    sppsv_obj->info );
    }
    if( sppsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sppsv is wrong\n", 
        sppsv_obj->inforef );
    }
    
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(sppsv_test, sppsv1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sppsv_obj->b_bufsize,
                           sppsv_obj->b, sppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(sppsv_test, sppsv2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sppsv_obj->b_bufsize,
                           sppsv_obj->b, sppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sppsv_test, sppsv3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sppsv_obj->b_bufsize,
                           sppsv_obj->b, sppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(sppsv_test, sppsv4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sppsv_obj->b_bufsize,
                           sppsv_obj->b, sppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin ppsv_scomplex_parameters  class definition */
class ppsv_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array ab contains the matrix A
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ppsv_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                    lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~ppsv_scomplex_parameters (); 
};  /* end of ppsv_scomplex_parameters  class definition */


/* Constructor ppsv_scomplex_parameters definition */
ppsv_scomplex_parameters:: ppsv_scomplex_parameters ( int matrix_layout_i, 
                              char uplo_i, lapack_int n_i,
                 lapack_int nrhs_i, lapack_int ldb_i) {
    int buf_siz;                                
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
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
        lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*(n + 1)/2) );
        lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "ppsv_scomplex_parameters object: malloc error. Exiting ";
       ppsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

ppsv_scomplex_parameters:: ~ppsv_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppsv_scomplex_parameters object: destructor invoked. \n");
#endif
//   ppsv_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class cppsv_test  : public  ::testing::Test {
public:
   ppsv_scomplex_parameters  *cppsv_obj;
   void SetUp();  
   void TearDown () { delete cppsv_obj; }
};

void cppsv_test::SetUp(){

    /* LAPACKE CPPSV prototype */
    typedef int (*Fptr_NL_LAPACKE_cppsv) ( int matrix_layout, 
         char uplo, lapack_int n, lapack_int nrhs,
               lapack_complex_float * a, lapack_complex_float * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_cppsv CPPSV;

    void *hModule, *dModule;
    cppsv_obj = new  ppsv_scomplex_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    CPPSV = (Fptr_NL_LAPACKE_cppsv)dlsym(hModule, "LAPACKE_cppsv");
    ASSERT_TRUE(CPPSV != NULL) << "failed to get the Netlib LAPACKE_cppsv symbol";
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    cppsv_obj->inforef = CPPSV( cppsv_obj->matrix_layout, 
          cppsv_obj->uplo, cppsv_obj->n,
         cppsv_obj->nrhs, ( lapack_complex_float *)cppsv_obj->aref, 
                          cppsv_obj->bref, cppsv_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    cppsv_obj->info = LAPACKE_cppsv( cppsv_obj->matrix_layout, 
               cppsv_obj->uplo, cppsv_obj->n,
                 cppsv_obj->nrhs, ( lapack_complex_float *)cppsv_obj->a, 
                                 cppsv_obj->b, cppsv_obj->ldb );

    if( cppsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cppsv is wrong\n", 
                    cppsv_obj->info );
    }
    if( cppsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cppsv is wrong\n", 
        cppsv_obj->inforef );
    }
    
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(cppsv_test, cppsv1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cppsv_obj->b_bufsize,
                           cppsv_obj->b, cppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cppsv_test, cppsv2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cppsv_obj->b_bufsize,
                           cppsv_obj->b, cppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cppsv_test, cppsv3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cppsv_obj->b_bufsize,
                           cppsv_obj->b, cppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cppsv_test, cppsv4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cppsv_obj->b_bufsize,
                           cppsv_obj->b, cppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin ppsv_dcomplex_parameters  class definition */
class ppsv_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array ab contains the matrix A
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      ppsv_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                    lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~ppsv_dcomplex_parameters (); 
};  /* end of ppsv_dcomplex_parameters  class definition */


/* Constructor ppsv_dcomplex_parameters definition */
ppsv_dcomplex_parameters:: ppsv_dcomplex_parameters ( int matrix_layout_i, 
                              char uplo_i, lapack_int n_i,
                 lapack_int nrhs_i, lapack_int ldb_i) {
    int buf_siz;                                
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
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
        lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*(n + 1)/2) );
        lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "ppsv_dcomplex_parameters object: malloc error. Exiting ";
       ppsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

ppsv_dcomplex_parameters:: ~ppsv_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppsv_dcomplex_parameters object: destructor invoked. \n");
#endif
//   ppsv_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class zppsv_test  : public  ::testing::Test {
public:
   ppsv_dcomplex_parameters  *zppsv_obj;
   void SetUp();  
   void TearDown () { delete zppsv_obj; }
};

void zppsv_test::SetUp(){

    /* LAPACKE ZPPSV prototype */
    typedef int (*Fptr_NL_LAPACKE_zppsv) ( int matrix_layout, 
         char uplo, lapack_int n, lapack_int nrhs,
               lapack_complex_double * a, lapack_complex_double * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_zppsv ZPPSV;

    void *hModule, *dModule;
    zppsv_obj = new  ppsv_dcomplex_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    ZPPSV = (Fptr_NL_LAPACKE_zppsv)dlsym(hModule, "LAPACKE_zppsv");
    ASSERT_TRUE(ZPPSV != NULL) << "failed to get the Netlib LAPACKE_zppsv symbol";
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    zppsv_obj->inforef = ZPPSV( zppsv_obj->matrix_layout, 
          zppsv_obj->uplo, zppsv_obj->n,
         zppsv_obj->nrhs, ( lapack_complex_double *)zppsv_obj->aref, 
                          zppsv_obj->bref, zppsv_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    zppsv_obj->info = LAPACKE_zppsv( zppsv_obj->matrix_layout, 
               zppsv_obj->uplo, zppsv_obj->n,
                 zppsv_obj->nrhs, ( lapack_complex_double *)zppsv_obj->a, 
                                 zppsv_obj->b, zppsv_obj->ldb );

    if( zppsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zppsv is wrong\n", 
                    zppsv_obj->info );
    }
    if( zppsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zppsv is wrong\n", 
        zppsv_obj->inforef );
    }
    
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(zppsv_test, zppsv1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zppsv_obj->b_bufsize,
                           zppsv_obj->b, zppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zppsv_test, zppsv2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zppsv_obj->b_bufsize,
                           zppsv_obj->b, zppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zppsv_test, zppsv3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zppsv_obj->b_bufsize,
                           zppsv_obj->b, zppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zppsv_test, zppsv4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zppsv_obj->b_bufsize,
                           zppsv_obj->b, zppsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}