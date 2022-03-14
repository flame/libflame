#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define pptrs_free() \
       if (b != NULL)    free (b   ); \
       if (bref != NULL) free (bref); \
       if (a != NULL)    free (a  ); \
       if (aref != NULL) free (aref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin pptrs_double_parameters  class definition */
class pptrs_double_parameters{
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
      pptrs_double_parameters ( int matrix_layout_i, char uplo_i, 
                                    lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~pptrs_double_parameters (); 
};  /* end of pptrs_double_parameters  class definition */


/* Constructor pptrs_double_parameters definition */
pptrs_double_parameters:: pptrs_double_parameters ( int matrix_layout_i, 
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
       EXPECT_FALSE( true) << "pptrs_double_parameters object: malloc error. Exiting ";
       pptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

pptrs_double_parameters:: ~pptrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptrs_double_parameters object: destructor invoked. \n");
#endif
//   pptrs_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class dpptrs_test  : public  ::testing::Test {
public:
   pptrs_double_parameters  *dpptrs_obj;
   void SetUp();  
   void TearDown () { delete dpptrs_obj; }
};

void dpptrs_test::SetUp(){

    /* LAPACKE DPPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dpptrs) ( int matrix_layout, 
         char uplo, lapack_int n, lapack_int nrhs,
              const double * a, double * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_dpptrs DPPTRS;
    
     /* LAPACKE DPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpptrf) ( int matrix_layout, 
           char uplo , lapack_int n , double * a );

    Fptr_NL_LAPACKE_dpptrf DPPTRF;

    void *hModule, *dModule;
    dpptrs_obj = new  pptrs_double_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    DPPTRS = (Fptr_NL_LAPACKE_dpptrs)dlsym(hModule, "LAPACKE_dpptrs");
    ASSERT_TRUE(DPPTRS != NULL) << "failed to get the Netlib LAPACKE_dpptrs symbol";

    DPPTRF = (Fptr_NL_LAPACKE_dpptrf)dlsym(hModule,"LAPACKE_dpptrf");
    ASSERT_TRUE(DPPTRF != NULL) << "failed to get the Netlib LAPACKE_dpptrf symbol";

    /* Pre condition: need to call pptrf - before calling pptrs function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dpptrs_obj->inforef = DPPTRF( dpptrs_obj->matrix_layout, 
          dpptrs_obj->uplo, dpptrs_obj->n,
         dpptrs_obj->aref);

    dpptrs_obj->inforef = DPPTRS( dpptrs_obj->matrix_layout, 
          dpptrs_obj->uplo, dpptrs_obj->n,
         dpptrs_obj->nrhs, (const double *)dpptrs_obj->aref, 
                          dpptrs_obj->bref, dpptrs_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    
    dpptrs_obj->info = LAPACKE_dpptrf( dpptrs_obj->matrix_layout, 
               dpptrs_obj->uplo, dpptrs_obj->n,
                 dpptrs_obj->a); 

    dpptrs_obj->info = LAPACKE_dpptrs( dpptrs_obj->matrix_layout, 
               dpptrs_obj->uplo, dpptrs_obj->n,
                 dpptrs_obj->nrhs, (const double *)dpptrs_obj->a, 
                                 dpptrs_obj->b, dpptrs_obj->ldb );

    if( dpptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpptrs is wrong\n", 
                    dpptrs_obj->info );
    }
    if( dpptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpptrs is wrong\n", 
        dpptrs_obj->inforef );
    }
    
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(dpptrs_test, dpptrs1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpptrs_obj->b_bufsize,
                           dpptrs_obj->b, dpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dpptrs_test, dpptrs2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpptrs_obj->b_bufsize,
                           dpptrs_obj->b, dpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpptrs_test, dpptrs3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpptrs_obj->b_bufsize,
                           dpptrs_obj->b, dpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dpptrs_test, dpptrs4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpptrs_obj->b_bufsize,
                           dpptrs_obj->b, dpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pptrs_float_parameters  class definition */
class pptrs_float_parameters{
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
      pptrs_float_parameters ( int matrix_layout_i, char uplo_i, 
                                    lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~pptrs_float_parameters (); 
};  /* end of pptrs_float_parameters  class definition */


/* Constructor pptrs_float_parameters definition */
pptrs_float_parameters:: pptrs_float_parameters ( int matrix_layout_i, 
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
       EXPECT_FALSE( true) << "pptrs_float_parameters object: malloc error. Exiting ";
       pptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

pptrs_float_parameters:: ~pptrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptrs_float_parameters object: destructor invoked. \n");
#endif
//   pptrs_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class spptrs_test  : public  ::testing::Test {
public:
   pptrs_float_parameters  *spptrs_obj;
   void SetUp();  
   void TearDown () { delete spptrs_obj; }
};

void spptrs_test::SetUp(){

    /* LAPACKE SPPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_spptrs) ( int matrix_layout, 
         char uplo, lapack_int n, lapack_int nrhs,
              const float * a, float * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_spptrs SPPTRS;
    
     /* LAPACKE SPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spptrf) ( int matrix_layout, 
           char uplo , lapack_int n , float * a );

    Fptr_NL_LAPACKE_spptrf SPPTRF;

    void *hModule, *dModule;
    spptrs_obj = new  pptrs_float_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    SPPTRS = (Fptr_NL_LAPACKE_spptrs)dlsym(hModule, "LAPACKE_spptrs");
    ASSERT_TRUE(SPPTRS != NULL) << "failed to get the Netlib LAPACKE_spptrs symbol";

    SPPTRF = (Fptr_NL_LAPACKE_spptrf)dlsym(hModule,"LAPACKE_spptrf");
    ASSERT_TRUE(SPPTRF != NULL) << "failed to get the Netlib LAPACKE_spptrf symbol";

    /* Pre condition: need to call pptrf - before calling pptrs function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    spptrs_obj->inforef = SPPTRF( spptrs_obj->matrix_layout, 
          spptrs_obj->uplo, spptrs_obj->n,
         spptrs_obj->aref); 

    spptrs_obj->inforef = SPPTRS( spptrs_obj->matrix_layout, 
          spptrs_obj->uplo, spptrs_obj->n,
         spptrs_obj->nrhs, (const float *)spptrs_obj->aref, 
                          spptrs_obj->bref, spptrs_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    
     spptrs_obj->info = LAPACKE_spptrf( spptrs_obj->matrix_layout, 
               spptrs_obj->uplo, spptrs_obj->n,
                 spptrs_obj->a);

    spptrs_obj->info = LAPACKE_spptrs( spptrs_obj->matrix_layout, 
               spptrs_obj->uplo, spptrs_obj->n,
                 spptrs_obj->nrhs, (const float *)spptrs_obj->a, 
                                 spptrs_obj->b, spptrs_obj->ldb );

    if( spptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spptrs is wrong\n", 
                    spptrs_obj->info );
    }
    if( spptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spptrs is wrong\n", 
        spptrs_obj->inforef );
    }
    
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(spptrs_test, spptrs1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spptrs_obj->b_bufsize,
                           spptrs_obj->b, spptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(spptrs_test, spptrs2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spptrs_obj->b_bufsize,
                           spptrs_obj->b, spptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spptrs_test, spptrs3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spptrs_obj->b_bufsize,
                           spptrs_obj->b, spptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(spptrs_test, spptrs4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spptrs_obj->b_bufsize,
                           spptrs_obj->b, spptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pptrs_scomplex_parameters  class definition */
class pptrs_scomplex_parameters{
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
      pptrs_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                    lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~pptrs_scomplex_parameters (); 
};  /* end of pptrs_scomplex_parameters  class definition */


/* Constructor pptrs_scomplex_parameters definition */
pptrs_scomplex_parameters:: pptrs_scomplex_parameters ( int matrix_layout_i, 
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
       EXPECT_FALSE( true) << "pptrs_scomplex_parameters object: malloc error. Exiting ";
       pptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

pptrs_scomplex_parameters:: ~pptrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptrs_scomplex_parameters object: destructor invoked. \n");
#endif
//   pptrs_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class cpptrs_test  : public  ::testing::Test {
public:
   pptrs_scomplex_parameters  *cpptrs_obj;
   void SetUp();  
   void TearDown () { delete cpptrs_obj; }
};

void cpptrs_test::SetUp(){

    /* LAPACKE CPPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_cpptrs) ( int matrix_layout, 
         char uplo, lapack_int n, lapack_int nrhs,
              const lapack_complex_float * a, lapack_complex_float * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_cpptrs CPPTRS;
    
     /* LAPACKE CPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpptrf) ( int matrix_layout, 
           char uplo , lapack_int n , lapack_complex_float * a );

    Fptr_NL_LAPACKE_cpptrf CPPTRF;

    void *hModule, *dModule;
    cpptrs_obj = new  pptrs_scomplex_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    CPPTRS = (Fptr_NL_LAPACKE_cpptrs)dlsym(hModule, "LAPACKE_cpptrs");
    ASSERT_TRUE(CPPTRS != NULL) << "failed to get the Netlib LAPACKE_cpptrs symbol";

    CPPTRF = (Fptr_NL_LAPACKE_cpptrf)dlsym(hModule,"LAPACKE_cpptrf");
    ASSERT_TRUE(CPPTRF != NULL) << "failed to get the Netlib LAPACKE_cpptrf symbol";

    /* Pre condition: need to call pptrf - before calling pptrs function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    cpptrs_obj->inforef = CPPTRF( cpptrs_obj->matrix_layout, 
          cpptrs_obj->uplo, cpptrs_obj->n,
         cpptrs_obj->aref); 

    cpptrs_obj->inforef = CPPTRS( cpptrs_obj->matrix_layout, 
          cpptrs_obj->uplo, cpptrs_obj->n,
         cpptrs_obj->nrhs, (const lapack_complex_float *)cpptrs_obj->aref, 
                          cpptrs_obj->bref, cpptrs_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    
     cpptrs_obj->info = LAPACKE_cpptrf( cpptrs_obj->matrix_layout, 
               cpptrs_obj->uplo, cpptrs_obj->n,
                 cpptrs_obj->a);

    cpptrs_obj->info = LAPACKE_cpptrs( cpptrs_obj->matrix_layout, 
               cpptrs_obj->uplo, cpptrs_obj->n,
                 cpptrs_obj->nrhs, (const lapack_complex_float *)cpptrs_obj->a, 
                                 cpptrs_obj->b, cpptrs_obj->ldb );

    if( cpptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpptrs is wrong\n", 
                    cpptrs_obj->info );
    }
    if( cpptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpptrs is wrong\n", 
        cpptrs_obj->inforef );
    }
    
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(cpptrs_test, cpptrs1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpptrs_obj->b_bufsize,
                           cpptrs_obj->b, cpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpptrs_test, cpptrs2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpptrs_obj->b_bufsize,
                           cpptrs_obj->b, cpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpptrs_test, cpptrs3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpptrs_obj->b_bufsize,
                           cpptrs_obj->b, cpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cpptrs_test, cpptrs4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpptrs_obj->b_bufsize,
                           cpptrs_obj->b, cpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pptrs_dcomplex_parameters  class definition */
class pptrs_dcomplex_parameters{
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
      pptrs_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                    lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~pptrs_dcomplex_parameters (); 
};  /* end of pptrs_dcomplex_parameters  class definition */


/* Constructor pptrs_dcomplex_parameters definition */
pptrs_dcomplex_parameters:: pptrs_dcomplex_parameters ( int matrix_layout_i, 
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
       EXPECT_FALSE( true) << "pptrs_dcomplex_parameters object: malloc error. Exiting ";
       pptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

pptrs_dcomplex_parameters:: ~pptrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptrs_dcomplex_parameters object: destructor invoked. \n");
#endif
//   pptrs_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class zpptrs_test  : public  ::testing::Test {
public:
   pptrs_dcomplex_parameters  *zpptrs_obj;
   void SetUp();  
   void TearDown () { delete zpptrs_obj; }
};

void zpptrs_test::SetUp(){

    /* LAPACKE ZPPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zpptrs) ( int matrix_layout, 
         char uplo, lapack_int n, lapack_int nrhs,
              const lapack_complex_double * a, lapack_complex_double * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_zpptrs ZPPTRS;
    
     /* LAPACKE ZPPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpptrf) ( int matrix_layout, 
           char uplo , lapack_int n , lapack_complex_double * a );

    Fptr_NL_LAPACKE_zpptrf ZPPTRF;

    void *hModule, *dModule;
    zpptrs_obj = new  pptrs_dcomplex_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    ZPPTRS = (Fptr_NL_LAPACKE_zpptrs)dlsym(hModule, "LAPACKE_zpptrs");
    ASSERT_TRUE(ZPPTRS != NULL) << "failed to get the Netlib LAPACKE_zpptrs symbol";

    ZPPTRF = (Fptr_NL_LAPACKE_zpptrf)dlsym(hModule,"LAPACKE_zpptrf");
    ASSERT_TRUE(ZPPTRF != NULL) << "failed to get the Netlib LAPACKE_zpptrf symbol";

    /* Pre condition: need to call pptrf - before calling pptrs function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    zpptrs_obj->inforef = ZPPTRF( zpptrs_obj->matrix_layout, 
          zpptrs_obj->uplo, zpptrs_obj->n,
         zpptrs_obj->aref); 

    zpptrs_obj->inforef = ZPPTRS( zpptrs_obj->matrix_layout, 
          zpptrs_obj->uplo, zpptrs_obj->n,
         zpptrs_obj->nrhs, (const lapack_complex_double *)zpptrs_obj->aref, 
                          zpptrs_obj->bref, zpptrs_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    
     zpptrs_obj->info = LAPACKE_zpptrf( zpptrs_obj->matrix_layout, 
               zpptrs_obj->uplo, zpptrs_obj->n,
                 zpptrs_obj->a);

    zpptrs_obj->info = LAPACKE_zpptrs( zpptrs_obj->matrix_layout, 
               zpptrs_obj->uplo, zpptrs_obj->n,
                 zpptrs_obj->nrhs, (const lapack_complex_double *)zpptrs_obj->a, 
                                 zpptrs_obj->b, zpptrs_obj->ldb );

    if( zpptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpptrs is wrong\n", 
                    zpptrs_obj->info );
    }
    if( zpptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpptrs is wrong\n", 
        zpptrs_obj->inforef );
    }
    
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(zpptrs_test, zpptrs1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpptrs_obj->b_bufsize,
                           zpptrs_obj->b, zpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zpptrs_test, zpptrs2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpptrs_obj->b_bufsize,
                           zpptrs_obj->b, zpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpptrs_test, zpptrs3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpptrs_obj->b_bufsize,
                           zpptrs_obj->b, zpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zpptrs_test, zpptrs4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpptrs_obj->b_bufsize,
                           zpptrs_obj->b, zpptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}