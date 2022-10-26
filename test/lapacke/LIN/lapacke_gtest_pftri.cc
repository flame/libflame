#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define pftri_free() \
       if (a != NULL)    free (a  ); \
       if (aref != NULL) free (aref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin pftri_float_parameters  class definition */
class pftri_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      float *a,*aref; //The array a contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pftri_float_parameters ( int matrix_layout_i, char uplo_i, 
                                    char trans_i, lapack_int n_i);
              
      ~pftri_float_parameters (); 
};  /* end of pftri_float_parameters  class definition */


/* Constructor pftri_float_parameters definition */
pftri_float_parameters:: pftri_float_parameters ( int matrix_layout_i, 
                              char trans_i, char uplo_i, lapack_int n_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
    diff = 0;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*(n + 1)/2) );

    if( (a==NULL) || (aref==NULL) ){
       EXPECT_FALSE( true) << "pftri_float_parameters object: malloc error. Exiting ";
       pftri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*(n + 1)/2);
    
   } /* end of Constructor  */

pftri_float_parameters:: ~pftri_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftri_float_parameters object: destructor invoked. \n");
#endif
   pftri_free();
}


//  Test fixture class definition
class spftri_test  : public  ::testing::Test {
public:
   pftri_float_parameters  *spftri_obj;
   void SetUp();  
   void TearDown () { delete spftri_obj; }
};

void spftri_test::SetUp(){

    /* LAPACKE SPFTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_spftri) ( int matrix_layout, 
        char transr, char uplo, lapack_int n, float * a  );

    Fptr_NL_LAPACKE_spftri SPFTRI;
    
     /* LAPACKE SPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spftrf) ( int matrix_layout, 
          char transr , char uplo , lapack_int n , float * a );

    Fptr_NL_LAPACKE_spftrf SPFTRF;

    spftri_obj = new  pftri_float_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    spftri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spftri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spftri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spftri_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SPFTRI = (Fptr_NL_LAPACKE_spftri)dlsym(spftri_obj->hModule, "LAPACKE_spftri");
    ASSERT_TRUE(SPFTRI != NULL) << "failed to get the Netlib LAPACKE_spftri symbol";

    SPFTRF = (Fptr_NL_LAPACKE_spftrf)dlsym(spftri_obj->hModule,"LAPACKE_spftrf");
    ASSERT_TRUE(SPFTRF != NULL) << "failed to get the Netlib LAPACKE_spftrf symbol";

    /* Pre condition: need to call pftrf - before calling pftri function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    spftri_obj->inforef = SPFTRF(   spftri_obj->matrix_layout, 
                                    spftri_obj->trans,
                                    spftri_obj->uplo,
                                    spftri_obj->n,
                                    spftri_obj->aref);

    spftri_obj->inforef = SPFTRI(   spftri_obj->matrix_layout, 
                                    spftri_obj->trans,
                                    spftri_obj->uplo,
                                    spftri_obj->n,
                                    spftri_obj->aref );

    /* Compute libflame's Lapacke o/p  */
    
    spftri_obj->info = LAPACKE_spftrf(  spftri_obj->matrix_layout, 
                                        spftri_obj->trans,
                                        spftri_obj->uplo,
                                        spftri_obj->n,
                                        spftri_obj->a); 

    spftri_obj->info = LAPACKE_spftri(  spftri_obj->matrix_layout, 
                                        spftri_obj->trans,
                                        spftri_obj->uplo,
                                        spftri_obj->n,
                                        spftri_obj->a );

    if( spftri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spftri is wrong\n", 
                    spftri_obj->info );
    }
    if( spftri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spftri is wrong\n", 
        spftri_obj->inforef );
    }
    
    /* Compute Difference between libflame and Netlib o/ps  */
    spftri_obj->diff =  computeDiff_s(  spftri_obj->n*(spftri_obj->n+1)/2,
                                        spftri_obj->a,
                                        spftri_obj->aref );
}


TEST_F(spftri_test, spftri1) {
    EXPECT_NEAR(0.0, spftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spftri_test, spftri2) {
    EXPECT_NEAR(0.0, spftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spftri_test, spftri3) {
    EXPECT_NEAR(0.0, spftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spftri_test, spftri4) {
    EXPECT_NEAR(0.0, spftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pftri_double_parameters  class definition */
class pftri_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      double *a,*aref; //The array a contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pftri_double_parameters ( int matrix_layout_i, char uplo_i, 
                                    char trans_i, lapack_int n_i);
              
      ~pftri_double_parameters (); 
};  /* end of pftri_double_parameters  class definition */


/* Constructor pftri_double_parameters definition */
pftri_double_parameters:: pftri_double_parameters ( int matrix_layout_i, 
                              char trans_i, char uplo_i, lapack_int n_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
    diff = 0;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*(n + 1)/2) );

    if( (a==NULL) || (aref==NULL) ){
       EXPECT_FALSE( true) << "pftri_double_parameters object: malloc error. Exiting ";
       pftri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*(n + 1)/2);
    
   } /* end of Constructor  */

pftri_double_parameters:: ~pftri_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftri_double_parameters object: destructor invoked. \n");
#endif
   pftri_free();
}


//  Test fixture class definition
class dpftri_test  : public  ::testing::Test {
public:
   pftri_double_parameters  *dpftri_obj;
   void SetUp();  
   void TearDown () { delete dpftri_obj; }
};

void dpftri_test::SetUp(){

    /* LAPACKE DPFTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_dpftri) ( int matrix_layout, 
        char transr, char uplo, lapack_int n, double * a  );

    Fptr_NL_LAPACKE_dpftri DPFTRI;
    
     /* LAPACKE DPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpftrf) ( int matrix_layout, 
          char transr , char uplo , lapack_int n , double * a );

    Fptr_NL_LAPACKE_dpftrf DPFTRF;

    dpftri_obj = new  pftri_double_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    dpftri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpftri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpftri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpftri_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DPFTRI = (Fptr_NL_LAPACKE_dpftri)dlsym(dpftri_obj->hModule, "LAPACKE_dpftri");
    ASSERT_TRUE(DPFTRI != NULL) << "failed to get the Netlib LAPACKE_dpftri symbol";

    DPFTRF = (Fptr_NL_LAPACKE_dpftrf)dlsym(dpftri_obj->hModule,"LAPACKE_dpftrf");
    ASSERT_TRUE(DPFTRF != NULL) << "failed to get the Netlib LAPACKE_dpftrf symbol";

    /* Pre condition: need to call pftrf - before calling pftri function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API 
    dpftri_obj->inforef = DPFTRF(   dpftri_obj->matrix_layout, 
                                    dpftri_obj->trans,
                                    dpftri_obj->uplo,
                                    dpftri_obj->n,
                                    dpftri_obj->aref);

    dpftri_obj->inforef = DPFTRI(   dpftri_obj->matrix_layout, 
                                    dpftri_obj->trans,
                                    dpftri_obj->uplo,
                                    dpftri_obj->n,
                                    dpftri_obj->aref );*/ 

    /* Compute libflame's Lapacke o/p  */
    
    dpftri_obj->info = LAPACKE_dpftrf(  dpftri_obj->matrix_layout, 
                                        dpftri_obj->trans,
                                        dpftri_obj->uplo,
                                        dpftri_obj->n,
                                        dpftri_obj->a); 

    dpftri_obj->info = LAPACKE_dpftri(  dpftri_obj->matrix_layout, 
                                        dpftri_obj->trans,
                                        dpftri_obj->uplo,
                                        dpftri_obj->n,
                                        dpftri_obj->a );

    if( dpftri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpftri is wrong\n", 
                    dpftri_obj->info );
    }
    if( dpftri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpftri is wrong\n", 
        dpftri_obj->inforef );
    }
    
    /* Compute Difference between libflame and Netlib o/ps  */
    dpftri_obj->diff =  computeDiff_d(  dpftri_obj->n*(dpftri_obj->n+1)/2,
                                        dpftri_obj->a,
                                        dpftri_obj->aref );
}


TEST_F(dpftri_test, dpftri1) {
    EXPECT_NEAR(0.0, dpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpftri_test, dpftri2) {
    EXPECT_NEAR(0.0, dpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpftri_test, dpftri3) {
    EXPECT_NEAR(0.0, dpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpftri_test, dpftri4) {
    EXPECT_NEAR(0.0, dpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pftri_scomplex_parameters  class definition */
class pftri_scomplex_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_float *a,*aref; //The array a contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pftri_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                    char trans_i, lapack_int n_i);
              
      ~pftri_scomplex_parameters (); 
};  /* end of pftri_scomplex_parameters  class definition */


/* Constructor pftri_scomplex_parameters definition */
pftri_scomplex_parameters:: pftri_scomplex_parameters ( int matrix_layout_i, 
                              char trans_i, char uplo_i, lapack_int n_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
    diff = 0;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*(n + 1)/2) );

    if( (a==NULL) || (aref==NULL) ){
       EXPECT_FALSE( true) << "pftri_scomplex_parameters object: malloc error. Exiting ";
       pftri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*(n + 1)/2);
    
   } /* end of Constructor  */

pftri_scomplex_parameters:: ~pftri_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftri_scomplex_parameters object: destructor invoked. \n");
#endif
   pftri_free();
}


//  Test fixture class definition
class cpftri_test  : public  ::testing::Test {
public:
   pftri_scomplex_parameters  *cpftri_obj;
   void SetUp();  
   void TearDown () { delete cpftri_obj; }
};

void cpftri_test::SetUp(){

    /* LAPACKE CPFTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_cpftri) ( int matrix_layout, 
        char transr, char uplo, lapack_int n, lapack_complex_float * a  );

    Fptr_NL_LAPACKE_cpftri CPFTRI;
    
     /* LAPACKE CPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpftrf) ( int matrix_layout, 
          char transr , char uplo , lapack_int n , lapack_complex_float * a );

    Fptr_NL_LAPACKE_cpftrf CPFTRF;

    cpftri_obj = new  pftri_scomplex_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    cpftri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpftri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpftri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpftri_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CPFTRI = (Fptr_NL_LAPACKE_cpftri)dlsym(cpftri_obj->hModule, "LAPACKE_cpftri");
    ASSERT_TRUE(CPFTRI != NULL) << "failed to get the Netlib LAPACKE_cpftri symbol";

    CPFTRF = (Fptr_NL_LAPACKE_cpftrf)dlsym(cpftri_obj->hModule,"LAPACKE_cpftrf");
    ASSERT_TRUE(CPFTRF != NULL) << "failed to get the Netlib LAPACKE_cpftrf symbol";

    /* Pre condition: need to call pftrf - before calling pftri function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
//#if 0
    cpftri_obj->inforef = CPFTRF(   cpftri_obj->matrix_layout, 
                                    cpftri_obj->trans,
                                    cpftri_obj->uplo,
                                    cpftri_obj->n,
                                    cpftri_obj->aref);

    cpftri_obj->inforef = CPFTRI(   cpftri_obj->matrix_layout, 
                                    cpftri_obj->trans,
                                    cpftri_obj->uplo,
                                    cpftri_obj->n,
                                    cpftri_obj->aref );
//#endif
    /* Compute libflame's Lapacke o/p  */
    
    cpftri_obj->info = LAPACKE_cpftrf(  cpftri_obj->matrix_layout, 
                                        cpftri_obj->trans,
                                        cpftri_obj->uplo,
                                        cpftri_obj->n,
                                        cpftri_obj->a); 

    cpftri_obj->info = LAPACKE_cpftri(  cpftri_obj->matrix_layout, 
                                        cpftri_obj->trans,
                                        cpftri_obj->uplo,
                                        cpftri_obj->n,
                                        cpftri_obj->a );

    if( cpftri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpftri is wrong\n", 
                    cpftri_obj->info );
    }
    if( cpftri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpftri is wrong\n", 
        cpftri_obj->inforef );
    }
    
    /* Compute Difference between libflame and Netlib o/ps  */
    int abuf_size = (cpftri_obj->n) * (cpftri_obj->n)/2;
    cpftri_obj->diff =  computeDiff_c(  abuf_size,
                                        cpftri_obj->a,
                                        cpftri_obj->aref );
}


TEST_F(cpftri_test, cpftri1) {
    EXPECT_NEAR(0.0, cpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpftri_test, cpftri2) {
    EXPECT_NEAR(0.0, cpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpftri_test, cpftri3) {
    EXPECT_NEAR(0.0, cpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpftri_test, cpftri4) {
    EXPECT_NEAR(0.0, cpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pftri_dcomplex_parameters  class definition */
class pftri_dcomplex_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_double *a,*aref; //The array a contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pftri_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                    char trans_i, lapack_int n_i);
              
      ~pftri_dcomplex_parameters (); 
};  /* end of pftri_dcomplex_parameters  class definition */


/* Constructor pftri_dcomplex_parameters definition */
pftri_dcomplex_parameters:: pftri_dcomplex_parameters ( int matrix_layout_i, 
                              char trans_i, char uplo_i, lapack_int n_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*(n + 1)/2) );

    if( (a==NULL) || (aref==NULL) ){
       EXPECT_FALSE( true) << "pftri_dcomplex_parameters object: malloc error. Exiting ";
       pftri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*(n + 1)/2);
    
   } /* end of Constructor  */

pftri_dcomplex_parameters:: ~pftri_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftri_dcomplex_parameters object: destructor invoked. \n");
#endif
   pftri_free();
}


//  Test fixture class definition
class zpftri_test  : public  ::testing::Test {
public:
   pftri_dcomplex_parameters  *zpftri_obj;
   void SetUp();  
   void TearDown () { delete zpftri_obj; }
};

void zpftri_test::SetUp(){

    /* LAPACKE ZPFTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_zpftri) ( int matrix_layout, 
        char transr, char uplo, lapack_int n, lapack_complex_double * a  );

    Fptr_NL_LAPACKE_zpftri ZPFTRI;
    
     /* LAPACKE ZPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpftrf) ( int matrix_layout, 
          char transr , char uplo , lapack_int n , lapack_complex_double * a );

    Fptr_NL_LAPACKE_zpftrf ZPFTRF;

    zpftri_obj = new  pftri_dcomplex_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    zpftri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpftri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpftri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpftri_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZPFTRI = (Fptr_NL_LAPACKE_zpftri)dlsym(zpftri_obj->hModule, "LAPACKE_zpftri");
    ASSERT_TRUE(ZPFTRI != NULL) << "failed to get the Netlib LAPACKE_zpftri symbol";

    ZPFTRF = (Fptr_NL_LAPACKE_zpftrf)dlsym(zpftri_obj->hModule,"LAPACKE_zpftrf");
    ASSERT_TRUE(ZPFTRF != NULL) << "failed to get the Netlib LAPACKE_zpftrf symbol";

    /* Pre condition: need to call pftrf - before calling pftri function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zpftri_obj->inforef = ZPFTRF(   zpftri_obj->matrix_layout, 
                                    zpftri_obj->trans,
                                    zpftri_obj->uplo,
                                    zpftri_obj->n,
                                    zpftri_obj->aref);

    zpftri_obj->inforef = ZPFTRI(   zpftri_obj->matrix_layout, 
                                    zpftri_obj->trans,
                                    zpftri_obj->uplo,
                                    zpftri_obj->n,
                                    zpftri_obj->aref );

    /* Compute libflame's Lapacke o/p  */
    
    zpftri_obj->info = LAPACKE_zpftrf(  zpftri_obj->matrix_layout, 
                                        zpftri_obj->trans,
                                        zpftri_obj->uplo,
                                        zpftri_obj->n,
                                        zpftri_obj->a); 

    zpftri_obj->info = LAPACKE_zpftri(  zpftri_obj->matrix_layout, 
                                        zpftri_obj->trans,
                                        zpftri_obj->uplo,
                                        zpftri_obj->n,
                                        zpftri_obj->a );

    if( zpftri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpftri is wrong\n", 
                    zpftri_obj->info );
    }
    if( zpftri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpftri is wrong\n", 
        zpftri_obj->inforef );
    }
    
    /* Compute Difference between libflame and Netlib o/ps  */
    zpftri_obj->diff =  computeDiff_z(  zpftri_obj->n*(zpftri_obj->n+1)/2,
                                        zpftri_obj->a,
                                        zpftri_obj->aref );
}


TEST_F(zpftri_test, zpftri1) {
    EXPECT_NEAR(0.0, zpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpftri_test, zpftri2) {
    EXPECT_NEAR(0.0, zpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpftri_test, zpftri3) {
    EXPECT_NEAR(0.0, zpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpftri_test, zpftri4) {
    EXPECT_NEAR(0.0, zpftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
