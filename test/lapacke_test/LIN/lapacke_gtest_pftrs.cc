#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define pftrs_free() \
       if (b != NULL)    free (b   ); \
       if (bref != NULL) free (bref); \
       if (a != NULL)    free (a  ); \
       if (aref != NULL) free (aref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin pftrs_double_parameters  class definition */
class pftrs_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
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
      pftrs_double_parameters ( int matrix_layout_i, char uplo_i, 
                                    char trans_i, lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~pftrs_double_parameters (); 
};  /* end of pftrs_double_parameters  class definition */


/* Constructor pftrs_double_parameters definition */
pftrs_double_parameters:: pftrs_double_parameters ( int matrix_layout_i, 
                              char trans_i, char uplo_i, lapack_int n_i,
                 lapack_int nrhs_i, lapack_int ldb_i) {
    int buf_siz;                                
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
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
       EXPECT_FALSE( true) << "pftrs_double_parameters object: malloc error. Exiting ";
       pftrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

pftrs_double_parameters:: ~pftrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftrs_double_parameters object: destructor invoked. \n");
#endif
//   pftrs_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class dpftrs_test  : public  ::testing::Test {
public:
   pftrs_double_parameters  *dpftrs_obj;
   void SetUp();  
   void TearDown () { delete dpftrs_obj; }
};

void dpftrs_test::SetUp(){

    /* LAPACKE DPFTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dpftrs) ( int matrix_layout, 
	    char transr, char uplo, lapack_int n, lapack_int nrhs,
	          const double * a, double * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_dpftrs DPFTRS;
	
     /* LAPACKE DPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpftrf) ( int matrix_layout, 
          char transr , char uplo , lapack_int n , double * a );

    Fptr_NL_LAPACKE_dpftrf DPFTRF;

    void *hModule, *dModule;
    dpftrs_obj = new  pftrs_double_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    DPFTRS = (Fptr_NL_LAPACKE_dpftrs)dlsym(hModule, "LAPACKE_dpftrs");
    ASSERT_TRUE(DPFTRS != NULL) << "failed to get the Netlib LAPACKE_dpftrs symbol";

    DPFTRF = (Fptr_NL_LAPACKE_dpftrf)dlsym(hModule,"LAPACKE_dpftrf");
    ASSERT_TRUE(DPFTRF != NULL) << "failed to get the Netlib LAPACKE_dpftrf symbol";

    /* Pre condition: need to call pftrf - before calling pftrs function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    dpftrs_obj->inforef = DPFTRF( dpftrs_obj->matrix_layout, 
	      dpftrs_obj->trans, dpftrs_obj->uplo, dpftrs_obj->n,
		 dpftrs_obj->aref);

    dpftrs_obj->inforef = DPFTRS( dpftrs_obj->matrix_layout, 
	      dpftrs_obj->trans, dpftrs_obj->uplo, dpftrs_obj->n,
		 dpftrs_obj->nrhs, (const double *)dpftrs_obj->aref, 
                          dpftrs_obj->bref, dpftrs_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
	
    dpftrs_obj->info = LAPACKE_dpftrf( dpftrs_obj->matrix_layout, 
               dpftrs_obj->trans, dpftrs_obj->uplo, dpftrs_obj->n,
				 dpftrs_obj->a); 

    dpftrs_obj->info = LAPACKE_dpftrs( dpftrs_obj->matrix_layout, 
               dpftrs_obj->trans, dpftrs_obj->uplo, dpftrs_obj->n,
				 dpftrs_obj->nrhs, (const double *)dpftrs_obj->a, 
                                 dpftrs_obj->b, dpftrs_obj->ldb );

    if( dpftrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpftrs is wrong\n", 
                    dpftrs_obj->info );
    }
    if( dpftrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpftrs is wrong\n", 
        dpftrs_obj->inforef );
    }
	
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(dpftrs_test, dpftrs1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpftrs_obj->b_bufsize,
                           dpftrs_obj->b, dpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dpftrs_test, dpftrs2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpftrs_obj->b_bufsize,
                           dpftrs_obj->b, dpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpftrs_test, dpftrs3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpftrs_obj->b_bufsize,
                           dpftrs_obj->b, dpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dpftrs_test, dpftrs4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpftrs_obj->b_bufsize,
                           dpftrs_obj->b, dpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pftrs_float_parameters  class definition */
class pftrs_float_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
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
      pftrs_float_parameters ( int matrix_layout_i, char uplo_i, 
                                    char trans_i, lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~pftrs_float_parameters (); 
};  /* end of pftrs_float_parameters  class definition */


/* Constructor pftrs_float_parameters definition */
pftrs_float_parameters:: pftrs_float_parameters ( int matrix_layout_i, 
                              char trans_i, char uplo_i, lapack_int n_i,
                 lapack_int nrhs_i, lapack_int ldb_i) {
    int buf_siz;                                
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
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
       EXPECT_FALSE( true) << "pftrs_float_parameters object: malloc error. Exiting ";
       pftrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

pftrs_float_parameters:: ~pftrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftrs_float_parameters object: destructor invoked. \n");
#endif
//   pftrs_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class spftrs_test  : public  ::testing::Test {
public:
   pftrs_float_parameters  *spftrs_obj;
   void SetUp();  
   void TearDown () { delete spftrs_obj; }
};

void spftrs_test::SetUp(){

    /* LAPACKE SPFTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_spftrs) ( int matrix_layout, 
	    char transr, char uplo, lapack_int n, lapack_int nrhs,
	          const float * a, float * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_spftrs SPFTRS;
	
     /* LAPACKE SPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spftrf) ( int matrix_layout, 
          char transr , char uplo , lapack_int n , float * a );

    Fptr_NL_LAPACKE_spftrf SPFTRF;

    void *hModule, *dModule;
    spftrs_obj = new  pftrs_float_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    SPFTRS = (Fptr_NL_LAPACKE_spftrs)dlsym(hModule, "LAPACKE_spftrs");
    ASSERT_TRUE(SPFTRS != NULL) << "failed to get the Netlib LAPACKE_spftrs symbol";

    SPFTRF = (Fptr_NL_LAPACKE_spftrf)dlsym(hModule,"LAPACKE_spftrf");
    ASSERT_TRUE(SPFTRF != NULL) << "failed to get the Netlib LAPACKE_spftrf symbol";

    /* Pre condition: need to call pftrf - before calling pftrs function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    spftrs_obj->inforef = SPFTRF( spftrs_obj->matrix_layout, 
	      spftrs_obj->trans, spftrs_obj->uplo, spftrs_obj->n,
		 spftrs_obj->aref); 

    spftrs_obj->inforef = SPFTRS( spftrs_obj->matrix_layout, 
	      spftrs_obj->trans, spftrs_obj->uplo, spftrs_obj->n,
		 spftrs_obj->nrhs, (const float *)spftrs_obj->aref, 
                          spftrs_obj->bref, spftrs_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
	
     spftrs_obj->info = LAPACKE_spftrf( spftrs_obj->matrix_layout, 
               spftrs_obj->trans, spftrs_obj->uplo, spftrs_obj->n,
				 spftrs_obj->a);

    spftrs_obj->info = LAPACKE_spftrs( spftrs_obj->matrix_layout, 
               spftrs_obj->trans, spftrs_obj->uplo, spftrs_obj->n,
				 spftrs_obj->nrhs, (const float *)spftrs_obj->a, 
                                 spftrs_obj->b, spftrs_obj->ldb );

    if( spftrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spftrs is wrong\n", 
                    spftrs_obj->info );
    }
    if( spftrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spftrs is wrong\n", 
        spftrs_obj->inforef );
    }
	
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(spftrs_test, spftrs1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spftrs_obj->b_bufsize,
                           spftrs_obj->b, spftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(spftrs_test, spftrs2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spftrs_obj->b_bufsize,
                           spftrs_obj->b, spftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spftrs_test, spftrs3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spftrs_obj->b_bufsize,
                           spftrs_obj->b, spftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(spftrs_test, spftrs4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spftrs_obj->b_bufsize,
                           spftrs_obj->b, spftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pftrs_scomplex_parameters  class definition */
class pftrs_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
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
      pftrs_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                    char trans_i, lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~pftrs_scomplex_parameters (); 
};  /* end of pftrs_scomplex_parameters  class definition */


/* Constructor pftrs_scomplex_parameters definition */
pftrs_scomplex_parameters:: pftrs_scomplex_parameters ( int matrix_layout_i, 
                              char trans_i, char uplo_i, lapack_int n_i,
                 lapack_int nrhs_i, lapack_int ldb_i) {
    int buf_siz;                                
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
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
       EXPECT_FALSE( true) << "pftrs_scomplex_parameters object: malloc error. Exiting ";
       pftrs_free();
	   exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

pftrs_scomplex_parameters:: ~pftrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftrs_scomplex_parameters object: destructor invoked. \n");
#endif
//   pftrs_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class cpftrs_test  : public  ::testing::Test {
public:
   pftrs_scomplex_parameters  *cpftrs_obj;
   void SetUp();  
   void TearDown () { delete cpftrs_obj; }
};

void cpftrs_test::SetUp(){

    /* LAPACKE CPFTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_cpftrs) ( int matrix_layout, 
	    char transr, char uplo, lapack_int n, lapack_int nrhs,
	          const lapack_complex_float * a, lapack_complex_float * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_cpftrs CPFTRS;
	
     /* LAPACKE CPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpftrf) ( int matrix_layout, 
          char transr , char uplo , lapack_int n , lapack_complex_float * a );

    Fptr_NL_LAPACKE_cpftrf CPFTRF;

    void *hModule, *dModule;
    cpftrs_obj = new  pftrs_scomplex_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    CPFTRS = (Fptr_NL_LAPACKE_cpftrs)dlsym(hModule, "LAPACKE_cpftrs");
    ASSERT_TRUE(CPFTRS != NULL) << "failed to get the Netlib LAPACKE_cpftrs symbol";

    CPFTRF = (Fptr_NL_LAPACKE_cpftrf)dlsym(hModule,"LAPACKE_cpftrf");
    ASSERT_TRUE(CPFTRF != NULL) << "failed to get the Netlib LAPACKE_cpftrf symbol";

    /* Pre condition: need to call pftrf - before calling pftrs function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    cpftrs_obj->inforef = CPFTRF( cpftrs_obj->matrix_layout, 
	      cpftrs_obj->trans, cpftrs_obj->uplo, cpftrs_obj->n,
		 cpftrs_obj->aref); 

    cpftrs_obj->inforef = CPFTRS( cpftrs_obj->matrix_layout, 
	      cpftrs_obj->trans, cpftrs_obj->uplo, cpftrs_obj->n,
		 cpftrs_obj->nrhs, (const lapack_complex_float *)cpftrs_obj->aref, 
                          cpftrs_obj->bref, cpftrs_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
	
     cpftrs_obj->info = LAPACKE_cpftrf( cpftrs_obj->matrix_layout, 
               cpftrs_obj->trans, cpftrs_obj->uplo, cpftrs_obj->n,
				 cpftrs_obj->a);

    cpftrs_obj->info = LAPACKE_cpftrs( cpftrs_obj->matrix_layout, 
               cpftrs_obj->trans, cpftrs_obj->uplo, cpftrs_obj->n,
				 cpftrs_obj->nrhs, (const lapack_complex_float *)cpftrs_obj->a, 
                                 cpftrs_obj->b, cpftrs_obj->ldb );

    if( cpftrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpftrs is wrong\n", 
                    cpftrs_obj->info );
    }
    if( cpftrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpftrs is wrong\n", 
        cpftrs_obj->inforef );
    }
	
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(cpftrs_test, cpftrs1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpftrs_obj->b_bufsize,
                           cpftrs_obj->b, cpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
#if 0
TEST_F(cpftrs_test, cpftrs2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpftrs_obj->b_bufsize,
                           cpftrs_obj->b, cpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpftrs_test, cpftrs3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpftrs_obj->b_bufsize,
                           cpftrs_obj->b, cpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(cpftrs_test, cpftrs4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpftrs_obj->b_bufsize,
                           cpftrs_obj->b, cpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
#endif

/* Begin pftrs_dcomplex_parameters  class definition */
class pftrs_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; // 'N' or 'T' or 'C'. Indicates the form of the equations:
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
      pftrs_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                    char trans_i, lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~pftrs_dcomplex_parameters (); 
};  /* end of pftrs_dcomplex_parameters  class definition */


/* Constructor pftrs_dcomplex_parameters definition */
pftrs_dcomplex_parameters:: pftrs_dcomplex_parameters ( int matrix_layout_i, 
                              char trans_i, char uplo_i, lapack_int n_i,
                 lapack_int nrhs_i, lapack_int ldb_i) {
    int buf_siz;                                
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
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
       EXPECT_FALSE( true) << "pftrs_dcomplex_parameters object: malloc error. Exiting ";
       pftrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*(n + 1)/2);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

pftrs_dcomplex_parameters:: ~pftrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pftrs_dcomplex_parameters object: destructor invoked. \n");
#endif
//   pftrs_free();
       if (b != NULL)    free (b   );
       if (bref != NULL) free (bref);
       if (a != NULL)    free (a  ); 
       if (aref != NULL) free (aref);

}


//  Test fixture class definition
class zpftrs_test  : public  ::testing::Test {
public:
   pftrs_dcomplex_parameters  *zpftrs_obj;
   void SetUp();  
   void TearDown () { delete zpftrs_obj; }
};

void zpftrs_test::SetUp(){

    /* LAPACKE ZPFTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zpftrs) ( int matrix_layout, 
	    char transr, char uplo, lapack_int n, lapack_int nrhs,
	          const lapack_complex_double * a, lapack_complex_double * b, lapack_int ldb  );
    Fptr_NL_LAPACKE_zpftrs ZPFTRS;
	
     /* LAPACKE ZPFTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpftrf) ( int matrix_layout, 
          char transr , char uplo , lapack_int n , lapack_complex_double * a );

    Fptr_NL_LAPACKE_zpftrf ZPFTRF;

    void *hModule, *dModule;
    zpftrs_obj = new  pftrs_dcomplex_parameters( lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].transr,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);
    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";

    ZPFTRS = (Fptr_NL_LAPACKE_zpftrs)dlsym(hModule, "LAPACKE_zpftrs");
    ASSERT_TRUE(ZPFTRS != NULL) << "failed to get the Netlib LAPACKE_zpftrs symbol";

    ZPFTRF = (Fptr_NL_LAPACKE_zpftrf)dlsym(hModule,"LAPACKE_zpftrf");
    ASSERT_TRUE(ZPFTRF != NULL) << "failed to get the Netlib LAPACKE_zpftrf symbol";

    /* Pre condition: need to call pftrf - before calling pftrs function */
    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    zpftrs_obj->inforef = ZPFTRF( zpftrs_obj->matrix_layout, 
	      zpftrs_obj->trans, zpftrs_obj->uplo, zpftrs_obj->n,
		 zpftrs_obj->aref); 

    zpftrs_obj->inforef = ZPFTRS( zpftrs_obj->matrix_layout, 
	      zpftrs_obj->trans, zpftrs_obj->uplo, zpftrs_obj->n,
		 zpftrs_obj->nrhs, (const lapack_complex_double *)zpftrs_obj->aref, 
                          zpftrs_obj->bref, zpftrs_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
	
     zpftrs_obj->info = LAPACKE_zpftrf( zpftrs_obj->matrix_layout, 
               zpftrs_obj->trans, zpftrs_obj->uplo, zpftrs_obj->n,
				 zpftrs_obj->a);

    zpftrs_obj->info = LAPACKE_zpftrs( zpftrs_obj->matrix_layout, 
               zpftrs_obj->trans, zpftrs_obj->uplo, zpftrs_obj->n,
				 zpftrs_obj->nrhs, (const lapack_complex_double *)zpftrs_obj->a, 
                                 zpftrs_obj->b, zpftrs_obj->ldb );

    if( zpftrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpftrs is wrong\n", 
                    zpftrs_obj->info );
    }
    if( zpftrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpftrs is wrong\n", 
        zpftrs_obj->inforef );
    }
	
    if( hModule != NULL) dlclose(hModule);
    if( dModule != NULL) dlclose(dModule);
}

TEST_F(zpftrs_test, zpftrs1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpftrs_obj->b_bufsize,
                           zpftrs_obj->b, zpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zpftrs_test, zpftrs2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpftrs_obj->b_bufsize,
                           zpftrs_obj->b, zpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpftrs_test, zpftrs3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpftrs_obj->b_bufsize,
                           zpftrs_obj->b, zpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zpftrs_test, zpftrs4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpftrs_obj->b_bufsize,
                           zpftrs_obj->b, zpftrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}