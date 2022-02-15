#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define pbtrs_free() \
       if (b != NULL)    free (b   ); \
       if (bref != NULL) free (bref); \
       if (ab != NULL)    free (ab  ); \
       if (abref != NULL) free (abref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin pbtrs_double_parameters  class definition */
class pbtrs_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; // storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.      
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double *ab, *abref; //The array 'a' contains the matrix A
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      pbtrs_double_parameters ( int matrix_layout_i, char uplo,
                             lapack_int n_i, lapack_int ldab_i, 
                            lapack_int kd_i, lapack_int nrhs_i,
                                             lapack_int ldb_i);
             
      ~pbtrs_double_parameters ();
};  /* end of pbtrs_double_parameters  class definition */


/* Constructor pbtrs_double_parameters definition */
pbtrs_double_parameters:: pbtrs_double_parameters( int matrix_layout_i,
                        char uplo_i, lapack_int n_i, lapack_int ldab_i,
                lapack_int kd_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    kd = kd_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pbtrs Double:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, ldab, ldb, nrhs);
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

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       pbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

pbtrs_double_parameters:: ~pbtrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbtrs_double_parameters object: destructor invoked. \n");
#endif
   pbtrs_free();
}


//  Test fixture class definition
class dpbtrs_test  : public  ::testing::Test {
public:
   pbtrs_double_parameters  *dpbtrs_obj;
   void SetUp();  
   void TearDown () { delete dpbtrs_obj; }
};


void dpbtrs_test::SetUp(){

    /* LAPACKE DPBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dpbtrs) ( int matrix_layout, 
                       char uplo, lapack_int n, lapack_int kd, 
                           lapack_int nrhs, const double * ab,
                lapack_int ldab,double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dpbtrs DPBTRS;

     /* LAPACKE DPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpbtrf) ( int matrix_layout,
                       char uplo, lapack_int n, lapack_int kd,
                                 double* ab,lapack_int ldab );

    Fptr_NL_LAPACKE_dpbtrf DPBTRF;

    dpbtrs_obj = new  pbtrs_double_parameters(
                  lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].ldab,
                           lin_solver_paramslist[idx].kd,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    dpbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DPBTRS = (Fptr_NL_LAPACKE_dpbtrs)dlsym(dpbtrs_obj->hModule, "LAPACKE_dpbtrs");
    ASSERT_TRUE(DPBTRS != NULL) << "failed to get the Netlib LAPACKE_dpbtrs symbol";

    DPBTRF = (Fptr_NL_LAPACKE_dpbtrf)dlsym(dpbtrs_obj->hModule,"LAPACKE_dpbtrf");
    ASSERT_TRUE(DPBTRF != NULL) << "failed to get the Netlib LAPACKE_dpbtrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    dpbtrs_obj->inforef = DPBTRF( dpbtrs_obj->matrix_layout, dpbtrs_obj->uplo,
                             dpbtrs_obj->n, dpbtrs_obj->kd, dpbtrs_obj->abref,
                                                           dpbtrs_obj->ldab );

    dpbtrs_obj->inforef = DPBTRS( dpbtrs_obj->matrix_layout, dpbtrs_obj->uplo,
                             dpbtrs_obj->n,dpbtrs_obj->kd, dpbtrs_obj->nrhs,
                           (const double *)dpbtrs_obj->abref,dpbtrs_obj->ldab,
                                          dpbtrs_obj->bref, dpbtrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dpbtrs_obj->info = LAPACKE_dpbtrf( dpbtrs_obj->matrix_layout,
                                dpbtrs_obj->uplo, dpbtrs_obj->n,
                                 dpbtrs_obj->kd, dpbtrs_obj->ab,
                                              dpbtrs_obj->ldab);

    dpbtrs_obj->info = LAPACKE_dpbtrs( dpbtrs_obj->matrix_layout,
                                dpbtrs_obj->uplo, dpbtrs_obj->n, 
                               dpbtrs_obj->kd,  dpbtrs_obj->nrhs,
                                  (const double *)dpbtrs_obj->ab,
               dpbtrs_obj->ldab, dpbtrs_obj->b, dpbtrs_obj->ldb );


    if( dpbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
               LAPACKE_dpbtrs is wrong\n",    dpbtrs_obj->info );
    }
    if( dpbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpbtrs is wrong\n",
        dpbtrs_obj->inforef );
    }
}

TEST_F(dpbtrs_test, dpbtrs1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpbtrs_obj->b_bufsize, dpbtrs_obj->b,
                                            dpbtrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbtrs_test, dpbtrs2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpbtrs_obj->b_bufsize, dpbtrs_obj->b,
                                            dpbtrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbtrs_test, dpbtrs3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpbtrs_obj->b_bufsize, dpbtrs_obj->b,
                                            dpbtrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbtrs_test, dpbtrs4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpbtrs_obj->b_bufsize, dpbtrs_obj->b,
                                            dpbtrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pbtrs_float_parameters  class definition */
class pbtrs_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.      
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float *ab, *abref; //The array 'a' contains the matrix A
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      pbtrs_float_parameters (int matrix_layout_i, char uplo, lapack_int n_i,
                                          lapack_int ldab_i, lapack_int kd_i,
                                        lapack_int nrhs_i, lapack_int ldb_i);
      ~pbtrs_float_parameters ();
};  /* end of pbtrs_float_parameters  class definition */


/* Constructor pbtrs_float_parameters definition */
pbtrs_float_parameters:: pbtrs_float_parameters ( int matrix_layout_i,
                      char uplo_i, lapack_int n_i, lapack_int ldab_i,
                lapack_int kd_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    kd = kd_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
    uplo = uplo_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pbtrs float:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, ldab, ldb, nrhs);
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

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       pbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( ab, abref, ldab*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

pbtrs_float_parameters:: ~pbtrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbtrs_float_parameters object: destructor invoked. \n");
#endif
   pbtrs_free();
}


//  Test fixture class definition
class spbtrs_test  : public  ::testing::Test {
public:
   pbtrs_float_parameters  *spbtrs_obj;
   void SetUp();  
   void TearDown () { delete spbtrs_obj; }
};


void spbtrs_test::SetUp(){

    /* LAPACKE SPBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_spbtrs) ( int matrix_layout, char uplo,
                            lapack_int n, lapack_int kd, lapack_int nrhs, 
                                        const float *ab, lapack_int ldab,
                                              float *b, lapack_int ldb );

    Fptr_NL_LAPACKE_spbtrs SPBTRS;

     /* LAPACKE SPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spbtrf) ( int matrix_layout,char uplo,
                  lapack_int n, lapack_int kd, float *ab,lapack_int ldab);

    Fptr_NL_LAPACKE_spbtrf SPBTRF;

    spbtrs_obj = new  pbtrs_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kd,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    spbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SPBTRS = (Fptr_NL_LAPACKE_spbtrs)dlsym(spbtrs_obj->hModule, "LAPACKE_spbtrs");
    ASSERT_TRUE(SPBTRS != NULL) << "failed to get the Netlib LAPACKE_spbtrs symbol";

    SPBTRF = (Fptr_NL_LAPACKE_spbtrf)dlsym(spbtrs_obj->hModule,"LAPACKE_spbtrf");
    ASSERT_TRUE(SPBTRF != NULL) << "failed to get the Netlib LAPACKE_spbtrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    spbtrs_obj->inforef = SPBTRF( spbtrs_obj->matrix_layout, spbtrs_obj->uplo,
                          spbtrs_obj->n, spbtrs_obj->kd, spbtrs_obj->abref,
                                                          spbtrs_obj->ldab );

    spbtrs_obj->inforef = SPBTRS( spbtrs_obj->matrix_layout,spbtrs_obj->uplo,
                             spbtrs_obj->n, spbtrs_obj->kd, spbtrs_obj->nrhs,
                             (const float *)spbtrs_obj->abref,
                         spbtrs_obj->ldab, spbtrs_obj->bref, spbtrs_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    spbtrs_obj->info = LAPACKE_spbtrf( spbtrs_obj->matrix_layout,
                        spbtrs_obj->uplo, spbtrs_obj->n,spbtrs_obj->kd,
                                     spbtrs_obj->ab, spbtrs_obj->ldab);

    spbtrs_obj->info = LAPACKE_spbtrs( spbtrs_obj->matrix_layout,
                                spbtrs_obj->uplo, spbtrs_obj->n, 
                              spbtrs_obj->kd, spbtrs_obj->nrhs,
                                  (const float *)spbtrs_obj->ab,
                    spbtrs_obj->ldab, spbtrs_obj->b, spbtrs_obj->ldb );
                    
    if( spbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spbtrs is wrong\n",
                    spbtrs_obj->info );
    }
    if( spbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spbtrs is wrong\n",
        spbtrs_obj->inforef );
    }
}

TEST_F(spbtrs_test, spbtrs1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spbtrs_obj->b_bufsize, spbtrs_obj->b,
                            spbtrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbtrs_test, spbtrs2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spbtrs_obj->b_bufsize, spbtrs_obj->b,
                            spbtrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbtrs_test, spbtrs3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spbtrs_obj->b_bufsize, spbtrs_obj->b,
                            spbtrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbtrs_test, spbtrs4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spbtrs_obj->b_bufsize, spbtrs_obj->b,
                            spbtrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin pbtrs_scomplex_parameters  class definition */
class pbtrs_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.      
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *ab, *abref; //The array 'a' contains the matrix A
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      pbtrs_scomplex_parameters ( int matrix_layout_i, char uplo,
                                lapack_int n_i, lapack_int ldab_i,
                                 lapack_int kd_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~pbtrs_scomplex_parameters ();
};  /* end of pbtrs_scomplex_parameters  class definition */


/* Constructor pbtrs_scomplex_parameters definition */
pbtrs_scomplex_parameters:: pbtrs_scomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, lapack_int ldab_i,
                                       lapack_int kd_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    kd = kd_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pbtrs scomplex:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, ldab, ldb, nrhs);
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

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       pbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

pbtrs_scomplex_parameters:: ~pbtrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbtrs_scomplex_parameters object: destructor invoked. \n");
#endif
   pbtrs_free();
}


//  Test fixture class definition
class cpbtrs_test  : public  ::testing::Test {
public:
   pbtrs_scomplex_parameters  *cpbtrs_obj;
   void SetUp();  
   void TearDown () { delete cpbtrs_obj; }
};


void cpbtrs_test::SetUp(){

    /* LAPACKE CPBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_cpbtrs) ( int matrix_layout, char uplo,
                                             lapack_int n, lapack_int kd,  
                         lapack_int nrhs, const lapack_complex_float * ab,
                                                          lapack_int ldab,
                              lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_cpbtrs CPBTRS;

     /* LAPACKE CPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpbtrf) ( int matrix_layout,char uplo,
                                            lapack_int n, lapack_int kd, 
                              lapack_complex_float* ab,lapack_int ldab);

    Fptr_NL_LAPACKE_cpbtrf CPBTRF;


    cpbtrs_obj = new  pbtrs_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kd,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    cpbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CPBTRS = (Fptr_NL_LAPACKE_cpbtrs)dlsym(cpbtrs_obj->hModule, "LAPACKE_cpbtrs");
    ASSERT_TRUE(CPBTRS != NULL) << "failed to get the Netlib LAPACKE_cpbtrs symbol";

    CPBTRF = (Fptr_NL_LAPACKE_cpbtrf)dlsym(cpbtrs_obj->hModule,"LAPACKE_cpbtrf");
    ASSERT_TRUE(CPBTRF != NULL) << "failed to get the Netlib LAPACKE_cpbtrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    cpbtrs_obj->inforef = CPBTRF( cpbtrs_obj->matrix_layout,
                            cpbtrs_obj->uplo, cpbtrs_obj->n,
                          cpbtrs_obj->kd, cpbtrs_obj->abref,
                                         cpbtrs_obj->ldab );

    cpbtrs_obj->inforef = CPBTRS( cpbtrs_obj->matrix_layout,
                                  cpbtrs_obj->uplo, cpbtrs_obj->n,
                                  cpbtrs_obj->kd, cpbtrs_obj->nrhs,
                   (const lapack_complex_float *)cpbtrs_obj->abref,
              cpbtrs_obj->ldab, cpbtrs_obj->bref, cpbtrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    cpbtrs_obj->info = LAPACKE_cpbtrf( cpbtrs_obj->matrix_layout,
                                    cpbtrs_obj->uplo, cpbtrs_obj->n,
                                    cpbtrs_obj->kd, cpbtrs_obj->ab,
                                                  cpbtrs_obj->ldab);

    cpbtrs_obj->info = LAPACKE_cpbtrs( cpbtrs_obj->matrix_layout,
                cpbtrs_obj->uplo, cpbtrs_obj->n, cpbtrs_obj->kd,
   cpbtrs_obj->nrhs, (const lapack_complex_float *)cpbtrs_obj->ab,
               cpbtrs_obj->ldab, cpbtrs_obj->b, cpbtrs_obj->ldb );

    if( cpbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cpbtrs is wrong\n", cpbtrs_obj->info );
    }
    if( cpbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpbtrs is wrong\n",
        cpbtrs_obj->inforef );
    }
}

TEST_F(cpbtrs_test, cpbtrs1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpbtrs_obj->b_bufsize,
                           cpbtrs_obj->b, cpbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbtrs_test, cpbtrs2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpbtrs_obj->b_bufsize,
                            cpbtrs_obj->b, cpbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbtrs_test, cpbtrs3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpbtrs_obj->b_bufsize,
                           cpbtrs_obj->b, cpbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbtrs_test, cpbtrs4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpbtrs_obj->b_bufsize,
                           cpbtrs_obj->b, cpbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin pbtrs_dcomplex_parameters  class definition */
class pbtrs_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.      
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int ldab;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *ab, *abref; //The array 'a' contains the matrix A
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      pbtrs_dcomplex_parameters ( int matrix_layout_i, char uplo,
                                lapack_int n_i, lapack_int ldab_i,
                                 lapack_int kd_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~pbtrs_dcomplex_parameters ();
};  /* end of pbtrs_dcomplex_parameters  class definition */


/* Constructor pbtrs_dcomplex_parameters definition */
pbtrs_dcomplex_parameters:: pbtrs_dcomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, lapack_int ldab_i,
                                                        lapack_int kd_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    kd = kd_i;
    ldab = ldab_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pbtrs DComplex:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, ldab, ldb, nrhs);
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

    if( (ab==NULL) || (abref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       pbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

pbtrs_dcomplex_parameters:: ~pbtrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbtrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   pbtrs_free();
}


//  Test fixture class definition
class zpbtrs_test  : public  ::testing::Test {
public:
   pbtrs_dcomplex_parameters  *zpbtrs_obj;
   void SetUp();  
   void TearDown () { delete zpbtrs_obj; }
};


void zpbtrs_test::SetUp(){

    /* LAPACKE ZPBTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zpbtrs) ( int matrix_layout, char uplo,
                                            lapack_int n, lapack_int kd, 
                        lapack_int nrhs, const lapack_complex_double * ab,
                                                        lapack_int ldab,
                              lapack_complex_double * b, lapack_int ldb );

    Fptr_NL_LAPACKE_zpbtrs ZPBTRS;

     /* LAPACKE ZPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpbtrf) ( int matrix_layout,char uplo,
                                            lapack_int n, lapack_int kd, 
                            lapack_complex_double* ab,lapack_int ldab);

    Fptr_NL_LAPACKE_zpbtrf ZPBTRF;


    zpbtrs_obj = new  pbtrs_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kd,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zpbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZPBTRS = (Fptr_NL_LAPACKE_zpbtrs)dlsym(zpbtrs_obj->hModule, "LAPACKE_zpbtrs");
    ASSERT_TRUE(ZPBTRS != NULL) << "failed to get the Netlib LAPACKE_zpbtrs symbol";

    ZPBTRF = (Fptr_NL_LAPACKE_zpbtrf)dlsym(zpbtrs_obj->hModule,"LAPACKE_zpbtrf");
    ASSERT_TRUE(ZPBTRF != NULL) << "failed to get the Netlib LAPACKE_zpbtrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zpbtrs_obj->inforef = ZPBTRF( zpbtrs_obj->matrix_layout, zpbtrs_obj->uplo,
         zpbtrs_obj->n, zpbtrs_obj->kd, zpbtrs_obj->abref, zpbtrs_obj->ldab );

    zpbtrs_obj->inforef = ZPBTRS( zpbtrs_obj->matrix_layout, zpbtrs_obj->uplo, 
                        zpbtrs_obj->n, zpbtrs_obj->kd, zpbtrs_obj->nrhs,
                        (const lapack_complex_double *)zpbtrs_obj->abref,
                        zpbtrs_obj->ldab, zpbtrs_obj->bref, zpbtrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zpbtrs_obj->info = LAPACKE_zpbtrf( zpbtrs_obj->matrix_layout,
                zpbtrs_obj->uplo, zpbtrs_obj->n, zpbtrs_obj->kd,
                               zpbtrs_obj->ab, zpbtrs_obj->ldab);

    zpbtrs_obj->info = LAPACKE_zpbtrs( zpbtrs_obj->matrix_layout,
                zpbtrs_obj->uplo, zpbtrs_obj->n, zpbtrs_obj->kd,
                zpbtrs_obj->nrhs, 
                (const lapack_complex_double *)zpbtrs_obj->ab,
                zpbtrs_obj->ldab, zpbtrs_obj->b, zpbtrs_obj->ldb );


    if( zpbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpbtrs is wrong\n",
                    zpbtrs_obj->info );
    }
    if( zpbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpbtrs is wrong\n",
        zpbtrs_obj->inforef );
    }
}

TEST_F(zpbtrs_test, zpbtrs1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpbtrs_obj->b_bufsize,
                           zpbtrs_obj->b, zpbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbtrs_test, zpbtrs2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpbtrs_obj->b_bufsize,
                           zpbtrs_obj->b, zpbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbtrs_test, zpbtrs3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpbtrs_obj->b_bufsize,
                           zpbtrs_obj->b, zpbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbtrs_test, zpbtrs4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpbtrs_obj->b_bufsize,
                           zpbtrs_obj->b, zpbtrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
