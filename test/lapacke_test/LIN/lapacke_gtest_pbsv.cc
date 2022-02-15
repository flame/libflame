#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define pbsv_free() \
       if (b != NULL)    free (b   ); \
       if (bref != NULL) free (bref); \
       if (ab != NULL)    free (ab  ); \
       if (abref != NULL) free (abref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin pbsv_double_parameters  class definition */
class pbsv_double_parameters{
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
      pbsv_double_parameters ( int matrix_layout_i, char uplo,
                             lapack_int n_i, lapack_int ldab_i, 
                            lapack_int kd_i, lapack_int nrhs_i,
                                             lapack_int ldb_i);
             
      ~pbsv_double_parameters ();
};  /* end of pbsv_double_parameters  class definition */


/* Constructor pbsv_double_parameters definition */
pbsv_double_parameters:: pbsv_double_parameters( int matrix_layout_i,
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
   printf(" \n pbsv Double:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "pbsv_double_parameters object: malloc error.";
       pbsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

pbsv_double_parameters:: ~pbsv_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbsv_double_parameters object: destructor invoked. \n");
#endif
   pbsv_free();
}


//  Test fixture class definition
class dpbsv_test  : public  ::testing::Test {
public:
   pbsv_double_parameters  *dpbsv_obj;
   void SetUp();  
   void TearDown () { delete dpbsv_obj; }
};


void dpbsv_test::SetUp(){

    /* LAPACKE DPBSV prototype */
    typedef int (*Fptr_NL_LAPACKE_dpbsv) ( int matrix_layout, 
                       char uplo, lapack_int n, lapack_int kd, 
                           lapack_int nrhs,  double * ab,
                lapack_int ldab,double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dpbsv DPBSV;
    dpbsv_obj = new  pbsv_double_parameters(
                  lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].ldab,
                           lin_solver_paramslist[idx].kd,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    dpbsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpbsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpbsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpbsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DPBSV = (Fptr_NL_LAPACKE_dpbsv)dlsym(dpbsv_obj->hModule, "LAPACKE_dpbsv");
    ASSERT_TRUE(DPBSV != NULL) << "failed to get the Netlib LAPACKE_dpbsv symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    dpbsv_obj->inforef = DPBSV( dpbsv_obj->matrix_layout, dpbsv_obj->uplo,
                             dpbsv_obj->n,dpbsv_obj->kd, dpbsv_obj->nrhs,
                           ( double *)dpbsv_obj->abref,dpbsv_obj->ldab,
                                          dpbsv_obj->bref, dpbsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dpbsv_obj->info = LAPACKE_dpbsv( dpbsv_obj->matrix_layout,
                                dpbsv_obj->uplo, dpbsv_obj->n, 
                               dpbsv_obj->kd,  dpbsv_obj->nrhs,
                                  ( double *)dpbsv_obj->ab,
               dpbsv_obj->ldab, dpbsv_obj->b, dpbsv_obj->ldb );


    if( dpbsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
               LAPACKE_dpbsv is wrong\n",    dpbsv_obj->info );
    }
    if( dpbsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpbsv is wrong\n",
        dpbsv_obj->inforef );
    }
}

TEST_F(dpbsv_test, dpbsv1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpbsv_obj->b_bufsize, dpbsv_obj->b,
                                            dpbsv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbsv_test, dpbsv2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpbsv_obj->b_bufsize, dpbsv_obj->b,
                                            dpbsv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbsv_test, dpbsv3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpbsv_obj->b_bufsize, dpbsv_obj->b,
                                            dpbsv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbsv_test, dpbsv4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpbsv_obj->b_bufsize, dpbsv_obj->b,
                                            dpbsv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pbsv_float_parameters  class definition */
class pbsv_float_parameters{

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
      pbsv_float_parameters (int matrix_layout_i, char uplo, lapack_int n_i,
                                          lapack_int ldab_i, lapack_int kd_i,
                                        lapack_int nrhs_i, lapack_int ldb_i);
      ~pbsv_float_parameters ();
};  /* end of pbsv_float_parameters  class definition */


/* Constructor pbsv_float_parameters definition */
pbsv_float_parameters:: pbsv_float_parameters ( int matrix_layout_i,
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
   printf(" \n pbsv float:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "pbsv_double_parameters object: malloc error.";
       pbsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( ab, abref, ldab*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

pbsv_float_parameters:: ~pbsv_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbsv_float_parameters object: destructor invoked. \n");
#endif
   pbsv_free();
}


//  Test fixture class definition
class spbsv_test  : public  ::testing::Test {
public:
   pbsv_float_parameters  *spbsv_obj;
   void SetUp();  
   void TearDown () { delete spbsv_obj; }
};


void spbsv_test::SetUp(){

    /* LAPACKE SPBSV prototype */
    typedef int (*Fptr_NL_LAPACKE_spbsv) ( int matrix_layout, char uplo,
                            lapack_int n, lapack_int kd, lapack_int nrhs, 
                                         float *ab, lapack_int ldab,
                                              float *b, lapack_int ldb );

    Fptr_NL_LAPACKE_spbsv SPBSV;
    spbsv_obj = new  pbsv_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kd,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    spbsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spbsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spbsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spbsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SPBSV = (Fptr_NL_LAPACKE_spbsv)dlsym(spbsv_obj->hModule, "LAPACKE_spbsv");
    ASSERT_TRUE(SPBSV != NULL) << "failed to get the Netlib LAPACKE_spbsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    spbsv_obj->inforef = SPBSV( spbsv_obj->matrix_layout,spbsv_obj->uplo,
                             spbsv_obj->n, spbsv_obj->kd, spbsv_obj->nrhs,
                             ( float *)spbsv_obj->abref,
                         spbsv_obj->ldab, spbsv_obj->bref, spbsv_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    spbsv_obj->info = LAPACKE_spbsv( spbsv_obj->matrix_layout,
                                spbsv_obj->uplo, spbsv_obj->n, 
                              spbsv_obj->kd, spbsv_obj->nrhs,
                                  ( float *)spbsv_obj->ab,
                    spbsv_obj->ldab, spbsv_obj->b, spbsv_obj->ldb );
                    
    if( spbsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spbsv is wrong\n",
                    spbsv_obj->info );
    }
    if( spbsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spbsv is wrong\n",
        spbsv_obj->inforef );
    }
}

TEST_F(spbsv_test, spbsv1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spbsv_obj->b_bufsize, spbsv_obj->b,
                            spbsv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbsv_test, spbsv2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spbsv_obj->b_bufsize, spbsv_obj->b,
                            spbsv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbsv_test, spbsv3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spbsv_obj->b_bufsize, spbsv_obj->b,
                            spbsv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbsv_test, spbsv4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spbsv_obj->b_bufsize, spbsv_obj->b,
                            spbsv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin pbsv_scomplex_parameters  class definition */
class pbsv_scomplex_parameters{
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
      pbsv_scomplex_parameters ( int matrix_layout_i, char uplo,
                                lapack_int n_i, lapack_int ldab_i,
                                 lapack_int kd_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~pbsv_scomplex_parameters ();
};  /* end of pbsv_scomplex_parameters  class definition */


/* Constructor pbsv_scomplex_parameters definition */
pbsv_scomplex_parameters:: pbsv_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n pbsv scomplex:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "pbsv_scomplex_parameters object: malloc error.";
       pbsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

pbsv_scomplex_parameters:: ~pbsv_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbsv_scomplex_parameters object: destructor invoked. \n");
#endif
   pbsv_free();
}


//  Test fixture class definition
class cpbsv_test  : public  ::testing::Test {
public:
   pbsv_scomplex_parameters  *cpbsv_obj;
   void SetUp();  
   void TearDown () { delete cpbsv_obj; }
};


void cpbsv_test::SetUp(){

    /* LAPACKE CPBSV prototype */
    typedef int (*Fptr_NL_LAPACKE_cpbsv) ( int matrix_layout, char uplo,
                                             lapack_int n, lapack_int kd,  
                         lapack_int nrhs,  lapack_complex_float * ab,
                                                          lapack_int ldab,
                              lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_cpbsv CPBSV;

    cpbsv_obj = new  pbsv_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kd,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    cpbsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpbsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpbsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpbsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CPBSV = (Fptr_NL_LAPACKE_cpbsv)dlsym(cpbsv_obj->hModule, "LAPACKE_cpbsv");
    ASSERT_TRUE(CPBSV != NULL) << "failed to get the Netlib LAPACKE_cpbsv symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    cpbsv_obj->inforef = CPBSV( cpbsv_obj->matrix_layout,
                                  cpbsv_obj->uplo, cpbsv_obj->n,
                                  cpbsv_obj->kd, cpbsv_obj->nrhs,
                   ( lapack_complex_float *)cpbsv_obj->abref,
              cpbsv_obj->ldab, cpbsv_obj->bref, cpbsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    cpbsv_obj->info = LAPACKE_cpbsv( cpbsv_obj->matrix_layout,
                cpbsv_obj->uplo, cpbsv_obj->n, cpbsv_obj->kd,
   cpbsv_obj->nrhs, ( lapack_complex_float *)cpbsv_obj->ab,
               cpbsv_obj->ldab, cpbsv_obj->b, cpbsv_obj->ldb );

    if( cpbsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cpbsv is wrong\n", cpbsv_obj->info );
    }
    if( cpbsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpbsv is wrong\n",
        cpbsv_obj->inforef );
    }
}

TEST_F(cpbsv_test, cpbsv1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpbsv_obj->b_bufsize,
                           cpbsv_obj->b, cpbsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbsv_test, cpbsv2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpbsv_obj->b_bufsize,
                            cpbsv_obj->b, cpbsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbsv_test, cpbsv3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpbsv_obj->b_bufsize,
                           cpbsv_obj->b, cpbsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbsv_test, cpbsv4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpbsv_obj->b_bufsize,
                           cpbsv_obj->b, cpbsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin pbsv_dcomplex_parameters  class definition */
class pbsv_dcomplex_parameters{
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
      pbsv_dcomplex_parameters ( int matrix_layout_i, char uplo,
                                lapack_int n_i, lapack_int ldab_i,
                                 lapack_int kd_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~pbsv_dcomplex_parameters ();
};  /* end of pbsv_dcomplex_parameters  class definition */


/* Constructor pbsv_dcomplex_parameters definition */
pbsv_dcomplex_parameters:: pbsv_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n pbsv DComplex:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "pbsv_dcomplex_parameters object: malloc error.";
       pbsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( ab, abref, n*ldab);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

pbsv_dcomplex_parameters:: ~pbsv_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbsv_dcomplex_parameters object: destructor invoked. \n");
#endif
   pbsv_free();
}


//  Test fixture class definition
class zpbsv_test  : public  ::testing::Test {
public:
   pbsv_dcomplex_parameters  *zpbsv_obj;
   void SetUp();  
   void TearDown () { delete zpbsv_obj; }
};


void zpbsv_test::SetUp(){

    /* LAPACKE ZPBSV prototype */
    typedef int (*Fptr_NL_LAPACKE_zpbsv) ( int matrix_layout, char uplo,
                                            lapack_int n, lapack_int kd, 
                        lapack_int nrhs,  lapack_complex_double * ab,
                                                        lapack_int ldab,
                              lapack_complex_double * b, lapack_int ldb );

    Fptr_NL_LAPACKE_zpbsv ZPBSV;

    zpbsv_obj = new  pbsv_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].ldab,
                                         lin_solver_paramslist[idx].kd,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zpbsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpbsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpbsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpbsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZPBSV = (Fptr_NL_LAPACKE_zpbsv)dlsym(zpbsv_obj->hModule, "LAPACKE_zpbsv");
    ASSERT_TRUE(ZPBSV != NULL) << "failed to get the Netlib LAPACKE_zpbsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    zpbsv_obj->inforef = ZPBSV( zpbsv_obj->matrix_layout, zpbsv_obj->uplo, 
                        zpbsv_obj->n, zpbsv_obj->kd, zpbsv_obj->nrhs,
                        ( lapack_complex_double *)zpbsv_obj->abref,
                        zpbsv_obj->ldab, zpbsv_obj->bref, zpbsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zpbsv_obj->info = LAPACKE_zpbsv( zpbsv_obj->matrix_layout,
                zpbsv_obj->uplo, zpbsv_obj->n, zpbsv_obj->kd,
                zpbsv_obj->nrhs, 
                ( lapack_complex_double *)zpbsv_obj->ab,
                zpbsv_obj->ldab, zpbsv_obj->b, zpbsv_obj->ldb );


    if( zpbsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpbsv is wrong\n",
                    zpbsv_obj->info );
    }
    if( zpbsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpbsv is wrong\n",
        zpbsv_obj->inforef );
    }
}

TEST_F(zpbsv_test, zpbsv1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpbsv_obj->b_bufsize,
                           zpbsv_obj->b, zpbsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbsv_test, zpbsv2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpbsv_obj->b_bufsize,
                           zpbsv_obj->b, zpbsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbsv_test, zpbsv3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpbsv_obj->b_bufsize,
                           zpbsv_obj->b, zpbsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbsv_test, zpbsv4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpbsv_obj->b_bufsize,
                           zpbsv_obj->b, zpbsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
