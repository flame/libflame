#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define pbrfs_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (x != NULL)    free (x  ); \
  if (xref != NULL) free (xref); \
  if (af != NULL)    free (af  ); \
  if (afref != NULL) free (afref); \
  if (ferr != NULL)    free (ferr  ); \
  if (ferrref != NULL) free (ferrref); \
  if (berr != NULL)    free (berr  ); \
  if (berrref != NULL) free (berrref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin pbrfs_float_parameters  class definition */
class pbrfs_float_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr, diff_s;
      float diff_a, diff_af;
      
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.
      float *a, *aref; //The array ab contains the matrix A
      float *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pbrfs_float_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i, lapack_int kd_i);
              
      ~pbrfs_float_parameters (); 
};  /* end of pbrfs_float_parameters  class definition */


/* Constructor pbrfs_float_parameters definition */
pbrfs_float_parameters:: pbrfs_float_parameters ( int matrix_layout_i, 
                                          char uplo_i, lapack_int n_i,
                                lapack_int nrhs_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    n = n_i;
    nrhs = nrhs_i;
    kd = kd_i;
    lda = n;
    ldaf = n;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n pbrfs float:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, uplo, lda, 
ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       pbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand ( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(float)));
    memcpy(afref, a, (n*n*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(float)));
    memcpy(xref, b, (b_bufsize*sizeof(float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

pbrfs_float_parameters:: ~pbrfs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbrfs_float_parameters object: destructor invoked. \n");
#endif
   pbrfs_free();
}

//  Test fixture class definition
class spbrfs_test  : public  ::testing::Test {
public:
   pbrfs_float_parameters  *spbrfs_obj;
   void SetUp();  
   void TearDown () { delete spbrfs_obj; }
};


void spbrfs_test::SetUp(){

    /* LAPACKE SPBRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_spbrfs) (int matrix_layout, char uplo,
        lapack_int n, lapack_int kd, lapack_int nrhs, const float* ab,
        lapack_int ldab, const float* afb, lapack_int ldafb,
        const float* b, lapack_int ldb, float* x, lapack_int ldx,
        float* ferr, float* berr);

    Fptr_NL_LAPACKE_spbrfs SPBRFS;

     /* LAPACKE SPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spbtrf) ( int matrix_layout , char uplo ,
                lapack_int n , lapack_int kd, float * a , lapack_int lda );

    Fptr_NL_LAPACKE_spbtrf SPBTRF;

    typedef int (*Fptr_NL_LAPACKE_spbtrs) ( int matrix_layout , char uplo ,
        lapack_int n , lapack_int kd , lapack_int nrhs , const float *ab ,
        lapack_int ldab , float *b , lapack_int ldb );

    Fptr_NL_LAPACKE_spbtrs SPBTRS;


    spbrfs_obj = new pbrfs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd);

    //idx = Circular_Increment_Index(idx);

    spbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SPBTRF = (Fptr_NL_LAPACKE_spbtrf)dlsym(spbrfs_obj->hModule,"LAPACKE_spbtrf");
    ASSERT_TRUE(SPBTRF != NULL) << "failed to pbt the Netlib LAPACKE_spbtrf symbol";
    
    SPBTRS = (Fptr_NL_LAPACKE_spbtrs)dlsym(spbrfs_obj->hModule, "LAPACKE_spbtrs");
    ASSERT_TRUE(SPBTRS != NULL) << "failed to pbt the Netlib LAPACKE_spbtrs symbol";
    
    SPBRFS = (Fptr_NL_LAPACKE_spbrfs)dlsym(spbrfs_obj->hModule, "LAPACKE_spbrfs");
    ASSERT_TRUE(SPBRFS != NULL) << "failed to pbt the Netlib LAPACKE_spbrfs symbol";

    /* Preparing the i/ps for pbrfs by invoking PBTRF, PBTRS API calls */
    spbrfs_obj->inforef = SPBTRF( spbrfs_obj->matrix_layout,
                                  spbrfs_obj->uplo, spbrfs_obj->n,
                                  spbrfs_obj->kd,
                                  spbrfs_obj->afref,
                                  spbrfs_obj->lda);
                           
    spbrfs_obj->info = LAPACKE_spbtrf( spbrfs_obj->matrix_layout,
                                       spbrfs_obj->uplo, spbrfs_obj->n,
                                       spbrfs_obj->kd,
                                       spbrfs_obj->af,
                                       spbrfs_obj->lda);

    spbrfs_obj->inforef = SPBTRS( spbrfs_obj->matrix_layout,
                                  spbrfs_obj->uplo,
                                  spbrfs_obj->n,
                                  spbrfs_obj->kd,
                                  spbrfs_obj->nrhs,
                                  spbrfs_obj->afref,
                                  spbrfs_obj->lda,
                                  spbrfs_obj->xref,
                                  spbrfs_obj->ldb
                                  );

    spbrfs_obj->info = LAPACKE_spbtrs( spbrfs_obj->matrix_layout,
                                  spbrfs_obj->uplo,
                                  spbrfs_obj->n,
                                  spbrfs_obj->kd,
                                  spbrfs_obj->nrhs,
                                  spbrfs_obj->af,
                                  spbrfs_obj->lda,
                                  spbrfs_obj->x,
                                  spbrfs_obj->ldb
                                  );


    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    spbrfs_obj->inforef = SPBRFS( spbrfs_obj->matrix_layout,
                                  spbrfs_obj->uplo,
                                  spbrfs_obj->n,
                                  spbrfs_obj->kd,
                                  spbrfs_obj->nrhs,
                                  (const float*)spbrfs_obj->aref,
                                  spbrfs_obj->lda,
                                  spbrfs_obj->afref, spbrfs_obj->ldaf,
                                  spbrfs_obj->bref, spbrfs_obj->ldb,
                                  spbrfs_obj->xref, spbrfs_obj->ldx,
                                  spbrfs_obj->ferrref,
                                  spbrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    spbrfs_obj->info = LAPACKE_spbrfs( spbrfs_obj->matrix_layout,
                                  spbrfs_obj->uplo,
                                  spbrfs_obj->n,
                                  spbrfs_obj->kd,
                                  spbrfs_obj->nrhs,
                                  (const float*)spbrfs_obj->a,
                                  spbrfs_obj->lda,
                                  spbrfs_obj->af, spbrfs_obj->ldaf,
                                  spbrfs_obj->b, spbrfs_obj->ldb,
                                  spbrfs_obj->x, spbrfs_obj->ldx,
                                  spbrfs_obj->ferr,
                                  spbrfs_obj->berr);

    if( spbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_spbrfs is wrong\n", spbrfs_obj->info );
    }
    if( spbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spbrfs is wrong\n", 
        spbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    spbrfs_obj->diff_xerr =  computeDiff_s( spbrfs_obj->x_bufsize, 
                spbrfs_obj->x, spbrfs_obj->xref );

    spbrfs_obj->diff_berr =  computeDiff_s( spbrfs_obj->nrhs, 
                spbrfs_obj->berr, spbrfs_obj->berrref );
                
    spbrfs_obj->diff_ferr =  computeDiff_s( spbrfs_obj->nrhs, 
                spbrfs_obj->ferr, spbrfs_obj->ferrref );

}

TEST_F(spbrfs_test, spbrfs1) {
    EXPECT_NEAR(0.0, spbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, spbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, spbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);

}

TEST_F(spbrfs_test, spbrfs2) {
    EXPECT_NEAR(0.0, spbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, spbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, spbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(spbrfs_test, spbrfs3) {
    EXPECT_NEAR(0.0, spbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, spbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, spbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(spbrfs_test, spbrfs4) {
    EXPECT_NEAR(0.0, spbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, spbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, spbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin pbrfs_double_parameters  class definition */
class pbrfs_double_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr, diff_s;
      double diff_a, diff_af;
      
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.
      double *a, *aref; //The array ab contains the matrix A
      double *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pbrfs_double_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i, lapack_int kd_i);
              
      ~pbrfs_double_parameters (); 
};  /* end of pbrfs_double_parameters  class definition */


/* Constructor pbrfs_double_parameters definition */
pbrfs_double_parameters:: pbrfs_double_parameters ( int matrix_layout_i, 
                                          char uplo_i, lapack_int n_i,
                                lapack_int nrhs_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    n = n_i;
    nrhs = nrhs_i;
    kd = kd_i;
    lda = n;
    ldaf = n;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n pbrfs double:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, uplo, lda, 
ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       pbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand ( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(double)));
    memcpy(afref, a, (n*n*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(double)));
    memcpy(xref, b, (b_bufsize*sizeof(double)));
    
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

pbrfs_double_parameters:: ~pbrfs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbrfs_double_parameters object: destructor invoked. \n");
#endif
   pbrfs_free();
}

//  Test fixture class definition
class dpbrfs_test  : public  ::testing::Test {
public:
   pbrfs_double_parameters  *dpbrfs_obj;
   void SetUp();  
   void TearDown () { delete dpbrfs_obj; }
};


void dpbrfs_test::SetUp(){

    /* LAPACKE DPBRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_dpbrfs) (int matrix_layout, char uplo,
        lapack_int n, lapack_int kd, lapack_int nrhs, const double* ab,
        lapack_int ldab, const double* afb, lapack_int ldafb,
        const double* b, lapack_int ldb, double* x, lapack_int ldx,
        double* ferr, double* berr);

    Fptr_NL_LAPACKE_dpbrfs DPBRFS;

     /* LAPACKE DPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpbtrf) ( int matrix_layout , char uplo ,
                lapack_int n , lapack_int kd, double * a , lapack_int lda );

    Fptr_NL_LAPACKE_dpbtrf DPBTRF;

    typedef int (*Fptr_NL_LAPACKE_dpbtrs) ( int matrix_layout , char uplo ,
        lapack_int n , lapack_int kd , lapack_int nrhs , const double *ab ,
        lapack_int ldab , double *b , lapack_int ldb );

    Fptr_NL_LAPACKE_dpbtrs DPBTRS;


    dpbrfs_obj = new pbrfs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd);

    //idx = Circular_Increment_Index(idx);

    dpbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DPBTRF = (Fptr_NL_LAPACKE_dpbtrf)dlsym(dpbrfs_obj->hModule,"LAPACKE_dpbtrf");
    ASSERT_TRUE(DPBTRF != NULL) << "failed to pbt the Netlib LAPACKE_dpbtrf symbol";
    
    DPBTRS = (Fptr_NL_LAPACKE_dpbtrs)dlsym(dpbrfs_obj->hModule, "LAPACKE_dpbtrs");
    ASSERT_TRUE(DPBTRS != NULL) << "failed to pbt the Netlib LAPACKE_dpbtrs symbol";
    
    DPBRFS = (Fptr_NL_LAPACKE_dpbrfs)dlsym(dpbrfs_obj->hModule, "LAPACKE_dpbrfs");
    ASSERT_TRUE(DPBRFS != NULL) << "failed to pbt the Netlib LAPACKE_dpbrfs symbol";

    /* Preparing the i/ps for pbrfs by invoking PBTRF, PBTRS API calls */
    dpbrfs_obj->inforef = DPBTRF( dpbrfs_obj->matrix_layout,
                                  dpbrfs_obj->uplo, dpbrfs_obj->n,
                                  dpbrfs_obj->kd,
                                  dpbrfs_obj->afref,
                                  dpbrfs_obj->lda);
                           
    dpbrfs_obj->info = LAPACKE_dpbtrf( dpbrfs_obj->matrix_layout,
                                       dpbrfs_obj->uplo, dpbrfs_obj->n,
                                       dpbrfs_obj->kd,
                                       dpbrfs_obj->af,
                                       dpbrfs_obj->lda);

    dpbrfs_obj->inforef = DPBTRS( dpbrfs_obj->matrix_layout,
                                  dpbrfs_obj->uplo,
                                  dpbrfs_obj->n,
                                  dpbrfs_obj->kd,
                                  dpbrfs_obj->nrhs,
                                  dpbrfs_obj->afref,
                                  dpbrfs_obj->lda,
                                  dpbrfs_obj->xref,
                                  dpbrfs_obj->ldb
                                  );

    dpbrfs_obj->info = LAPACKE_dpbtrs( dpbrfs_obj->matrix_layout,
                                  dpbrfs_obj->uplo,
                                  dpbrfs_obj->n,
                                  dpbrfs_obj->kd,
                                  dpbrfs_obj->nrhs,
                                  dpbrfs_obj->af,
                                  dpbrfs_obj->lda,
                                  dpbrfs_obj->x,
                                  dpbrfs_obj->ldb
                                  );


    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    dpbrfs_obj->inforef = DPBRFS( dpbrfs_obj->matrix_layout,
                                  dpbrfs_obj->uplo,
                                  dpbrfs_obj->n,
                                  dpbrfs_obj->kd,
                                  dpbrfs_obj->nrhs,
                                  (const double*)dpbrfs_obj->aref,
                                  dpbrfs_obj->lda,
                                  dpbrfs_obj->afref, dpbrfs_obj->ldaf,
                                  dpbrfs_obj->bref, dpbrfs_obj->ldb,
                                  dpbrfs_obj->xref, dpbrfs_obj->ldx,
                                  dpbrfs_obj->ferrref,
                                  dpbrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    dpbrfs_obj->info = LAPACKE_dpbrfs( dpbrfs_obj->matrix_layout,
                                  dpbrfs_obj->uplo,
                                  dpbrfs_obj->n,
                                  dpbrfs_obj->kd,
                                  dpbrfs_obj->nrhs,
                                  (const double*)dpbrfs_obj->a,
                                  dpbrfs_obj->lda,
                                  dpbrfs_obj->af, dpbrfs_obj->ldaf,
                                  dpbrfs_obj->b, dpbrfs_obj->ldb,
                                  dpbrfs_obj->x, dpbrfs_obj->ldx,
                                  dpbrfs_obj->ferr,
                                  dpbrfs_obj->berr);

    if( dpbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dpbrfs is wrong\n", dpbrfs_obj->info );
    }
    if( dpbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpbrfs is wrong\n", 
        dpbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dpbrfs_obj->diff_xerr =  computeDiff_d( dpbrfs_obj->x_bufsize, 
                dpbrfs_obj->x, dpbrfs_obj->xref );

    dpbrfs_obj->diff_berr =  computeDiff_d( dpbrfs_obj->nrhs, 
                dpbrfs_obj->berr, dpbrfs_obj->berrref );
                
    dpbrfs_obj->diff_ferr =  computeDiff_d( dpbrfs_obj->nrhs, 
                dpbrfs_obj->ferr, dpbrfs_obj->ferrref );

}

TEST_F(dpbrfs_test, dpbrfs1) {
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);

}

TEST_F(dpbrfs_test, dpbrfs2) {
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dpbrfs_test, dpbrfs3) {
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dpbrfs_test, dpbrfs4) {
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin pbrfs_scomplex_parameters  class definition */
class pbrfs_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_float *a, *aref; //The array ab contains the matrix A
      lapack_complex_float *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pbrfs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i, lapack_int kd_i);
              
      ~pbrfs_scomplex_parameters (); 
};  /* end of pbrfs_scomplex_parameters  class definition */


/* Constructor pbrfs_scomplex_parameters definition */
pbrfs_scomplex_parameters:: pbrfs_scomplex_parameters ( int matrix_layout_i, 
                                          char uplo_i, lapack_int n_i,
                                lapack_int nrhs_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    n = n_i;
    nrhs = nrhs_i;
    kd = kd_i;
    lda = n;
    ldaf = n;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n pbrfs lapack_complex_float:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, uplo, lda, 
ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       pbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand ( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(lapack_complex_float)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, (b_bufsize*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

pbrfs_scomplex_parameters:: ~pbrfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbrfs_scomplex_parameters object: destructor invoked. \n");
#endif
   pbrfs_free();
}

//  Test fixture class definition
class cpbrfs_test  : public  ::testing::Test {
public:
   pbrfs_scomplex_parameters  *cpbrfs_obj;
   void SetUp();  
   void TearDown () { delete cpbrfs_obj; }
};


void cpbrfs_test::SetUp(){

    /* LAPACKE CPBRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_cpbrfs) (int matrix_layout, char uplo,
        lapack_int n, lapack_int kd, lapack_int nrhs, const lapack_complex_float* ab,
        lapack_int ldab, const lapack_complex_float* afb, lapack_int ldafb,
        const lapack_complex_float* b, lapack_int ldb, lapack_complex_float* x, lapack_int ldx,
        float* ferr, float* berr);

    Fptr_NL_LAPACKE_cpbrfs CPBRFS;

     /* LAPACKE CPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpbtrf) ( int matrix_layout , char uplo ,
                lapack_int n , lapack_int kd, lapack_complex_float * a , lapack_int lda );

    Fptr_NL_LAPACKE_cpbtrf CPBTRF;

    typedef int (*Fptr_NL_LAPACKE_cpbtrs) ( int matrix_layout , char uplo ,
        lapack_int n , lapack_int kd , lapack_int nrhs , const lapack_complex_float *ab ,
        lapack_int ldab , lapack_complex_float *b , lapack_int ldb );

    Fptr_NL_LAPACKE_cpbtrs CPBTRS;


    cpbrfs_obj = new pbrfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd);

    //idx = Circular_Increment_Index(idx);

    cpbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CPBTRF = (Fptr_NL_LAPACKE_cpbtrf)dlsym(cpbrfs_obj->hModule,"LAPACKE_cpbtrf");
    ASSERT_TRUE(CPBTRF != NULL) << "failed to pbt the Netlib LAPACKE_cpbtrf symbol";
    
    CPBTRS = (Fptr_NL_LAPACKE_cpbtrs)dlsym(cpbrfs_obj->hModule, "LAPACKE_cpbtrs");
    ASSERT_TRUE(CPBTRS != NULL) << "failed to pbt the Netlib LAPACKE_cpbtrs symbol";
    
    CPBRFS = (Fptr_NL_LAPACKE_cpbrfs)dlsym(cpbrfs_obj->hModule, "LAPACKE_cpbrfs");
    ASSERT_TRUE(CPBRFS != NULL) << "failed to pbt the Netlib LAPACKE_cpbrfs symbol";

    /* Preparing the i/ps for pbrfs by invoking PBTRF, PBTRS API calls */
    cpbrfs_obj->inforef = CPBTRF( cpbrfs_obj->matrix_layout,
                                  cpbrfs_obj->uplo, cpbrfs_obj->n,
                                  cpbrfs_obj->kd,
                                  cpbrfs_obj->afref,
                                  cpbrfs_obj->lda);
                           
    cpbrfs_obj->info = LAPACKE_cpbtrf( cpbrfs_obj->matrix_layout,
                                       cpbrfs_obj->uplo, cpbrfs_obj->n,
                                       cpbrfs_obj->kd,
                                       cpbrfs_obj->af,
                                       cpbrfs_obj->lda);

    cpbrfs_obj->inforef = CPBTRS( cpbrfs_obj->matrix_layout,
                                  cpbrfs_obj->uplo,
                                  cpbrfs_obj->n,
                                  cpbrfs_obj->kd,
                                  cpbrfs_obj->nrhs,
                                  cpbrfs_obj->afref,
                                  cpbrfs_obj->lda,
                                  cpbrfs_obj->xref,
                                  cpbrfs_obj->ldb
                                  );

    cpbrfs_obj->info = LAPACKE_cpbtrs( cpbrfs_obj->matrix_layout,
                                  cpbrfs_obj->uplo,
                                  cpbrfs_obj->n,
                                  cpbrfs_obj->kd,
                                  cpbrfs_obj->nrhs,
                                  cpbrfs_obj->af,
                                  cpbrfs_obj->lda,
                                  cpbrfs_obj->x,
                                  cpbrfs_obj->ldb
                                  );


    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    cpbrfs_obj->inforef = CPBRFS( cpbrfs_obj->matrix_layout,
                                  cpbrfs_obj->uplo,
                                  cpbrfs_obj->n,
                                  cpbrfs_obj->kd,
                                  cpbrfs_obj->nrhs,
                                  (const lapack_complex_float*)cpbrfs_obj->aref,
                                  cpbrfs_obj->lda,
                                  cpbrfs_obj->afref, cpbrfs_obj->ldaf,
                                  cpbrfs_obj->bref, cpbrfs_obj->ldb,
                                  cpbrfs_obj->xref, cpbrfs_obj->ldx,
                                  cpbrfs_obj->ferrref,
                                  cpbrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    cpbrfs_obj->info = LAPACKE_cpbrfs( cpbrfs_obj->matrix_layout,
                                  cpbrfs_obj->uplo,
                                  cpbrfs_obj->n,
                                  cpbrfs_obj->kd,
                                  cpbrfs_obj->nrhs,
                                  (const lapack_complex_float*)cpbrfs_obj->a,
                                  cpbrfs_obj->lda,
                                  cpbrfs_obj->af, cpbrfs_obj->ldaf,
                                  cpbrfs_obj->b, cpbrfs_obj->ldb,
                                  cpbrfs_obj->x, cpbrfs_obj->ldx,
                                  cpbrfs_obj->ferr,
                                  cpbrfs_obj->berr);

    if( cpbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cpbrfs is wrong\n", cpbrfs_obj->info );
    }
    if( cpbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpbrfs is wrong\n", 
        cpbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cpbrfs_obj->diff_xerr =  computeDiff_c( cpbrfs_obj->x_bufsize, 
                cpbrfs_obj->x, cpbrfs_obj->xref );

    cpbrfs_obj->diff_berr =  computeDiff_s( cpbrfs_obj->nrhs, 
                cpbrfs_obj->berr, cpbrfs_obj->berrref );
                
    cpbrfs_obj->diff_ferr =  computeDiff_s( cpbrfs_obj->nrhs, 
                cpbrfs_obj->ferr, cpbrfs_obj->ferrref );

}

TEST_F(cpbrfs_test, cpbrfs1) {
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);

}

TEST_F(cpbrfs_test, cpbrfs2) {
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(cpbrfs_test, cpbrfs3) {
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(cpbrfs_test, cpbrfs4) {
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, cpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin pbrfs_dcomplex_parameters  class definition */
class pbrfs_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int kd;// The number of subdiagonals within the band of A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldaf;  //  leading dimension of 'af'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_double *a, *aref; //The array ab contains the matrix A
      lapack_complex_double *af, *afref; //contains the factored form of the matrix A
      
      /* Output parameters */
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pbrfs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i, lapack_int kd_i);
              
      ~pbrfs_dcomplex_parameters (); 
};  /* end of pbrfs_dcomplex_parameters  class definition */


/* Constructor pbrfs_dcomplex_parameters definition */
pbrfs_dcomplex_parameters:: pbrfs_dcomplex_parameters ( int matrix_layout_i, 
                                          char uplo_i, lapack_int n_i,
                                lapack_int nrhs_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    n = n_i;
    nrhs = nrhs_i;
    kd = kd_i;
    lda = n;
    ldaf = n;
    
    b_bufsize = n*nrhs;
    x_bufsize = b_bufsize;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
        ldx = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
        ldx = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

#if LAPACKE_TEST_VERBOSE
   printf(" \n pbrfs lapack_complex_double:  n: %d, uplo: %c  lda: %d  \
ldb: %d nrhs: %d  ldaf: %d  ldx: %d \n",  n, uplo, lda, 
ldb, nrhs, ldaf, ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &af, &afref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       pbrfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand ( a, aref, n*n);
    memcpy(af, a, (n*n*sizeof(lapack_complex_double)));
    memcpy(afref, a, (n*n*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, (b_bufsize*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);

    
   } /* end of Constructor  */

pbrfs_dcomplex_parameters:: ~pbrfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbrfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   pbrfs_free();
}

//  Test fixture class definition
class zpbrfs_test  : public  ::testing::Test {
public:
   pbrfs_dcomplex_parameters  *zpbrfs_obj;
   void SetUp();  
   void TearDown () { delete zpbrfs_obj; }
};


void zpbrfs_test::SetUp(){

    /* LAPACKE ZPBRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_zpbrfs) (int matrix_layout, char uplo,
        lapack_int n, lapack_int kd, lapack_int nrhs, const lapack_complex_double* ab,
        lapack_int ldab, const lapack_complex_double* afb, lapack_int ldafb,
        const lapack_complex_double* b, lapack_int ldb, lapack_complex_double* x, lapack_int ldx,
        double* ferr, double* berr);

    Fptr_NL_LAPACKE_zpbrfs ZPBRFS;

     /* LAPACKE ZPBTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpbtrf) ( int matrix_layout , char uplo ,
                lapack_int n , lapack_int kd, lapack_complex_double * a , lapack_int lda );

    Fptr_NL_LAPACKE_zpbtrf ZPBTRF;

    typedef int (*Fptr_NL_LAPACKE_zpbtrs) ( int matrix_layout , char uplo ,
        lapack_int n , lapack_int kd , lapack_int nrhs , const lapack_complex_double *ab ,
        lapack_int ldab , lapack_complex_double *b , lapack_int ldb );

    Fptr_NL_LAPACKE_zpbtrs ZPBTRS;


    zpbrfs_obj = new pbrfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].kd);

    //idx = Circular_Increment_Index(idx);

    zpbrfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpbrfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpbrfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpbrfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZPBTRF = (Fptr_NL_LAPACKE_zpbtrf)dlsym(zpbrfs_obj->hModule,"LAPACKE_zpbtrf");
    ASSERT_TRUE(ZPBTRF != NULL) << "failed to pbt the Netlib LAPACKE_zpbtrf symbol";
    
    ZPBTRS = (Fptr_NL_LAPACKE_zpbtrs)dlsym(zpbrfs_obj->hModule, "LAPACKE_zpbtrs");
    ASSERT_TRUE(ZPBTRS != NULL) << "failed to pbt the Netlib LAPACKE_zpbtrs symbol";
    
    ZPBRFS = (Fptr_NL_LAPACKE_zpbrfs)dlsym(zpbrfs_obj->hModule, "LAPACKE_zpbrfs");
    ASSERT_TRUE(ZPBRFS != NULL) << "failed to pbt the Netlib LAPACKE_zpbrfs symbol";

    /* Preparing the i/ps for pbrfs by invoking PBTRF, PBTRS API calls */
    zpbrfs_obj->inforef = ZPBTRF( zpbrfs_obj->matrix_layout,
                                  zpbrfs_obj->uplo, zpbrfs_obj->n,
                                  zpbrfs_obj->kd,
                                  zpbrfs_obj->afref,
                                  zpbrfs_obj->lda);
                           
    zpbrfs_obj->info = LAPACKE_zpbtrf( zpbrfs_obj->matrix_layout,
                                       zpbrfs_obj->uplo, zpbrfs_obj->n,
                                       zpbrfs_obj->kd,
                                       zpbrfs_obj->af,
                                       zpbrfs_obj->lda);

    zpbrfs_obj->inforef = ZPBTRS( zpbrfs_obj->matrix_layout,
                                  zpbrfs_obj->uplo,
                                  zpbrfs_obj->n,
                                  zpbrfs_obj->kd,
                                  zpbrfs_obj->nrhs,
                                  zpbrfs_obj->afref,
                                  zpbrfs_obj->lda,
                                  zpbrfs_obj->xref,
                                  zpbrfs_obj->ldb
                                  );

    zpbrfs_obj->info = LAPACKE_zpbtrs( zpbrfs_obj->matrix_layout,
                                  zpbrfs_obj->uplo,
                                  zpbrfs_obj->n,
                                  zpbrfs_obj->kd,
                                  zpbrfs_obj->nrhs,
                                  zpbrfs_obj->af,
                                  zpbrfs_obj->lda,
                                  zpbrfs_obj->x,
                                  zpbrfs_obj->ldb
                                  );


    /* Compute the reference o/p by invoking Netlib-Lapack's API */
    zpbrfs_obj->inforef = ZPBRFS( zpbrfs_obj->matrix_layout,
                                  zpbrfs_obj->uplo,
                                  zpbrfs_obj->n,
                                  zpbrfs_obj->kd,
                                  zpbrfs_obj->nrhs,
                                  (const lapack_complex_double*)zpbrfs_obj->aref,
                                  zpbrfs_obj->lda,
                                  zpbrfs_obj->afref, zpbrfs_obj->ldaf,
                                  zpbrfs_obj->bref, zpbrfs_obj->ldb,
                                  zpbrfs_obj->xref, zpbrfs_obj->ldx,
                                  zpbrfs_obj->ferrref,
                                  zpbrfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zpbrfs_obj->info = LAPACKE_zpbrfs( zpbrfs_obj->matrix_layout,
                                  zpbrfs_obj->uplo,
                                  zpbrfs_obj->n,
                                  zpbrfs_obj->kd,
                                  zpbrfs_obj->nrhs,
                                  (const lapack_complex_double*)zpbrfs_obj->a,
                                  zpbrfs_obj->lda,
                                  zpbrfs_obj->af, zpbrfs_obj->ldaf,
                                  zpbrfs_obj->b, zpbrfs_obj->ldb,
                                  zpbrfs_obj->x, zpbrfs_obj->ldx,
                                  zpbrfs_obj->ferr,
                                  zpbrfs_obj->berr);

    if( zpbrfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zpbrfs is wrong\n", zpbrfs_obj->info );
    }
    if( zpbrfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpbrfs is wrong\n", 
        zpbrfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zpbrfs_obj->diff_xerr =  computeDiff_z( zpbrfs_obj->x_bufsize, 
                zpbrfs_obj->x, zpbrfs_obj->xref );

    zpbrfs_obj->diff_berr =  computeDiff_d( zpbrfs_obj->nrhs, 
                zpbrfs_obj->berr, zpbrfs_obj->berrref );
                
    zpbrfs_obj->diff_ferr =  computeDiff_d( zpbrfs_obj->nrhs, 
                zpbrfs_obj->ferr, zpbrfs_obj->ferrref );

}

TEST_F(zpbrfs_test, zpbrfs1) {
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);

}

TEST_F(zpbrfs_test, zpbrfs2) {
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zpbrfs_test, zpbrfs3) {
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zpbrfs_test, zpbrfs4) {
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_xerr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_berr, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zpbrfs_obj->diff_ferr, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}
