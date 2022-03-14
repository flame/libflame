
#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define hprfs_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (x != NULL)    free (x  ); \
  if (xref != NULL) free (xref); \
  if (af != NULL)    free (af  ); \
  if (afref != NULL) free (afref); \
  if (ipiv != NULL)    free (ipiv  ); \
  if (ipivref != NULL) free (ipivref); \
  if (ferr != NULL)    free (ferr  ); \
  if (ferrref != NULL) free (ferrref); \
  if (berr != NULL)    free (berr  ); \
  if (berrref != NULL) free (berrref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin hprfs_scomplex_parameters  class definition */
class hprfs_scomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      float diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
      float threshold;
      
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_float *a, *aref; //The array ab contains the matrix A
      lapack_complex_float *af, *afref; //contains the factored form of the matrix A
      lapack_int *ipiv, *ipivref; //  leading dimension of 'x'
      
      /* Output parameters */
      lapack_complex_float* x, *xref;
      float* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      float* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      hprfs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~hprfs_scomplex_parameters (); 
};  /* end of hprfs_scomplex_parameters  class definition */


/* Constructor hprfs_scomplex_parameters definition */
hprfs_scomplex_parameters:: hprfs_scomplex_parameters ( int matrix_layout_i, 
                char uplo_i, lapack_int n_i, lapack_int nrhs_i) {

    int a_buf_size =  n_i*(n_i+1)/2;    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;

    n = n_i;
    nrhs = nrhs_i;
    
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
    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;
    diff = 0;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n hprfs lapack_complex_float:  n: %d,  uplo: %c   \
ldb: %d nrhs: %d   ldx: %d \n",  n, uplo, ldb, nrhs,  ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_buf_size);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &af, &afref, a_buf_size);
    lapacke_gtest_alloc_float_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_float_buffer_pair( &berr, &berrref, nrhs);
    lapacke_gtest_alloc_int_buffer_pair( &ipiv, &ipivref, n);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       hprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, a_buf_size);
    memcpy(af, a, (a_buf_size*sizeof(lapack_complex_float)));
    memcpy(afref, a, (a_buf_size*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(lapack_complex_float)));
    memcpy(xref, b, (b_bufsize*sizeof(lapack_complex_float)));
    
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);

    
   } /* end of Constructor  */

hprfs_scomplex_parameters:: ~hprfs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hprfs_scomplex_parameters object: destructor invoked. \n");
#endif
   hprfs_free();
}

//  Test fixture class definition
class chprfs_test  : public  ::testing::Test {
public:
   hprfs_scomplex_parameters  *chprfs_obj;
   void SetUp();  
   void TearDown () { delete chprfs_obj; }
};


void chprfs_test::SetUp(){

    /* LAPACKE CHPRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_chprfs) (int matrix_layout, char uplo,
        lapack_int n, lapack_int nrhs, const lapack_complex_float* ap,
        const lapack_complex_float* afp, const lapack_int* ipiv,
        const lapack_complex_float* b, lapack_int ldb,
        lapack_complex_float* x, lapack_int ldx, float* ferr, float* berr);

    Fptr_NL_LAPACKE_chprfs CHPRFS;

    /* LAPACKE CHPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_chptrs) ( int matrix_layout , char uplo ,
        lapack_int n , lapack_int nrhs , const lapack_complex_float * ap ,
        const lapack_int * ipiv , lapack_complex_float * b , lapack_int ldb );
        
    Fptr_NL_LAPACKE_chptrs CHPTRS;
    
     /* LAPACKE CHPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chptrf) ( int matrix_layout , char uplo ,
            lapack_int n , lapack_complex_float * ap , lapack_int * ipiv );

    Fptr_NL_LAPACKE_chptrf CHPTRF;

    chprfs_obj = new hprfs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    chprfs_obj->threshold = lin_solver_paramslist[idx].solver_threhold;

    idx = Circular_Increment_Index(idx);

    chprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CHPRFS = (Fptr_NL_LAPACKE_chprfs)dlsym(chprfs_obj->hModule, "LAPACKE_chprfs");
    ASSERT_TRUE(CHPRFS != NULL) << "failed to ppt the Netlib LAPACKE_chprfs symbol";

    CHPTRS = (Fptr_NL_LAPACKE_chptrs)dlsym(chprfs_obj->hModule, "LAPACKE_chptrs");
    ASSERT_TRUE(CHPTRS != NULL) << "failed to get the Netlib LAPACKE_chptrs symbol";

    CHPTRF = (Fptr_NL_LAPACKE_chptrf)dlsym(chprfs_obj->hModule,"LAPACKE_chptrf");
    ASSERT_TRUE(CHPTRF != NULL) << "failed to get the Netlib LAPACKE_chptrf symbol";
    
    /* invoke the hptrf, hptrs APIs to compute i/ps for HPRFS calls.  */
    chprfs_obj->inforef = CHPTRF( chprfs_obj->matrix_layout,
                            chprfs_obj->uplo, chprfs_obj->n,
                            chprfs_obj->afref,
							chprfs_obj->ipivref);
                               
    chprfs_obj->info = LAPACKE_chptrf( chprfs_obj->matrix_layout,
                         chprfs_obj->uplo, chprfs_obj->n,
                               chprfs_obj->af,
							chprfs_obj->ipiv);

    chprfs_obj->inforef = CHPTRS( chprfs_obj->matrix_layout, 
          chprfs_obj->uplo, chprfs_obj->n, chprfs_obj->nrhs,
			 (const lapack_complex_float *)chprfs_obj->aref, 
       chprfs_obj->ipivref, chprfs_obj->xref, chprfs_obj->ldx);

    chprfs_obj->info = LAPACKE_chptrs( chprfs_obj->matrix_layout, 
               chprfs_obj->uplo, chprfs_obj->n, chprfs_obj->nrhs,
					 (const lapack_complex_float *)chprfs_obj->a, 
               chprfs_obj->ipiv, chprfs_obj->x, chprfs_obj->ldb );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    chprfs_obj->inforef = CHPRFS( chprfs_obj->matrix_layout,
                                  chprfs_obj->uplo, chprfs_obj->n,
                                  chprfs_obj->nrhs,
                                  chprfs_obj->aref, 
                                  chprfs_obj->afref,
                                  chprfs_obj->ipivref,
                                  chprfs_obj->bref, chprfs_obj->ldb,
                                  chprfs_obj->xref, chprfs_obj->ldx,
                                  chprfs_obj->ferrref,
                                  chprfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    chprfs_obj->info = LAPACKE_chprfs( chprfs_obj->matrix_layout,
                                  chprfs_obj->uplo, chprfs_obj->n,
                                  chprfs_obj->nrhs,
                                  chprfs_obj->a, 
                                  chprfs_obj->af,
                                  chprfs_obj->ipiv,
                                  chprfs_obj->b, chprfs_obj->ldb,
                                  chprfs_obj->x, chprfs_obj->ldx,
                                  chprfs_obj->ferr,
                                  chprfs_obj->berr);

    if( chprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chprfs is wrong\n", chprfs_obj->info );
    }
    if( chprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chprfs is wrong\n", 
        chprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chprfs_obj->diff_xerr =  computeDiff_c( chprfs_obj->x_bufsize, 
                chprfs_obj->x, chprfs_obj->xref );

    chprfs_obj->diff_berr =  computeDiff_s( chprfs_obj->nrhs, 
                chprfs_obj->berr, chprfs_obj->berrref );
                
    chprfs_obj->diff_ferr =  computeDiff_s( chprfs_obj->nrhs, 
                chprfs_obj->ferr, chprfs_obj->ferrref );
    
}

TEST_F(chprfs_test, chprfs1) {
    EXPECT_NEAR(0.0, chprfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chprfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chprfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chprfs_test, chprfs2) {
    EXPECT_NEAR(0.0, chprfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chprfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chprfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chprfs_test, chprfs3) {
    EXPECT_NEAR(0.0, chprfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chprfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chprfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chprfs_test, chprfs4) {
    EXPECT_NEAR(0.0, chprfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chprfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, chprfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}


/* Begin hprfs_dcomplex_parameters  class definition */
class hprfs_dcomplex_parameters{
   public:
      int b_bufsize;
      int x_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      double diff_berr, diff_ferr, diff_xerr;
      void *hModule, *dModule;
      float threshold;
      
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.

      lapack_int n; // The order of A (Square matrix).
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_int ldx; //  leading dimension of 'x'

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.
      lapack_complex_double *a, *aref; //The array ab contains the matrix A
      lapack_complex_double *af, *afref; //contains the factored form of the matrix A
      lapack_int *ipiv, *ipivref; //  leading dimension of 'x'
      
      /* Output parameters */
      lapack_complex_double* x, *xref;
      double* ferr, *ferrref; // estimated forward error bound for each solution vector x.
      double* berr, *berrref; // component-wise relative backward error for each solution vector xj.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      hprfs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
              
      ~hprfs_dcomplex_parameters (); 
};  /* end of hprfs_dcomplex_parameters  class definition */


/* Constructor hprfs_dcomplex_parameters definition */
hprfs_dcomplex_parameters:: hprfs_dcomplex_parameters ( int matrix_layout_i, 
                char uplo_i, lapack_int n_i, lapack_int nrhs_i) {

    int a_buf_size =  n_i*(n_i+1)/2;    
    matrix_layout = matrix_layout_i;
    uplo = uplo_i;

    n = n_i;
    nrhs = nrhs_i;
    
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
    hModule = NULL;
    dModule = NULL;
    diff_berr = 0;
    diff_ferr = 0;
    diff_xerr = 0;
    diff = 0;
    
#if LAPACKE_TEST_VERBOSE
   printf(" \n hprfs lapack_complex_double:  n: %d,  uplo: %c   \
ldb: %d nrhs: %d   ldx: %d \n",  n, uplo, ldb, nrhs,  ldx);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_buf_size);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &x, &xref, x_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &af, &afref, a_buf_size);
    lapacke_gtest_alloc_double_buffer_pair( &ferr, &ferrref, nrhs);
    lapacke_gtest_alloc_double_buffer_pair( &berr, &berrref, nrhs);
    lapacke_gtest_alloc_int_buffer_pair( &ipiv, &ipivref, n);
    
    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (x==NULL) || (xref==NULL) || \
        (af==NULL) || (afref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) || \
        (ferr==NULL) || (ferrref==NULL) || \
        (berr==NULL) || (berrref==NULL) ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       hprfs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, a_buf_size);
    memcpy(af, a, (a_buf_size*sizeof(lapack_complex_double)));
    memcpy(afref, a, (a_buf_size*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    memcpy(x, b, (b_bufsize*sizeof(lapack_complex_double)));
    memcpy(xref, b, (b_bufsize*sizeof(lapack_complex_double)));
    
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(x, xref, x_bufsize, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(ferr, ferrref, nrhs, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant(berr, berrref, nrhs, 0.0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);

    
   } /* end of Constructor  */

hprfs_dcomplex_parameters:: ~hprfs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hprfs_dcomplex_parameters object: destructor invoked. \n");
#endif
   hprfs_free();
}

//  Test fixture class definition
class zhprfs_test  : public  ::testing::Test {
public:
   hprfs_dcomplex_parameters  *zhprfs_obj;
   void SetUp();  
   void TearDown () { delete zhprfs_obj; }
};


void zhprfs_test::SetUp(){

    /* LAPACKE ZHPRFS prototype */
    typedef int (*Fptr_NL_LAPACKE_zhprfs) (int matrix_layout, char uplo,
        lapack_int n, lapack_int nrhs, const lapack_complex_double* ap,
        const lapack_complex_double* afp, const lapack_int* ipiv,
        const lapack_complex_double* b, lapack_int ldb,
        lapack_complex_double* x, lapack_int ldx, double* ferr, double* berr);

    Fptr_NL_LAPACKE_zhprfs ZHPRFS;

    /* LAPACKE ZHPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zhptrs) ( int matrix_layout , char uplo ,
        lapack_int n , lapack_int nrhs , const lapack_complex_double * ap ,
        const lapack_int * ipiv , lapack_complex_double * b , lapack_int ldb );
        
    Fptr_NL_LAPACKE_zhptrs ZHPTRS;
    
     /* LAPACKE ZHPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhptrf) ( int matrix_layout , char uplo ,
            lapack_int n , lapack_complex_double * ap , lapack_int * ipiv );

    Fptr_NL_LAPACKE_zhptrf ZHPTRF;

    zhprfs_obj = new hprfs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs);

    zhprfs_obj->threshold = lin_solver_paramslist[idx].solver_threhold;

    idx = Circular_Increment_Index(idx);

    zhprfs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhprfs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhprfs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhprfs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZHPRFS = (Fptr_NL_LAPACKE_zhprfs)dlsym(zhprfs_obj->hModule, "LAPACKE_zhprfs");
    ASSERT_TRUE(ZHPRFS != NULL) << "failed to ppt the Netlib LAPACKE_zhprfs symbol";

    ZHPTRS = (Fptr_NL_LAPACKE_zhptrs)dlsym(zhprfs_obj->hModule, "LAPACKE_zhptrs");
    ASSERT_TRUE(ZHPTRS != NULL) << "failed to get the Netlib LAPACKE_zhptrs symbol";

    ZHPTRF = (Fptr_NL_LAPACKE_zhptrf)dlsym(zhprfs_obj->hModule,"LAPACKE_zhptrf");
    ASSERT_TRUE(ZHPTRF != NULL) << "failed to get the Netlib LAPACKE_zhptrf symbol";
    
    /* invoke the hptrf, hptrs APIs to compute i/ps for HPRFS calls.  */
    zhprfs_obj->inforef = ZHPTRF( zhprfs_obj->matrix_layout,
                            zhprfs_obj->uplo, zhprfs_obj->n,
                            zhprfs_obj->afref,
							zhprfs_obj->ipivref);
                               
    //zhprfs_obj->info = ZHPTRF( zhprfs_obj->matrix_layout,
    zhprfs_obj->info = LAPACKE_zhptrf( zhprfs_obj->matrix_layout,
                         zhprfs_obj->uplo, zhprfs_obj->n,
                               zhprfs_obj->af,
							zhprfs_obj->ipiv);

    zhprfs_obj->inforef = ZHPTRS( zhprfs_obj->matrix_layout, 
          zhprfs_obj->uplo, zhprfs_obj->n, zhprfs_obj->nrhs,
			 (const lapack_complex_double *)zhprfs_obj->aref, 
       zhprfs_obj->ipivref, zhprfs_obj->xref, zhprfs_obj->ldx);

    zhprfs_obj->info = LAPACKE_zhptrs( zhprfs_obj->matrix_layout, 
               zhprfs_obj->uplo, zhprfs_obj->n, zhprfs_obj->nrhs,
					 (const lapack_complex_double *)zhprfs_obj->a, 
               zhprfs_obj->ipiv, zhprfs_obj->x, zhprfs_obj->ldb );

    /* Compute the reference o/p by invoking Netlib-Lapack's API */ 
    zhprfs_obj->inforef = ZHPRFS( zhprfs_obj->matrix_layout,
                                  zhprfs_obj->uplo, zhprfs_obj->n,
                                  zhprfs_obj->nrhs,
                                  zhprfs_obj->aref, 
                                  zhprfs_obj->afref,
                                  zhprfs_obj->ipivref,
                                  zhprfs_obj->bref, zhprfs_obj->ldb,
                                  zhprfs_obj->xref, zhprfs_obj->ldx,
                                  zhprfs_obj->ferrref,
                                  zhprfs_obj->berrref);

    /* Compute libflame's Lapacke o/p  */
    zhprfs_obj->info = LAPACKE_zhprfs( zhprfs_obj->matrix_layout,
                                  zhprfs_obj->uplo, zhprfs_obj->n,
                                  zhprfs_obj->nrhs,
                                  zhprfs_obj->a, 
                                  zhprfs_obj->af,
                                  zhprfs_obj->ipiv,
                                  zhprfs_obj->b, zhprfs_obj->ldb,
                                  zhprfs_obj->x, zhprfs_obj->ldx,
                                  zhprfs_obj->ferr,
                                  zhprfs_obj->berr);

    if( zhprfs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhprfs is wrong\n", zhprfs_obj->info );
    }
    if( zhprfs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhprfs is wrong\n", 
        zhprfs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhprfs_obj->diff_xerr =  computeDiff_z( zhprfs_obj->x_bufsize, 
                zhprfs_obj->x, zhprfs_obj->xref );

    zhprfs_obj->diff_berr =  computeDiff_d( zhprfs_obj->nrhs, 
                zhprfs_obj->berr, zhprfs_obj->berrref );
                
    zhprfs_obj->diff_ferr =  computeDiff_d( zhprfs_obj->nrhs, 
                zhprfs_obj->ferr, zhprfs_obj->ferrref );
    
}

TEST_F(zhprfs_test, zhprfs1) {
    EXPECT_NEAR(0.0, zhprfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhprfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhprfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhprfs_test, zhprfs2) {
    EXPECT_NEAR(0.0, zhprfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhprfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhprfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhprfs_test, zhprfs3) {
    EXPECT_NEAR(0.0, zhprfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhprfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhprfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhprfs_test, zhprfs4) {
    EXPECT_NEAR(0.0, zhprfs_obj->diff_xerr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhprfs_obj->diff_berr, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(0.0, zhprfs_obj->diff_ferr, LAPACKE_GTEST_THRESHOLD);
}
