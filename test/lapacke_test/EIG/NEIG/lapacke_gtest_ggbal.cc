#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define ggbal_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (lscale!=NULL)   free(lscale); \
    if (lscaleref!=NULL)  free(lscaleref); \
    if (rscale!=NULL)        free(rscale); \
    if (rscaleref!=NULL)     free(rscaleref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ggbal_float_common_parameters  class definition */
class ggbal_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_lscale, diff_rscale;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'N' or 'P' or 'S' or 'B'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b

    /* Input / Output parameters */
    float* a, *aref; // contains the n-by-n general matrix A.
    float* b, *bref; // contains the n-by-n upper triangular matrix B.
    float *lscale, *lscaleref;
    float *rscale, *rscaleref;

    /* Input / Output parameters */
    lapack_int ilo, iloref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi, ihiref;

    /*Return Values */
    int info, inforef;

   public:
      ggbal_float_parameters (int matrix_layout_i, char job, lapack_int n);
      ~ggbal_float_parameters ();
};

/* Constructor definition  float_common_parameters */
ggbal_float_parameters:: ggbal_float_parameters (int matrix_layout_i,
                                          char job_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    n = n_i;

    lda = n;
    ldb = n;

    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_b = 0;
    diff_lscale = 0;
    diff_rscale = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
    iloref = 0;
    ihiref = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbal float: matrix_layout: %d n: %d  job: %c \n",
                                 matrix_layout, n, job);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_float_buffer_pair( &rscale, &rscaleref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) ){
       EXPECT_FALSE( true) << "ggbal_float_parameters object: malloc error. Exiting ";
       ggbal_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    //lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'ggbal_float_common_parameters' */
ggbal_float_parameters :: ~ggbal_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggbal_free();
} 

//  Test fixture class definition
class sggbal_test  : public  ::testing::Test {
public:
   ggbal_float_parameters  *sggbal_obj;
   void SetUp();  
   void TearDown () { delete sggbal_obj; }
};

void sggbal_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_sggbal) ( int matrix_layout, char job,
                    lapack_int n, float* a, lapack_int lda, float* b,
                    lapack_int ldb, lapack_int* ilo, lapack_int* ihi,
                    float* lscale, float* rscale);
                 
    Fptr_NL_LAPACKE_sggbal SGGBAL;

    sggbal_obj = new  ggbal_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].n );
                                         
    idx = Circular_Increment_Index(idx);
    sggbal_obj->threshold = eig_non_sym_paramslist[idx].ggbal_threshold;
    sggbal_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sggbal_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sggbal_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sggbal_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGBAL = (Fptr_NL_LAPACKE_sggbal)dlsym(sggbal_obj->hModule, "LAPACKE_sggbal");
    ASSERT_TRUE(SGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_sggbal symbol";


    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sggbal_obj->inforef = SGGBAL(   sggbal_obj->matrix_layout,
                                    sggbal_obj->job,
                                    sggbal_obj->n,
                                    sggbal_obj->aref,
                                    sggbal_obj->lda,
                                    sggbal_obj->bref,
                                    sggbal_obj->ldb,
                                    &sggbal_obj->iloref,
                                    &sggbal_obj->ihiref,
                                    sggbal_obj->lscaleref,
                                    sggbal_obj->rscaleref
                                    );

    /* Compute libflame's Lapacke o/p  */
    sggbal_obj->info = LAPACKE_sggbal(  sggbal_obj->matrix_layout,
                                    sggbal_obj->job,
                                    sggbal_obj->n,
                                    sggbal_obj->a,
                                    sggbal_obj->lda,
                                    sggbal_obj->b,
                                    sggbal_obj->ldb,
                                    &sggbal_obj->ilo,
                                    &sggbal_obj->ihi,
                                    sggbal_obj->lscale,
                                    sggbal_obj->rscale
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sggbal_obj->diff_a =  computeDiff_s( (sggbal_obj->lda)*(sggbal_obj->n), 
                sggbal_obj->a, sggbal_obj->aref );

    sggbal_obj->diff_b =  computeDiff_s( (sggbal_obj->ldb)*(sggbal_obj->n), 
                sggbal_obj->b, sggbal_obj->bref );

    sggbal_obj->diff_lscale =  computeDiff_s( sggbal_obj->n, 
                sggbal_obj->lscale, sggbal_obj->lscaleref );

    sggbal_obj->diff_rscale =  computeDiff_s( sggbal_obj->n, 
                sggbal_obj->rscale, sggbal_obj->rscaleref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbal float: \n diff_a: %f \n diff_b: %f \n \
diff_lscale: %f \n diff_rscale: %f \n ilo: %d \n iloref: %d \
diff_ilo: %d \n ihi: %d \n ihiref: %d \n diff_ihi: %d \n",
       sggbal_obj->diff_a, sggbal_obj->diff_b, sggbal_obj->diff_lscale,
       sggbal_obj->diff_rscale, sggbal_obj->ilo, sggbal_obj->iloref,
       (sggbal_obj->ilo - sggbal_obj->iloref),
       sggbal_obj->ihi, sggbal_obj->ihiref,
       (sggbal_obj->ihi - sggbal_obj->ihiref) );
#endif
}

TEST_F(sggbal_test, sggbal1) {
    EXPECT_NEAR(0.0, sggbal_obj->diff_a, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_b, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_lscale, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_rscale, sggbal_obj->threshold);
    EXPECT_EQ(sggbal_obj->ilo, sggbal_obj->iloref);
    EXPECT_EQ(sggbal_obj->ihi, sggbal_obj->ihiref);
}

TEST_F(sggbal_test, sggbal2) {
    EXPECT_NEAR(0.0, sggbal_obj->diff_a, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_b, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_lscale, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_rscale, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, (sggbal_obj->ilo - sggbal_obj->iloref), sggbal_obj->threshold);
    EXPECT_NEAR(0.0, (sggbal_obj->ihi - sggbal_obj->ihiref), sggbal_obj->threshold);
}

TEST_F(sggbal_test, sggbal3) {
    EXPECT_NEAR(0.0, sggbal_obj->diff_a, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_b, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_lscale, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_rscale, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, (sggbal_obj->ilo - sggbal_obj->iloref), sggbal_obj->threshold);
    EXPECT_NEAR(0.0, (sggbal_obj->ihi - sggbal_obj->ihiref), sggbal_obj->threshold);
}

TEST_F(sggbal_test, sggbal4) {
    EXPECT_NEAR(0.0, sggbal_obj->diff_a, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_b, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_lscale, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, sggbal_obj->diff_rscale, sggbal_obj->threshold);
    EXPECT_NEAR(0.0, (sggbal_obj->ilo - sggbal_obj->iloref), sggbal_obj->threshold);
    EXPECT_NEAR(0.0, (sggbal_obj->ihi - sggbal_obj->ihiref), sggbal_obj->threshold);
}

/* Begin ggbal_double_common_parameters  class definition */
class ggbal_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b, diff_lscale, diff_rscale;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'N' or 'P' or 'S' or 'B'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b

    /* Input / Output parameters */
    double* a, *aref; // contains the n-by-n general matrix A.
    double* b, *bref; // contains the n-by-n upper triangular matrix B.
    double *lscale, *lscaleref;
    double *rscale, *rscaleref;

    /* Input / Output parameters */
    lapack_int ilo, iloref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi, ihiref;

    /*Return Values */
    int info, inforef;

   public:
      ggbal_double_parameters (int matrix_layout_i, char job, lapack_int n);
      ~ggbal_double_parameters ();
};

/* Constructor definition  double_common_parameters */
ggbal_double_parameters:: ggbal_double_parameters (int matrix_layout_i,
                                          char job_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    n = n_i;

    lda = n;
    ldb = n;

    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_b = 0;
    diff_lscale = 0;
    diff_rscale = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
    iloref = 0;
    ihiref = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbal double: matrix_layout: %d n: %d  job: %c \n",
                                 matrix_layout, n, job);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_double_buffer_pair( &rscale, &rscaleref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) ){
       EXPECT_FALSE( true) << "ggbal_double_parameters object: malloc error. Exiting ";
       ggbal_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, lda*n );
    //lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'ggbal_double_common_parameters' */
ggbal_double_parameters :: ~ggbal_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggbal_free();
} 

//  Test fixture class definition
class dggbal_test  : public  ::testing::Test {
public:
   ggbal_double_parameters  *dggbal_obj;
   void SetUp();  
   void TearDown () { delete dggbal_obj; }
};

void dggbal_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_dggbal) ( int matrix_layout, char job,
                    lapack_int n, double* a, lapack_int lda, double* b,
                    lapack_int ldb, lapack_int* ilo, lapack_int* ihi,
                    double* lscale, double* rscale);
                 
    Fptr_NL_LAPACKE_dggbal DGGBAL;

    dggbal_obj = new ggbal_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].n );
                                         
    idx = Circular_Increment_Index(idx);
    dggbal_obj->threshold = eig_non_sym_paramslist[idx].ggbal_threshold;
    dggbal_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dggbal_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dggbal_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dggbal_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGBAL = (Fptr_NL_LAPACKE_dggbal)dlsym(dggbal_obj->hModule, "LAPACKE_dggbal");
    ASSERT_TRUE(DGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_dggbal symbol";


    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dggbal_obj->inforef = DGGBAL(   dggbal_obj->matrix_layout,
                                    dggbal_obj->job,
                                    dggbal_obj->n,
                                    dggbal_obj->aref,
                                    dggbal_obj->lda,
                                    dggbal_obj->bref,
                                    dggbal_obj->ldb,
                                    &dggbal_obj->iloref,
                                    &dggbal_obj->ihiref,
                                    dggbal_obj->lscaleref,
                                    dggbal_obj->rscaleref
                                    );

    /* Compute libflame's Lapacke o/p  */
    dggbal_obj->info = LAPACKE_dggbal(  dggbal_obj->matrix_layout,
                                    dggbal_obj->job,
                                    dggbal_obj->n,
                                    dggbal_obj->a,
                                    dggbal_obj->lda,
                                    dggbal_obj->b,
                                    dggbal_obj->ldb,
                                    &dggbal_obj->ilo,
                                    &dggbal_obj->ihi,
                                    dggbal_obj->lscale,
                                    dggbal_obj->rscale
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    dggbal_obj->diff_a =  computeDiff_d( (dggbal_obj->lda)*(dggbal_obj->n), 
                dggbal_obj->a, dggbal_obj->aref );

    dggbal_obj->diff_b =  computeDiff_d( (dggbal_obj->ldb)*(dggbal_obj->n), 
                dggbal_obj->b, dggbal_obj->bref );

    dggbal_obj->diff_lscale =  computeDiff_d( dggbal_obj->n, 
                dggbal_obj->lscale, dggbal_obj->lscaleref );

    dggbal_obj->diff_rscale =  computeDiff_d( dggbal_obj->n, 
                dggbal_obj->rscale, dggbal_obj->rscaleref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbal double: \n diff_a: %f \n diff_b: %f \n \
diff_lscale: %f \n diff_rscale: %f \n ilo: %d \n iloref: %d \
diff_ilo: %d \n ihi: %d \n ihiref: %d \n diff_ihi: %d \n",
       dggbal_obj->diff_a, dggbal_obj->diff_b, dggbal_obj->diff_lscale,
       dggbal_obj->diff_rscale, dggbal_obj->ilo, dggbal_obj->iloref,
       (dggbal_obj->ilo - dggbal_obj->iloref),
       dggbal_obj->ihi, dggbal_obj->ihiref,
       (dggbal_obj->ihi - dggbal_obj->ihiref) );
#endif
}

TEST_F(dggbal_test, dggbal1) {
    EXPECT_NEAR(0.0, dggbal_obj->diff_a, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_b, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_lscale, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_rscale, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, (dggbal_obj->ilo - dggbal_obj->iloref), dggbal_obj->threshold);
    EXPECT_NEAR(0.0, (dggbal_obj->ihi - dggbal_obj->ihiref), dggbal_obj->threshold);
}

TEST_F(dggbal_test, dggbal2) {
    EXPECT_NEAR(0.0, dggbal_obj->diff_a, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_b, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_lscale, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_rscale, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, (dggbal_obj->ilo - dggbal_obj->iloref), dggbal_obj->threshold);
    EXPECT_NEAR(0.0, (dggbal_obj->ihi - dggbal_obj->ihiref), dggbal_obj->threshold);
}

TEST_F(dggbal_test, dggbal3) {
    EXPECT_NEAR(0.0, dggbal_obj->diff_a, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_b, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_lscale, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_rscale, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, (dggbal_obj->ilo - dggbal_obj->iloref), dggbal_obj->threshold);
    EXPECT_NEAR(0.0, (dggbal_obj->ihi - dggbal_obj->ihiref), dggbal_obj->threshold);
}

TEST_F(dggbal_test, dggbal4) {
    EXPECT_NEAR(0.0, dggbal_obj->diff_a, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_b, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_lscale, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, dggbal_obj->diff_rscale, dggbal_obj->threshold);
    EXPECT_NEAR(0.0, (dggbal_obj->ilo - dggbal_obj->iloref), dggbal_obj->threshold);
    EXPECT_NEAR(0.0, (dggbal_obj->ihi - dggbal_obj->ihiref), dggbal_obj->threshold);
}

/* Begin ggbal_scomplex_common_parameters  class definition */
class ggbal_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_lscale, diff_rscale;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'N' or 'P' or 'S' or 'B'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b

    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_float* b, *bref; // contains the n-by-n upper triangular matrix B.
    float *lscale, *lscaleref;
    float *rscale, *rscaleref;

    /* Input / Output parameters */
    lapack_int ilo, iloref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi, ihiref;

    /*Return Values */
    int info, inforef;

   public:
      ggbal_scomplex_parameters (int matrix_layout_i, char job, lapack_int n);
      ~ggbal_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
ggbal_scomplex_parameters:: ggbal_scomplex_parameters (int matrix_layout_i,
                                          char job_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    n = n_i;

    lda = n;
    ldb = n;

    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_b = 0;
    diff_lscale = 0;
    diff_rscale = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
    iloref = 0;
    ihiref = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbal lapack_complex_float: matrix_layout: %d n: %d  job: %c \n",
                                 matrix_layout, n, job);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_float_buffer_pair( &rscale, &rscaleref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) ){
       EXPECT_FALSE( true) << "ggbal_scomplex_parameters object: malloc error. Exiting ";
       ggbal_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'ggbal_scomplex_common_parameters' */
ggbal_scomplex_parameters :: ~ggbal_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggbal_free();
} 

//  Test fixture class definition
class cggbal_test  : public  ::testing::Test {
public:
   ggbal_scomplex_parameters  *cggbal_obj;
   void SetUp();  
   void TearDown () { delete cggbal_obj; }
};

void cggbal_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_cggbal) ( int matrix_layout, char job,
                    lapack_int n, lapack_complex_float* a, lapack_int lda, lapack_complex_float* b,
                    lapack_int ldb, lapack_int* ilo, lapack_int* ihi,
                    float* lscale, float* rscale);
                 
    Fptr_NL_LAPACKE_cggbal CGGBAL;

    cggbal_obj = new  ggbal_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].n );
                                         
    idx = Circular_Increment_Index(idx);
    cggbal_obj->threshold = eig_non_sym_paramslist[idx].ggbal_threshold;
    cggbal_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cggbal_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cggbal_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cggbal_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGBAL = (Fptr_NL_LAPACKE_cggbal)dlsym(cggbal_obj->hModule, "LAPACKE_cggbal");
    ASSERT_TRUE(CGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_cggbal symbol";


    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cggbal_obj->inforef = CGGBAL(   cggbal_obj->matrix_layout,
                                    cggbal_obj->job,
                                    cggbal_obj->n,
                                    cggbal_obj->aref,
                                    cggbal_obj->lda,
                                    cggbal_obj->bref,
                                    cggbal_obj->ldb,
                                    &cggbal_obj->iloref,
                                    &cggbal_obj->ihiref,
                                    cggbal_obj->lscaleref,
                                    cggbal_obj->rscaleref
                                    );

    /* Compute libflame's Lapacke o/p  */
    cggbal_obj->info = LAPACKE_cggbal(  cggbal_obj->matrix_layout,
                                    cggbal_obj->job,
                                    cggbal_obj->n,
                                    cggbal_obj->a,
                                    cggbal_obj->lda,
                                    cggbal_obj->b,
                                    cggbal_obj->ldb,
                                    &cggbal_obj->ilo,
                                    &cggbal_obj->ihi,
                                    cggbal_obj->lscale,
                                    cggbal_obj->rscale
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    cggbal_obj->diff_a =  computeDiff_c( (cggbal_obj->lda)*(cggbal_obj->n), 
                cggbal_obj->a, cggbal_obj->aref );

    cggbal_obj->diff_b =  computeDiff_c( (cggbal_obj->ldb)*(cggbal_obj->n), 
                cggbal_obj->b, cggbal_obj->bref );

    cggbal_obj->diff_lscale =  computeDiff_s( cggbal_obj->n, 
                cggbal_obj->lscale, cggbal_obj->lscaleref );

    cggbal_obj->diff_rscale =  computeDiff_s( cggbal_obj->n, 
                cggbal_obj->rscale, cggbal_obj->rscaleref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbal lapack_complex_float: \n diff_a: %f \n diff_b: %f \n \
diff_lscale: %f \n diff_rscale: %f \n ilo: %d \n iloref: %d \
diff_ilo: %d \n ihi: %d \n ihiref: %d \n diff_ihi: %d \n",
       cggbal_obj->diff_a, cggbal_obj->diff_b, cggbal_obj->diff_lscale,
       cggbal_obj->diff_rscale, cggbal_obj->ilo, cggbal_obj->iloref,
       (cggbal_obj->ilo - cggbal_obj->iloref),
       cggbal_obj->ihi, cggbal_obj->ihiref,
       (cggbal_obj->ihi - cggbal_obj->ihiref) );
#endif
}

TEST_F(cggbal_test, cggbal1) {
    EXPECT_NEAR(0.0, cggbal_obj->diff_a, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_b, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_lscale, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_rscale, cggbal_obj->threshold);
    EXPECT_EQ(cggbal_obj->ilo, cggbal_obj->iloref);
    EXPECT_EQ(cggbal_obj->ihi, cggbal_obj->ihiref);
}

TEST_F(cggbal_test, cggbal2) {
    EXPECT_NEAR(0.0, cggbal_obj->diff_a, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_b, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_lscale, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_rscale, cggbal_obj->threshold);
    EXPECT_EQ(cggbal_obj->ilo, cggbal_obj->iloref);
    EXPECT_EQ(cggbal_obj->ihi, cggbal_obj->ihiref);
}

TEST_F(cggbal_test, cggbal3) {
    EXPECT_NEAR(0.0, cggbal_obj->diff_a, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_b, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_lscale, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_rscale, cggbal_obj->threshold);
    EXPECT_EQ(cggbal_obj->ilo, cggbal_obj->iloref);
    EXPECT_EQ(cggbal_obj->ihi, cggbal_obj->ihiref);
}

TEST_F(cggbal_test, cggbal4) {
    EXPECT_NEAR(0.0, cggbal_obj->diff_a, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_b, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_lscale, cggbal_obj->threshold);
    EXPECT_NEAR(0.0, cggbal_obj->diff_rscale, cggbal_obj->threshold);
    EXPECT_EQ(cggbal_obj->ilo, cggbal_obj->iloref);
    EXPECT_EQ(cggbal_obj->ihi, cggbal_obj->ihiref);
}

/* Begin ggbal_dcomplex_common_parameters  class definition */
class ggbal_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b, diff_lscale, diff_rscale;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'N' or 'P' or 'S' or 'B'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b

    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_double* b, *bref; // contains the n-by-n upper triangular matrix B.
    double *lscale, *lscaleref;
    double *rscale, *rscaleref;

    /* Input / Output parameters */
    lapack_int ilo, iloref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi, ihiref;

    /*Return Values */
    int info, inforef;

   public:
      ggbal_dcomplex_parameters (int matrix_layout_i, char job, lapack_int n);
      ~ggbal_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
ggbal_dcomplex_parameters:: ggbal_dcomplex_parameters (int matrix_layout_i,
                                          char job_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
    n = n_i;

    lda = n;
    ldb = n;

    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_b = 0;
    diff_lscale = 0;
    diff_rscale = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
    iloref = 0;
    ihiref = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbal lapack_complex_double: matrix_layout: %d n: %d  job: %c \n",
                                 matrix_layout, n, job);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_double_buffer_pair( &rscale, &rscaleref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) ){
       EXPECT_FALSE( true) << "ggbal_dcomplex_parameters object: malloc error. Exiting ";
       ggbal_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'ggbal_dcomplex_common_parameters' */
ggbal_dcomplex_parameters :: ~ggbal_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggbal_free();
} 

//  Test fixture class definition
class zggbal_test  : public  ::testing::Test {
public:
   ggbal_dcomplex_parameters  *zggbal_obj;
   void SetUp();  
   void TearDown () { delete zggbal_obj; }
};

void zggbal_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_zggbal) ( int matrix_layout, char job,
                    lapack_int n, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b,
                    lapack_int ldb, lapack_int* ilo, lapack_int* ihi,
                    double* lscale, double* rscale);
                 
    Fptr_NL_LAPACKE_zggbal ZGGBAL;

    zggbal_obj = new  ggbal_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].n );
                                         
    idx = Circular_Increment_Index(idx);
    zggbal_obj->threshold = eig_non_sym_paramslist[idx].ggbal_threshold;
    zggbal_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zggbal_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zggbal_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zggbal_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGBAL = (Fptr_NL_LAPACKE_zggbal)dlsym(zggbal_obj->hModule, "LAPACKE_zggbal");
    ASSERT_TRUE(ZGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_zggbal symbol";


    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zggbal_obj->inforef = ZGGBAL(   zggbal_obj->matrix_layout,
                                    zggbal_obj->job,
                                    zggbal_obj->n,
                                    zggbal_obj->aref,
                                    zggbal_obj->lda,
                                    zggbal_obj->bref,
                                    zggbal_obj->ldb,
                                    &zggbal_obj->iloref,
                                    &zggbal_obj->ihiref,
                                    zggbal_obj->lscaleref,
                                    zggbal_obj->rscaleref
                                    );

    /* Compute libflame's Lapacke o/p  */
    zggbal_obj->info = LAPACKE_zggbal(  zggbal_obj->matrix_layout,
                                    zggbal_obj->job,
                                    zggbal_obj->n,
                                    zggbal_obj->a,
                                    zggbal_obj->lda,
                                    zggbal_obj->b,
                                    zggbal_obj->ldb,
                                    &zggbal_obj->ilo,
                                    &zggbal_obj->ihi,
                                    zggbal_obj->lscale,
                                    zggbal_obj->rscale
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zggbal_obj->diff_a =  computeDiff_z( (zggbal_obj->lda)*(zggbal_obj->n), 
                zggbal_obj->a, zggbal_obj->aref );

    zggbal_obj->diff_b =  computeDiff_z( (zggbal_obj->ldb)*(zggbal_obj->n), 
                zggbal_obj->b, zggbal_obj->bref );

    zggbal_obj->diff_lscale =  computeDiff_d( zggbal_obj->n, 
                zggbal_obj->lscale, zggbal_obj->lscaleref );

    zggbal_obj->diff_rscale =  computeDiff_d( zggbal_obj->n, 
                zggbal_obj->rscale, zggbal_obj->rscaleref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbal lapack_complex_double: \n diff_a: %f \n diff_b: %f \n \
diff_lscale: %f \n diff_rscale: %f \n ilo: %d \n iloref: %d \
diff_ilo: %d \n ihi: %d \n ihiref: %d \n diff_ihi: %d \n",
       zggbal_obj->diff_a, zggbal_obj->diff_b, zggbal_obj->diff_lscale,
       zggbal_obj->diff_rscale, zggbal_obj->ilo, zggbal_obj->iloref,
       (zggbal_obj->ilo - zggbal_obj->iloref),
       zggbal_obj->ihi, zggbal_obj->ihiref,
       (zggbal_obj->ihi - zggbal_obj->ihiref) );
#endif
}

TEST_F(zggbal_test, zggbal1) {
    EXPECT_NEAR(0.0, zggbal_obj->diff_a, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_b, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_lscale, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_rscale, zggbal_obj->threshold);
    EXPECT_EQ(zggbal_obj->ilo, zggbal_obj->iloref);
    EXPECT_EQ(zggbal_obj->ihi, zggbal_obj->ihiref);
}

TEST_F(zggbal_test, zggbal2) {
    EXPECT_NEAR(0.0, zggbal_obj->diff_a, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_b, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_lscale, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_rscale, zggbal_obj->threshold);
    EXPECT_EQ(zggbal_obj->ilo, zggbal_obj->iloref);
    EXPECT_EQ(zggbal_obj->ihi, zggbal_obj->ihiref);
}

TEST_F(zggbal_test, zggbal3) {
    EXPECT_NEAR(0.0, zggbal_obj->diff_a, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_b, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_lscale, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_rscale, zggbal_obj->threshold);
    EXPECT_EQ(zggbal_obj->ilo, zggbal_obj->iloref);
    EXPECT_EQ(zggbal_obj->ihi, zggbal_obj->ihiref);
}

TEST_F(zggbal_test, zggbal4) {
    EXPECT_NEAR(0.0, zggbal_obj->diff_a, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_b, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_lscale, zggbal_obj->threshold);
    EXPECT_NEAR(0.0, zggbal_obj->diff_rscale, zggbal_obj->threshold);
    EXPECT_EQ(zggbal_obj->ilo, zggbal_obj->iloref);
    EXPECT_EQ(zggbal_obj->ihi, zggbal_obj->ihiref);
}