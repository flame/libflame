#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define gghd3_free() \
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

/* Begin gghd3_float_common_parameters  class definition */
class gghd3_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_lscale, diff_rscale;
    float threshold;
    void *hModule, *dModule;
	
lapack_int LAPACKE_sgghd3 (int matrix_layout, char compq, char compz, lapack_int n, lapack_int ilo, lapack_int ihi, float * a, lapack_int lda, float * b, lapack_int ldb, float * q, lapack_int ldq, float * z, lapack_int ldz);
    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'N' or 'P' or 'S' or 'B'.
    char compq; // Must be 'N', 'I', or 'V'.
    char compz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b

    /* Input parameters */
    lapack_int ilo, iloref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi, ihiref;
    float *lscale, *lscaleref;
    float *rscale, *rscaleref;
	
    /* Input / Output parameters */
    float* a, *aref; // contains the n-by-n general matrix A.
    float* b, *bref; // contains the n-by-n upper triangular matrix B.
    float* q, *qref; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    float* z, *zref; //


    /*Return Values */
    int info, inforef;

   public:
      gghd3_float_parameters int matrix_layout_i, char compq, char compz,
                              lapack_int n);
      ~gghd3_float_parameters ();
};

/* Constructor definition  float_common_parameters */
gghd3_float_parameters:: gghd3_float_parameters (int matrix_layout_i,
                            char compq_i, char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    compq = compq_i;
    compz = compz_i;
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
   printf(" \n gghd3 float: matrix_layout: %d n: %d  job: %c \n",
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
       EXPECT_FALSE( true) << "gghd3_float_parameters object: malloc error. Exiting ";
       gghd3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    //lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'gghd3_float_common_parameters' */
gghd3_float_parameters :: ~gghd3_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gghd3_free();
} 

//  Test fixture class definition
class sgghd3_test  : public  ::testing::Test {
public:
   gghd3_float_parameters  *sgghd3_obj;
   void SetUp();  
   void TearDown () { delete sgghd3_obj; }
};

void sgghd3_test::SetUp()
{
    /* LAPACKE SGGBAL prototype */
    typedef int (*Fptr_NL_LAPACKE_sggbal) ( int matrix_layout, char job,
                    lapack_int n, float* a, lapack_int lda, float* b,
                    lapack_int ldb, lapack_int* ilo, lapack_int* ihi,
                    float* lscale, float* rscale);
                 
    Fptr_NL_LAPACKE_sggbal SGGBAL;

    sgghd3_obj = new  gghd3_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].n );
                                         
    idx = Circular_Increment_Index(idx);
    sgghd3_obj->threshold = eig_non_sym_paramslist[idx].gghd3_threshold;
    sgghd3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgghd3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgghd3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgghd3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGBAL = (Fptr_NL_LAPACKE_sggbal)dlsym(sgghd3_obj->hModule, "LAPACKE_sggbal");
    ASSERT_TRUE(SGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_sggbal symbol";


    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgghd3_obj->inforef = SGGBAL(   sgghd3_obj->matrix_layout,
                                    sgghd3_obj->job,
                                    sgghd3_obj->n,
                                    sgghd3_obj->aref,
                                    sgghd3_obj->lda,
                                    sgghd3_obj->bref,
                                    sgghd3_obj->ldb,
                                    &sgghd3_obj->iloref,
                                    &sgghd3_obj->ihiref,
                                    sgghd3_obj->lscaleref,
                                    sgghd3_obj->rscaleref
                                    );

    /* Compute libflame's Lapacke o/p  */
    sgghd3_obj->info = LAPACKE_sggbal(  sgghd3_obj->matrix_layout,
                                    sgghd3_obj->job,
                                    sgghd3_obj->n,
                                    sgghd3_obj->a,
                                    sgghd3_obj->lda,
                                    sgghd3_obj->b,
                                    sgghd3_obj->ldb,
                                    &sgghd3_obj->ilo,
                                    &sgghd3_obj->ihi,
                                    sgghd3_obj->lscale,
                                    sgghd3_obj->rscale
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sgghd3_obj->diff_a =  computeDiff_s( (sgghd3_obj->lda)*(sgghd3_obj->n), 
                sgghd3_obj->a, sgghd3_obj->aref );

    sgghd3_obj->diff_b =  computeDiff_s( (sgghd3_obj->ldb)*(sgghd3_obj->n), 
                sgghd3_obj->b, sgghd3_obj->bref );

    sgghd3_obj->diff_lscale =  computeDiff_s( sgghd3_obj->n, 
                sgghd3_obj->lscale, sgghd3_obj->lscaleref );

    sgghd3_obj->diff_rscale =  computeDiff_s( sgghd3_obj->n, 
                sgghd3_obj->rscale, sgghd3_obj->rscaleref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gghd3 float: \n diff_a: %f \n diff_b: %f \n \
diff_lscale: %f \n diff_rscale: %f \n ilo: %d \n iloref: %d \
diff_ilo: %d \n ihi: %d \n ihiref: %d \n diff_ihi: %d \n",
       sgghd3_obj->diff_a, sgghd3_obj->diff_b, sgghd3_obj->diff_lscale,
       sgghd3_obj->diff_rscale, sgghd3_obj->ilo, sgghd3_obj->iloref,
       (sgghd3_obj->ilo - sgghd3_obj->iloref),
       sgghd3_obj->ihi, sgghd3_obj->ihiref,
       (sgghd3_obj->ihi - sgghd3_obj->ihiref) );
#endif
}

TEST_F(sgghd3_test, sgghd31) {
    EXPECT_NEAR(0.0, sgghd3_obj->diff_a, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_b, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_lscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_rscale, sgghd3_obj->threshold);
    EXPECT_EQ(sgghd3_obj->ilo, sgghd3_obj->iloref);
    EXPECT_EQ(sgghd3_obj->ihi, sgghd3_obj->ihiref);
}

TEST_F(sgghd3_test, sgghd32) {
    EXPECT_NEAR(0.0, sgghd3_obj->diff_a, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_b, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_lscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_rscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ilo - sgghd3_obj->iloref), sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ihi - sgghd3_obj->ihiref), sgghd3_obj->threshold);
}

TEST_F(sgghd3_test, sgghd33) {
    EXPECT_NEAR(0.0, sgghd3_obj->diff_a, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_b, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_lscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_rscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ilo - sgghd3_obj->iloref), sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ihi - sgghd3_obj->ihiref), sgghd3_obj->threshold);
}

TEST_F(sgghd3_test, sgghd34) {
    EXPECT_NEAR(0.0, sgghd3_obj->diff_a, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_b, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_lscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_rscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ilo - sgghd3_obj->iloref), sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ihi - sgghd3_obj->ihiref), sgghd3_obj->threshold);
}

    if (rscaleref!=NULL)     free(rscaleref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin gghd3_float_common_parameters  class definition */
class gghd3_float_parameters{

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
      gghd3_float_parameters (int matrix_layout_i, char job, lapack_int n);
      ~gghd3_float_parameters ();
};

/* Constructor definition  float_common_parameters */
gghd3_float_parameters:: gghd3_float_parameters (int matrix_layout_i,
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
   printf(" \n gghd3 float: matrix_layout: %d n: %d  job: %c \n",
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
       EXPECT_FALSE( true) << "gghd3_float_parameters object: malloc error. Exiting ";
       gghd3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    //lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'gghd3_float_common_parameters' */
gghd3_float_parameters :: ~gghd3_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   gghd3_free();
} 

//  Test fixture class definition
class sgghd3_test  : public  ::testing::Test {
public:
   gghd3_float_parameters  *sgghd3_obj;
   void SetUp();  
   void TearDown () { delete sgghd3_obj; }
};

void sgghd3_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_sgghd3) ( int matrix_layout, char job,
                    lapack_int n, float* a, lapack_int lda, float* b,
                    lapack_int ldb, lapack_int* ilo, lapack_int* ihi,
                    float* lscale, float* rscale);
                 
    Fptr_NL_LAPACKE_sgghd3 SGGHD3;

    sgghd3_obj = new  gghd3_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].n );
                                         
    idx = Circular_Increment_Index(idx);
    sgghd3_obj->threshold = eig_non_sym_paramslist[idx].gghd3_threshold;
    sgghd3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgghd3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgghd3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgghd3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGHD3 = (Fptr_NL_LAPACKE_sgghd3)dlsym(sgghd3_obj->hModule, "LAPACKE_sgghd3");
    ASSERT_TRUE(SGGHD3 != NULL) << "failed to ppt the Netlib LAPACKE_sgghd3 symbol";


    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sgghd3_obj->inforef = SGGHD3(   sgghd3_obj->matrix_layout,
                                    sgghd3_obj->job,
                                    sgghd3_obj->n,
                                    sgghd3_obj->aref,
                                    sgghd3_obj->lda,
                                    sgghd3_obj->bref,
                                    sgghd3_obj->ldb,
                                    &sgghd3_obj->iloref,
                                    &sgghd3_obj->ihiref,
                                    sgghd3_obj->lscaleref,
                                    sgghd3_obj->rscaleref
                                    );

    /* Compute libflame's Lapacke o/p  */
    sgghd3_obj->info = LAPACKE_sgghd3(  sgghd3_obj->matrix_layout,
                                    sgghd3_obj->job,
                                    sgghd3_obj->n,
                                    sgghd3_obj->a,
                                    sgghd3_obj->lda,
                                    sgghd3_obj->b,
                                    sgghd3_obj->ldb,
                                    &sgghd3_obj->ilo,
                                    &sgghd3_obj->ihi,
                                    sgghd3_obj->lscale,
                                    sgghd3_obj->rscale
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sgghd3_obj->diff_a =  computeDiff_s( (sgghd3_obj->lda)*(sgghd3_obj->n), 
                sgghd3_obj->a, sgghd3_obj->aref );

    sgghd3_obj->diff_b =  computeDiff_s( (sgghd3_obj->ldb)*(sgghd3_obj->n), 
                sgghd3_obj->b, sgghd3_obj->bref );

    sgghd3_obj->diff_lscale =  computeDiff_s( sgghd3_obj->n, 
                sgghd3_obj->lscale, sgghd3_obj->lscaleref );

    sgghd3_obj->diff_rscale =  computeDiff_s( sgghd3_obj->n, 
                sgghd3_obj->rscale, sgghd3_obj->rscaleref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n gghd3 float: \n diff_a: %f \n diff_b: %f \n \
diff_lscale: %f \n diff_rscale: %f \n ilo: %d \n iloref: %d \
diff_ilo: %d \n ihi: %d \n ihiref: %d \n diff_ihi: %d \n",
       sgghd3_obj->diff_a, sgghd3_obj->diff_b, sgghd3_obj->diff_lscale,
       sgghd3_obj->diff_rscale, sgghd3_obj->ilo, sgghd3_obj->iloref,
       (sgghd3_obj->ilo - sgghd3_obj->iloref),
       sgghd3_obj->ihi, sgghd3_obj->ihiref,
       (sgghd3_obj->ihi - sgghd3_obj->ihiref) );
#endif
}

TEST_F(sgghd3_test, sgghd31) {
    EXPECT_NEAR(0.0, sgghd3_obj->diff_a, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_b, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_lscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_rscale, sgghd3_obj->threshold);
    EXPECT_EQ(sgghd3_obj->ilo, sgghd3_obj->iloref);
    EXPECT_EQ(sgghd3_obj->ihi, sgghd3_obj->ihiref);
}

TEST_F(sgghd3_test, sgghd32) {
    EXPECT_NEAR(0.0, sgghd3_obj->diff_a, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_b, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_lscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_rscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ilo - sgghd3_obj->iloref), sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ihi - sgghd3_obj->ihiref), sgghd3_obj->threshold);
}

TEST_F(sgghd3_test, sgghd33) {
    EXPECT_NEAR(0.0, sgghd3_obj->diff_a, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_b, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_lscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_rscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ilo - sgghd3_obj->iloref), sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ihi - sgghd3_obj->ihiref), sgghd3_obj->threshold);
}

TEST_F(sgghd3_test, sgghd34) {
    EXPECT_NEAR(0.0, sgghd3_obj->diff_a, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_b, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_lscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, sgghd3_obj->diff_rscale, sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ilo - sgghd3_obj->iloref), sgghd3_obj->threshold);
    EXPECT_NEAR(0.0, (sgghd3_obj->ihi - sgghd3_obj->ihiref), sgghd3_obj->threshold);
}
