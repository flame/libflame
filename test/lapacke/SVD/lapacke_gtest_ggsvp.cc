#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define ggsvp_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (u!=NULL)        free(u); \
    if (uref!=NULL)     free(uref); \
    if (v!=NULL)        free(v); \
    if (vref!=NULL)     free(vref); \
    if (q!=NULL)        free(q); \
    if (qref!=NULL)     free(qref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ggsvp_float_common_parameters  class definition */
class ggsvp_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_u, diff_v, diff_q;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char jobu; // Must be 'U' or 'N'.
    char jobv; // Must be 'V' or 'N'.
    char jobq; // Must be 'Q' or 'N'.
    lapack_int m; // The number of rows of the matrix A
    lapack_int p; // The number of rows of the matrix B
    lapack_int n; // The number of columns of the matrices A and B
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b

    // thresholds to determine the effective numerical rank of matrix B and a subblock of A.
    float tola;
    float tolb;
    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p) 
    lapack_int ldq; // The leading dimension of the output array q . ldq≥ max(1, n)
    /* Input / Output parameters */
    float* a, *aref; // contains m-by-n matrix A.
    float* b, *bref; // contains  p-by-n matrix B.

    /* Output parameters */
    lapack_int k, kref; // the dimension of subblocks.
    lapack_int l, lref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    float* u, *uref;
    float* v, *vref;
    float* q, *qref;
    /*Return Values */
    int info, inforef;

   public:
      ggsvp_float_parameters (int matrix_layout_i, char jobu, char jobv,
                                  char jobq, lapack_int m, lapack_int p, 
                                  lapack_int n, float tola, float tolb);
      ~ggsvp_float_parameters ();
};

/* Constructor definition  float_common_parameters */
ggsvp_float_parameters:: ggsvp_float_parameters (int matrix_layout_i,
                    char jobu_i, char jobv_i, char jobq_i, lapack_int m_i,  
            lapack_int p_i, lapack_int n_i, float tola_i, float tolb_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    jobq = jobq_i;
    n = n_i;
    m = m_i;
    p = p_i;
    tola = tola_i;
    tolb = tolb_i;

    if (matrix_layout == LAPACK_ROW_MAJOR)
    {
        lda = n;
        ldb = n;
    }   else
    {
        lda = m;
        ldb = p;
    }
    ldq = n;
    ldu = m;
    ldv = p;

    // Initialize with default values
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_b = 0;
    diff_u = 0;
    diff_v = 0;
    diff_q = 0;
    k = kref = 0;
    l = lref = 0;


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggsvp float: matrix_layout: %d   jobu: %c \t \
jobv: %c \t jobq: %c  \n m: %d \t p: %d \t n: %d \n \
tola: %f  \t tolb: %f \n", matrix_layout, jobu_i, jobv_i, jobq_i, m_i,  
p_i, n_i, tola_i, tolb_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, p*n );
    lapacke_gtest_alloc_float_buffer_pair( &u, &uref, ldu*m );
    lapacke_gtest_alloc_float_buffer_pair( &v, &vref, ldv*p );
    lapacke_gtest_alloc_float_buffer_pair( &q, &qref, ldq*n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (v==NULL) || (vref==NULL) || \
        (q==NULL) || (qref==NULL)  ){
       EXPECT_FALSE( true) << "ggsvp_float_parameters object: malloc error. Exiting ";
       ggsvp_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, p*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( u, uref, ldu*m, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( v, vref, ldv*p, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( q, qref, ldq*n, 0);

   } /* end of Constructor  */

/* Destructor definition  'ggsvp_float_common_parameters' */
ggsvp_float_parameters :: ~ggsvp_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggsvp_free();
} 

//  Test fixture class definition
class sggsvp_test  : public  ::testing::Test {
public:
   ggsvp_float_parameters  *sggsvp_obj;
   void SetUp();  
   void TearDown () { delete sggsvp_obj; }
};

void sggsvp_test::SetUp()
{
    /* LAPACKE sggsvp prototype */
    typedef int (*Fptr_NL_LAPACKE_sggsvp) ( int matrix_layout, char jobu, 
            char jobv, char jobq, lapack_int m, lapack_int p, lapack_int n, 
            float* a, lapack_int lda, float* b, lapack_int ldb, float tola,
            float tolb, lapack_int* k, lapack_int* l, float* u, lapack_int ldu,
            float* v, lapack_int ldv, float* q, lapack_int ldqe);

    Fptr_NL_LAPACKE_sggsvp SGGSVP;

    sggsvp_obj = new  ggsvp_float_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu,
                                         svd_paramslist[idx].jobv,
                                         svd_paramslist[idx].jobq,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].p,
                                         svd_paramslist[idx].n,
                                         svd_paramslist[idx].tola,
                                         svd_paramslist[idx].tolb );

    idx = Circular_Increment_Index(idx);
    sggsvp_obj->threshold = svd_paramslist[idx].svd_threshold;
    sggsvp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sggsvp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sggsvp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sggsvp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGSVP = (Fptr_NL_LAPACKE_sggsvp)dlsym(sggsvp_obj->hModule, "LAPACKE_sggsvp3");
    ASSERT_TRUE(SGGSVP != NULL) << "failed to ppt the Netlib LAPACKE_sggsvp symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sggsvp_obj->inforef = SGGSVP(   sggsvp_obj->matrix_layout,
                                    sggsvp_obj->jobu,
                                    sggsvp_obj->jobv,
                                    sggsvp_obj->jobq,
                                    sggsvp_obj->m,
                                    sggsvp_obj->p,
                                    sggsvp_obj->n,
                                    sggsvp_obj->aref,
                                    sggsvp_obj->lda,
                                    sggsvp_obj->bref,
                                    sggsvp_obj->ldb,
                                    sggsvp_obj->tola,
                                    sggsvp_obj->tolb,
                                    &sggsvp_obj->kref,
                                    &sggsvp_obj->lref,
                                    sggsvp_obj->uref,
                                    sggsvp_obj->ldu,
                                    sggsvp_obj->vref,
                                    sggsvp_obj->ldv,
                                    sggsvp_obj->qref,
                                    sggsvp_obj->ldq
                                    );

    /* Compute libflame's Lapacke o/p  */
    sggsvp_obj->info = LAPACKE_sggsvp3(  sggsvp_obj->matrix_layout,
                                    sggsvp_obj->jobu,
                                    sggsvp_obj->jobv,
                                    sggsvp_obj->jobq,
                                    sggsvp_obj->m,
                                    sggsvp_obj->p,
                                    sggsvp_obj->n,
                                    sggsvp_obj->a,
                                    sggsvp_obj->lda,
                                    sggsvp_obj->b,
                                    sggsvp_obj->ldb,
                                    sggsvp_obj->tola,
                                    sggsvp_obj->tolb,
                                    &sggsvp_obj->k,
                                    &sggsvp_obj->l,
                                    sggsvp_obj->u,
                                    sggsvp_obj->ldu,
                                    sggsvp_obj->v,
                                    sggsvp_obj->ldv,
                                    sggsvp_obj->q,
                                    sggsvp_obj->ldq
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    sggsvp_obj->diff_a =  computeDiff_s( (sggsvp_obj->m)*(sggsvp_obj->n), 
                sggsvp_obj->a, sggsvp_obj->aref );

    sggsvp_obj->diff_b =  computeDiff_s( (sggsvp_obj->p)*(sggsvp_obj->n), 
                sggsvp_obj->b, sggsvp_obj->bref );

    sggsvp_obj->diff_u =  computeDiff_s( (sggsvp_obj->ldu)*(sggsvp_obj->m), 
                sggsvp_obj->u, sggsvp_obj->uref );

    sggsvp_obj->diff_v =  computeDiff_s( (sggsvp_obj->p)*(sggsvp_obj->ldv), 
                sggsvp_obj->v, sggsvp_obj->vref );

    sggsvp_obj->diff_q =  computeDiff_s( (sggsvp_obj->ldq)*(sggsvp_obj->n), 
                sggsvp_obj->q, sggsvp_obj->qref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggsvp float: \n diff_a: %f \n diff_b: %f \n \
diff_u: %f \n diff_v: %f \n k: %d \n kref: %d \
diff_q: %f \n l: %d \n lref: %d \n ",
       sggsvp_obj->diff_a, sggsvp_obj->diff_b, sggsvp_obj->diff_u,
       sggsvp_obj->diff_v, sggsvp_obj->k, sggsvp_obj->kref, sggsvp_obj->diff_q,
       sggsvp_obj->l, sggsvp_obj->lref );
#endif
}

TEST_F(sggsvp_test, sggsvp1) {
    EXPECT_NEAR(0.0, sggsvp_obj->diff_a, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_b, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_q, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_u, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_v, sggsvp_obj->threshold);
    EXPECT_EQ(sggsvp_obj->k, sggsvp_obj->kref);
    EXPECT_EQ(sggsvp_obj->l, sggsvp_obj->lref);
}

TEST_F(sggsvp_test, sggsvp2) {
    EXPECT_NEAR(0.0, sggsvp_obj->diff_a, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_b, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_q, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_u, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_v, sggsvp_obj->threshold);
    EXPECT_EQ(sggsvp_obj->k, sggsvp_obj->kref);
    EXPECT_EQ(sggsvp_obj->l, sggsvp_obj->lref);
}

TEST_F(sggsvp_test, sggsvp3) {
    EXPECT_NEAR(0.0, sggsvp_obj->diff_a, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_b, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_q, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_u, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_v, sggsvp_obj->threshold);
    EXPECT_EQ(sggsvp_obj->k, sggsvp_obj->kref);
    EXPECT_EQ(sggsvp_obj->l, sggsvp_obj->lref);
}

TEST_F(sggsvp_test, sggsvp4) {
    EXPECT_NEAR(0.0, sggsvp_obj->diff_a, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_b, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_q, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_u, sggsvp_obj->threshold);
    EXPECT_NEAR(0.0, sggsvp_obj->diff_v, sggsvp_obj->threshold);
    EXPECT_EQ(sggsvp_obj->k, sggsvp_obj->kref);
    EXPECT_EQ(sggsvp_obj->l, sggsvp_obj->lref);
}
