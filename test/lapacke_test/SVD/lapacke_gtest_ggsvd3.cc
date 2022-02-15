#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define ggsvd3_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (alpha!=NULL)        free(alpha); \
    if (alpharef!=NULL)     free(alpharef); \
    if (beta!=NULL)        free(beta); \
    if (betaref!=NULL)     free(betaref); \
    if (iwork!=NULL)        free(iwork); \
    if (iworkref!=NULL)     free(iworkref); \
    if (u!=NULL)        free(u); \
    if (uref!=NULL)     free(uref); \
    if (v!=NULL)        free(v); \
    if (vref!=NULL)     free(vref); \
    if (q!=NULL)        free(q); \
    if (qref!=NULL)     free(qref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ggsvd3_float_common_parameters  class definition */
class ggsvd3_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_u, diff_v, diff_q, diff_alpha, diff_beta;
	int diff_iwork;
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
	
	/*  generalized singular value pairs of a and b; */
    float* alpha, *alpharef;
    float* beta, *betaref;
	
	/*  stores the sorting information.  */
	lapack_int *iwork, *iworkref;
    /*Return Values */
    int info, inforef;

   public:
      ggsvd3_float_parameters (int matrix_layout_i, char jobu, char jobv,
                                  char jobq, lapack_int m, lapack_int p, 
                                  lapack_int n);
      ~ggsvd3_float_parameters ();
};

/* Constructor definition  ggsvd3 float_common_parameters */
ggsvd3_float_parameters:: ggsvd3_float_parameters (int matrix_layout_i,
                    char jobu_i, char jobv_i, char jobq_i, lapack_int m_i,  
            lapack_int p_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    jobq = jobq_i;
    n = n_i;
    m = m_i;
    p = p_i;

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
	diff_alpha =0;
	diff_beta = 0;
	diff_iwork = 0;
    k = kref = 0;
    l = lref = 0;


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggsvd3 float: matrix_layout: %d   jobu: %c \t \
jobv: %c \t jobq: %c  \n m: %d \t p: %d \t n: %d \n",
matrix_layout, jobu_i, jobv_i, jobq_i, m_i, p_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_int_buffer_pair( &iwork, &iworkref, n );
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, p*n );
    lapacke_gtest_alloc_float_buffer_pair( &u, &uref, ldu*m );
    lapacke_gtest_alloc_float_buffer_pair( &v, &vref, ldv*p );
    lapacke_gtest_alloc_float_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_float_buffer_pair( &alpha, &alpharef, n );
    lapacke_gtest_alloc_float_buffer_pair( &beta, &betaref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (v==NULL) || (vref==NULL) || \
        (iwork==NULL) || (iworkref==NULL) || \
        (alpha==NULL) || (alpharef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (q==NULL) || (qref==NULL)  ){
       EXPECT_FALSE( true) << "ggsvd3_float_parameters object: malloc error. Exiting ";
       ggsvd3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, p*n );
    lapacke_gtest_init_float_buffer_pair_with_constant( u, uref, ldu*m, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( v, vref, ldv*p, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( q, qref, ldq*n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( alpha, alpharef, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant( iwork, iworkref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'ggsvd3_float_common_parameters' */
ggsvd3_float_parameters :: ~ggsvd3_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggsvd3_free();
} 

//  Test fixture class definition
class sggsvd3_test  : public  ::testing::Test {
public:
   ggsvd3_float_parameters  *sggsvd3_obj;
   void SetUp();  
   void TearDown () { delete sggsvd3_obj; }
};

void sggsvd3_test::SetUp()
{
    /* LAPACKE STRSNA prototype */
    typedef int (*Fptr_NL_LAPACKE_sggsvd3) ( int matrix_layout, char jobu, 
			char jobv, char jobq, lapack_int m, lapack_int n, lapack_int p, 
			lapack_int * k, lapack_int * l, float * a, lapack_int lda, 
			float * b, lapack_int ldb, float * alpha, float * beta, float * u, 
			lapack_int ldu, float * v, lapack_int ldv, float * q, 
			lapack_int ldq, lapack_int * iwork);
			
    Fptr_NL_LAPACKE_sggsvd3 SGGSVD3;

    sggsvd3_obj = new  ggsvd3_float_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu,
                                         svd_paramslist[idx].jobv,
                                         svd_paramslist[idx].jobq,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].p,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    sggsvd3_obj->threshold = svd_paramslist[idx].svd_threshold;
    sggsvd3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sggsvd3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sggsvd3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sggsvd3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGSVD3 = (Fptr_NL_LAPACKE_sggsvd3)dlsym(sggsvd3_obj->hModule, "LAPACKE_sggsvd3");
    ASSERT_TRUE(SGGSVD3 != NULL) << "failed to ppt the Netlib LAPACKE_sggsvd3 symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sggsvd3_obj->inforef = SGGSVD3( sggsvd3_obj->matrix_layout,
                                    sggsvd3_obj->jobu,
                                    sggsvd3_obj->jobv,
                                    sggsvd3_obj->jobq,
                                    sggsvd3_obj->m,
                                    sggsvd3_obj->n,
                                    sggsvd3_obj->p,
                                    &sggsvd3_obj->kref,
                                    &sggsvd3_obj->lref,
                                    sggsvd3_obj->aref,
                                    sggsvd3_obj->lda,
                                    sggsvd3_obj->bref,
                                    sggsvd3_obj->ldb,
                                    sggsvd3_obj->alpha,
                                    sggsvd3_obj->beta,
                                    sggsvd3_obj->uref,
                                    sggsvd3_obj->ldu,
                                    sggsvd3_obj->vref,
                                    sggsvd3_obj->ldv,
                                    sggsvd3_obj->qref,
                                    sggsvd3_obj->ldq,
                                    sggsvd3_obj->iworkref
                                    );

    /* Compute libflame's Lapacke o/p  */
    sggsvd3_obj->info = LAPACKE_sggsvd3(  sggsvd3_obj->matrix_layout,
                                    sggsvd3_obj->jobu,
                                    sggsvd3_obj->jobv,
                                    sggsvd3_obj->jobq,
                                    sggsvd3_obj->m,
                                    sggsvd3_obj->n,
                                    sggsvd3_obj->p,
                                    &sggsvd3_obj->k,
                                    &sggsvd3_obj->l,
                                    sggsvd3_obj->a,
                                    sggsvd3_obj->lda,
                                    sggsvd3_obj->b,
                                    sggsvd3_obj->ldb,
                                    sggsvd3_obj->alpha,
                                    sggsvd3_obj->beta,
                                    sggsvd3_obj->u,
                                    sggsvd3_obj->ldu,
                                    sggsvd3_obj->v,
                                    sggsvd3_obj->ldv,
                                    sggsvd3_obj->q,
                                    sggsvd3_obj->ldq,
                                    sggsvd3_obj->iwork
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    sggsvd3_obj->diff_a =  computeDiff_s( (sggsvd3_obj->m)*(sggsvd3_obj->n), 
                sggsvd3_obj->a, sggsvd3_obj->aref );

    sggsvd3_obj->diff_b =  computeDiff_s( (sggsvd3_obj->p)*(sggsvd3_obj->n), 
                sggsvd3_obj->b, sggsvd3_obj->bref );

    sggsvd3_obj->diff_u =  computeDiff_s( (sggsvd3_obj->ldu)*(sggsvd3_obj->m), 
                sggsvd3_obj->u, sggsvd3_obj->uref );

    sggsvd3_obj->diff_v =  computeDiff_s( (sggsvd3_obj->p)*(sggsvd3_obj->ldv), 
                sggsvd3_obj->v, sggsvd3_obj->vref );

    sggsvd3_obj->diff_q =  computeDiff_s( (sggsvd3_obj->ldq)*(sggsvd3_obj->n), 
                sggsvd3_obj->q, sggsvd3_obj->qref );

    sggsvd3_obj->diff_alpha =  computeDiff_s( sggsvd3_obj->n, 
                sggsvd3_obj->alpha, sggsvd3_obj->alpharef );

    sggsvd3_obj->diff_beta =  computeDiff_s( sggsvd3_obj->n, 
                sggsvd3_obj->beta, sggsvd3_obj->betaref );
				
    sggsvd3_obj->diff_iwork =  computeDiff_i( sggsvd3_obj->n, 
                sggsvd3_obj->iwork, sggsvd3_obj->iworkref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggsvd3 float: \n diff_a: %f \n diff_b: %f \n \
diff_u: %f \n diff_v: %f \n k: %d \n kref: %d \
diff_q: %f \n l: %d \n lref: %d \n  diff_alpha: %f \n \
diff_beta: %f \n diff_iwork: %d  \n ",
       sggsvd3_obj->diff_a, sggsvd3_obj->diff_b, sggsvd3_obj->diff_u,
       sggsvd3_obj->diff_v, sggsvd3_obj->k, sggsvd3_obj->kref, sggsvd3_obj->diff_q,
       sggsvd3_obj->l, sggsvd3_obj->lref, sggsvd3_obj->diff_alpha,
	   sggsvd3_obj->diff_beta,  sggsvd3_obj->diff_iwork);
#endif
}

TEST_F(sggsvd3_test, sggsvd3_1) {
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_a, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_b, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_q, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_u, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_v, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_alpha, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_beta, sggsvd3_obj->threshold);
    EXPECT_EQ(0, sggsvd3_obj->diff_iwork);
    EXPECT_EQ(sggsvd3_obj->k, sggsvd3_obj->kref);
    EXPECT_EQ(sggsvd3_obj->l, sggsvd3_obj->lref);
}

TEST_F(sggsvd3_test, sggsvd3_2) {
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_a, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_b, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_q, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_u, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_v, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_alpha, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_beta, sggsvd3_obj->threshold);
    EXPECT_EQ(0, sggsvd3_obj->diff_iwork);
    EXPECT_EQ(sggsvd3_obj->k, sggsvd3_obj->kref);
    EXPECT_EQ(sggsvd3_obj->l, sggsvd3_obj->lref);
}

TEST_F(sggsvd3_test, sggsvd3_3) {
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_a, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_b, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_q, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_u, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_v, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_alpha, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_beta, sggsvd3_obj->threshold);
    //EXPECT_EQ(0, sggsvd3_obj->diff_iwork);
    EXPECT_EQ(sggsvd3_obj->k, sggsvd3_obj->kref);
    EXPECT_EQ(sggsvd3_obj->l, sggsvd3_obj->lref);
}

TEST_F(sggsvd3_test, sggsvd3_4) {
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_a, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_b, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_q, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_u, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_v, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_alpha, sggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, sggsvd3_obj->diff_beta, sggsvd3_obj->threshold);
    //EXPECT_EQ(0, sggsvd3_obj->diff_iwork);
    EXPECT_EQ(sggsvd3_obj->k, sggsvd3_obj->kref);
    EXPECT_EQ(sggsvd3_obj->l, sggsvd3_obj->lref);
}


/* Begin ggsvd3_double_common_parameters  class definition */
class ggsvd3_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b, diff_u, diff_v, diff_q, diff_alpha, diff_beta;
	int diff_iwork;
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

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p) 
    lapack_int ldq; // The leading dimension of the output array q . ldq≥ max(1, n)
    /* Input / Output parameters */
    double* a, *aref; // contains m-by-n matrix A.
    double* b, *bref; // contains  p-by-n matrix B.

    /* Output parameters */
    lapack_int k, kref; // the dimension of subblocks.
    lapack_int l, lref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    double* u, *uref;
    double* v, *vref;
    double* q, *qref;
	
	/*  generalized singular value pairs of a and b; */
    double* alpha, *alpharef;
    double* beta, *betaref;
	
	/*  stores the sorting information.  */
	lapack_int *iwork, *iworkref;
    /*Return Values */
    int info, inforef;

   public:
      ggsvd3_double_parameters (int matrix_layout_i, char jobu, char jobv,
                                  char jobq, lapack_int m, lapack_int p, 
                                  lapack_int n);
      ~ggsvd3_double_parameters ();
};

/* Constructor definition  ggsvd3 double_common_parameters */
ggsvd3_double_parameters:: ggsvd3_double_parameters (int matrix_layout_i,
                    char jobu_i, char jobv_i, char jobq_i, lapack_int m_i,  
            lapack_int p_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    jobq = jobq_i;
    n = n_i;
    m = m_i;
    p = p_i;

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
	diff_alpha =0;
	diff_beta = 0;
	diff_iwork = 0;
    k = kref = 0;
    l = lref = 0;


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggsvd3 double: matrix_layout: %d   jobu: %c \t \
jobv: %c \t jobq: %c  \n m: %d \t p: %d \t n: %d \n",
matrix_layout, jobu_i, jobv_i, jobq_i, m_i, p_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_int_buffer_pair( &iwork, &iworkref, n );
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, p*n );
    lapacke_gtest_alloc_double_buffer_pair( &u, &uref, ldu*m );
    lapacke_gtest_alloc_double_buffer_pair( &v, &vref, ldv*p );
    lapacke_gtest_alloc_double_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_double_buffer_pair( &alpha, &alpharef, n );
    lapacke_gtest_alloc_double_buffer_pair( &beta, &betaref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (v==NULL) || (vref==NULL) || \
        (iwork==NULL) || (iworkref==NULL) || \
        (alpha==NULL) || (alpharef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (q==NULL) || (qref==NULL)  ){
       EXPECT_FALSE( true) << "ggsvd3_double_parameters object: malloc error. Exiting ";
       ggsvd3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, p*n );
    lapacke_gtest_init_double_buffer_pair_with_constant( u, uref, ldu*m, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( v, vref, ldv*p, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( q, qref, ldq*n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( alpha, alpharef, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant( iwork, iworkref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'ggsvd3_double_common_parameters' */
ggsvd3_double_parameters :: ~ggsvd3_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggsvd3_free();
} 

//  Test fixture class definition
class dggsvd3_test  : public  ::testing::Test {
public:
   ggsvd3_double_parameters  *dggsvd3_obj;
   void SetUp();  
   void TearDown () { delete dggsvd3_obj; }
};

void dggsvd3_test::SetUp()
{
    /* LAPACKE_dggsvd3 prototype */
    typedef int (*Fptr_NL_LAPACKE_dggsvd3) ( int matrix_layout, char jobu, 
			char jobv, char jobq, lapack_int m, lapack_int n, lapack_int p, 
			lapack_int * k, lapack_int * l, double * a, lapack_int lda, 
			double * b, lapack_int ldb, double * alpha, double * beta, double * u, 
			lapack_int ldu, double * v, lapack_int ldv, double * q, 
			lapack_int ldq, lapack_int * iwork);
			
    Fptr_NL_LAPACKE_dggsvd3 DGGSVD3;

    dggsvd3_obj = new  ggsvd3_double_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu,
                                         svd_paramslist[idx].jobv,
                                         svd_paramslist[idx].jobq,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].p,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    dggsvd3_obj->threshold = svd_paramslist[idx].svd_threshold;
    dggsvd3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dggsvd3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dggsvd3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dggsvd3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGSVD3 = (Fptr_NL_LAPACKE_dggsvd3)dlsym(dggsvd3_obj->hModule, "LAPACKE_dggsvd3");
    ASSERT_TRUE(DGGSVD3 != NULL) << "failed to ppt the Netlib LAPACKE_dggsvd3 symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dggsvd3_obj->inforef = DGGSVD3( dggsvd3_obj->matrix_layout,
                                    dggsvd3_obj->jobu,
                                    dggsvd3_obj->jobv,
                                    dggsvd3_obj->jobq,
                                    dggsvd3_obj->m,
                                    dggsvd3_obj->n,
                                    dggsvd3_obj->p,
                                    &dggsvd3_obj->kref,
                                    &dggsvd3_obj->lref,
                                    dggsvd3_obj->aref,
                                    dggsvd3_obj->lda,
                                    dggsvd3_obj->bref,
                                    dggsvd3_obj->ldb,
                                    dggsvd3_obj->alpha,
                                    dggsvd3_obj->beta,
                                    dggsvd3_obj->uref,
                                    dggsvd3_obj->ldu,
                                    dggsvd3_obj->vref,
                                    dggsvd3_obj->ldv,
                                    dggsvd3_obj->qref,
                                    dggsvd3_obj->ldq,
                                    dggsvd3_obj->iworkref
                                    );

    /* Compute libflame's Lapacke o/p  */
    dggsvd3_obj->info = LAPACKE_dggsvd3(  dggsvd3_obj->matrix_layout,
                                    dggsvd3_obj->jobu,
                                    dggsvd3_obj->jobv,
                                    dggsvd3_obj->jobq,
                                    dggsvd3_obj->m,
                                    dggsvd3_obj->n,
                                    dggsvd3_obj->p,
                                    &dggsvd3_obj->k,
                                    &dggsvd3_obj->l,
                                    dggsvd3_obj->a,
                                    dggsvd3_obj->lda,
                                    dggsvd3_obj->b,
                                    dggsvd3_obj->ldb,
                                    dggsvd3_obj->alpha,
                                    dggsvd3_obj->beta,
                                    dggsvd3_obj->u,
                                    dggsvd3_obj->ldu,
                                    dggsvd3_obj->v,
                                    dggsvd3_obj->ldv,
                                    dggsvd3_obj->q,
                                    dggsvd3_obj->ldq,
                                    dggsvd3_obj->iwork
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    dggsvd3_obj->diff_a =  computeDiff_d( (dggsvd3_obj->m)*(dggsvd3_obj->n), 
                dggsvd3_obj->a, dggsvd3_obj->aref );

    dggsvd3_obj->diff_b =  computeDiff_d( (dggsvd3_obj->p)*(dggsvd3_obj->n), 
                dggsvd3_obj->b, dggsvd3_obj->bref );

    dggsvd3_obj->diff_u =  computeDiff_d( (dggsvd3_obj->ldu)*(dggsvd3_obj->m), 
                dggsvd3_obj->u, dggsvd3_obj->uref );

    dggsvd3_obj->diff_v =  computeDiff_d( (dggsvd3_obj->p)*(dggsvd3_obj->ldv), 
                dggsvd3_obj->v, dggsvd3_obj->vref );

    dggsvd3_obj->diff_q =  computeDiff_d( (dggsvd3_obj->ldq)*(dggsvd3_obj->n), 
                dggsvd3_obj->q, dggsvd3_obj->qref );

    dggsvd3_obj->diff_alpha =  computeDiff_d( dggsvd3_obj->n, 
                dggsvd3_obj->alpha, dggsvd3_obj->alpharef );

    dggsvd3_obj->diff_beta =  computeDiff_d( dggsvd3_obj->n, 
                dggsvd3_obj->beta, dggsvd3_obj->betaref );
				
    dggsvd3_obj->diff_iwork =  computeDiff_i( dggsvd3_obj->n, 
                dggsvd3_obj->iwork, dggsvd3_obj->iworkref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggsvd3 double: \n diff_a: %f \n diff_b: %f \n \
diff_u: %f \n diff_v: %f \n k: %d \n kref: %d \
diff_q: %f \n l: %d \n lref: %d \n  diff_alpha: %f \n \
diff_beta: %f \n diff_iwork: %d  \n ",
       dggsvd3_obj->diff_a, dggsvd3_obj->diff_b, dggsvd3_obj->diff_u,
       dggsvd3_obj->diff_v, dggsvd3_obj->k, dggsvd3_obj->kref, dggsvd3_obj->diff_q,
       dggsvd3_obj->l, dggsvd3_obj->lref, dggsvd3_obj->diff_alpha,
	   dggsvd3_obj->diff_beta,  dggsvd3_obj->diff_iwork);
#endif
}

TEST_F(dggsvd3_test, dggsvd3_1) {
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_a, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_b, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_q, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_u, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_v, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_alpha, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_beta, dggsvd3_obj->threshold);
    EXPECT_EQ(0, dggsvd3_obj->diff_iwork);
    EXPECT_EQ(dggsvd3_obj->k, dggsvd3_obj->kref);
    EXPECT_EQ(dggsvd3_obj->l, dggsvd3_obj->lref);
}

TEST_F(dggsvd3_test, dggsvd3_2) {
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_a, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_b, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_q, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_u, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_v, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_alpha, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_beta, dggsvd3_obj->threshold);
    EXPECT_EQ(0, dggsvd3_obj->diff_iwork);
    EXPECT_EQ(dggsvd3_obj->k, dggsvd3_obj->kref);
    EXPECT_EQ(dggsvd3_obj->l, dggsvd3_obj->lref);
}

TEST_F(dggsvd3_test, dggsvd3_3) {
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_a, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_b, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_q, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_u, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_v, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_alpha, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_beta, dggsvd3_obj->threshold);
    //EXPECT_EQ(0, dggsvd3_obj->diff_iwork);
    EXPECT_EQ(dggsvd3_obj->k, dggsvd3_obj->kref);
    EXPECT_EQ(dggsvd3_obj->l, dggsvd3_obj->lref);
}

TEST_F(dggsvd3_test, dggsvd3_4) {
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_a, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_b, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_q, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_u, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_v, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_alpha, dggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, dggsvd3_obj->diff_beta, dggsvd3_obj->threshold);
    //EXPECT_EQ(0, dggsvd3_obj->diff_iwork);
    EXPECT_EQ(dggsvd3_obj->k, dggsvd3_obj->kref);
    EXPECT_EQ(dggsvd3_obj->l, dggsvd3_obj->lref);
}

/* Begin ggsvd3_scomplex_common_parameters  class definition */
class ggsvd3_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_u, diff_v, diff_q;
	float diff_alpha, diff_beta;
	int diff_iwork;
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

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p) 
    lapack_int ldq; // The leading dimension of the output array q . ldq≥ max(1, n)
    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains m-by-n matrix A.
    lapack_complex_float* b, *bref; // contains  p-by-n matrix B.

    /* Output parameters */
    lapack_int k, kref; // the dimension of subblocks.
    lapack_int l, lref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    lapack_complex_float* u, *uref;
    lapack_complex_float* v, *vref;
    lapack_complex_float* q, *qref;
	
	/*  generalized singular value pairs of a and b; */
    float* alpha, *alpharef;
    float* beta, *betaref;
	
	/*  stores the sorting information.  */
	lapack_int *iwork, *iworkref;
    /*Return Values */
    int info, inforef;

   public:
      ggsvd3_scomplex_parameters (int matrix_layout_i, char jobu, char jobv,
                                  char jobq, lapack_int m, lapack_int p, 
                                  lapack_int n);
      ~ggsvd3_scomplex_parameters ();
};

/* Constructor definition  ggsvd3 lapack_complex_float_common_parameters */
ggsvd3_scomplex_parameters:: ggsvd3_scomplex_parameters (int matrix_layout_i,
                    char jobu_i, char jobv_i, char jobq_i, lapack_int m_i,  
            lapack_int p_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    jobq = jobq_i;
    n = n_i;
    m = m_i;
    p = p_i;

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
	diff_alpha =0;
	diff_beta = 0;
	diff_iwork = 0;
    k = kref = 0;
    l = lref = 0;


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggsvd3 lapack_complex_float: matrix_layout: %d   jobu: %c \t \
jobv: %c \t jobq: %c  \n m: %d \t p: %d \t n: %d \n",
matrix_layout, jobu_i, jobv_i, jobq_i, m_i, p_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_int_buffer_pair( &iwork, &iworkref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, p*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &u, &uref, ldu*m );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &v, &vref, ldv*p );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_float_buffer_pair( &alpha, &alpharef, n );
    lapacke_gtest_alloc_float_buffer_pair( &beta, &betaref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (v==NULL) || (vref==NULL) || \
        (iwork==NULL) || (iworkref==NULL) || \
        (alpha==NULL) || (alpharef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (q==NULL) || (qref==NULL)  ){
       EXPECT_FALSE( true) << "ggsvd3_scomplex_parameters object: malloc error. Exiting ";
       ggsvd3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, p*n );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( u, uref, ldu*m, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( v, vref, ldv*p, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( q, qref, ldq*n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( alpha, alpharef, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant( iwork, iworkref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'ggsvd3_scomplex_common_parameters' */
ggsvd3_scomplex_parameters :: ~ggsvd3_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggsvd3_free();
} 

//  Test fixture class definition
class cggsvd3_test  : public  ::testing::Test {
public:
   ggsvd3_scomplex_parameters  *cggsvd3_obj;
   void SetUp();  
   void TearDown () { delete cggsvd3_obj; }
};

void cggsvd3_test::SetUp()
{
    /* LAPACKE_cggsvd3 prototype */
    typedef int (*Fptr_NL_LAPACKE_cggsvd3) ( int matrix_layout, char jobu, 
			char jobv, char jobq, lapack_int m, lapack_int n, lapack_int p, 
			lapack_int * k, lapack_int * l, lapack_complex_float * a, lapack_int lda, 
			lapack_complex_float * b, lapack_int ldb, float * alpha, float * beta, lapack_complex_float * u, 
			lapack_int ldu, lapack_complex_float * v, lapack_int ldv, lapack_complex_float * q, 
			lapack_int ldq, lapack_int * iwork);
			
    Fptr_NL_LAPACKE_cggsvd3 CGGSVD3;

    cggsvd3_obj = new  ggsvd3_scomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu,
                                         svd_paramslist[idx].jobv,
                                         svd_paramslist[idx].jobq,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].p,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    cggsvd3_obj->threshold = svd_paramslist[idx].svd_threshold;
    cggsvd3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cggsvd3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cggsvd3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cggsvd3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGSVD3 = (Fptr_NL_LAPACKE_cggsvd3)dlsym(cggsvd3_obj->hModule, "LAPACKE_cggsvd3");
    ASSERT_TRUE(CGGSVD3 != NULL) << "failed to ppt the Netlib LAPACKE_cggsvd3 symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cggsvd3_obj->inforef = CGGSVD3( cggsvd3_obj->matrix_layout,
                                    cggsvd3_obj->jobu,
                                    cggsvd3_obj->jobv,
                                    cggsvd3_obj->jobq,
                                    cggsvd3_obj->m,
                                    cggsvd3_obj->n,
                                    cggsvd3_obj->p,
                                    &cggsvd3_obj->kref,
                                    &cggsvd3_obj->lref,
                                    cggsvd3_obj->aref,
                                    cggsvd3_obj->lda,
                                    cggsvd3_obj->bref,
                                    cggsvd3_obj->ldb,
                                    cggsvd3_obj->alpha,
                                    cggsvd3_obj->beta,
                                    cggsvd3_obj->uref,
                                    cggsvd3_obj->ldu,
                                    cggsvd3_obj->vref,
                                    cggsvd3_obj->ldv,
                                    cggsvd3_obj->qref,
                                    cggsvd3_obj->ldq,
                                    cggsvd3_obj->iworkref
                                    );

    /* Compute libflame's Lapacke o/p  */
    cggsvd3_obj->info = LAPACKE_cggsvd3(  cggsvd3_obj->matrix_layout,
                                    cggsvd3_obj->jobu,
                                    cggsvd3_obj->jobv,
                                    cggsvd3_obj->jobq,
                                    cggsvd3_obj->m,
                                    cggsvd3_obj->n,
                                    cggsvd3_obj->p,
                                    &cggsvd3_obj->k,
                                    &cggsvd3_obj->l,
                                    cggsvd3_obj->a,
                                    cggsvd3_obj->lda,
                                    cggsvd3_obj->b,
                                    cggsvd3_obj->ldb,
                                    cggsvd3_obj->alpha,
                                    cggsvd3_obj->beta,
                                    cggsvd3_obj->u,
                                    cggsvd3_obj->ldu,
                                    cggsvd3_obj->v,
                                    cggsvd3_obj->ldv,
                                    cggsvd3_obj->q,
                                    cggsvd3_obj->ldq,
                                    cggsvd3_obj->iwork
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    cggsvd3_obj->diff_a =  computeDiff_c( (cggsvd3_obj->m)*(cggsvd3_obj->n), 
                cggsvd3_obj->a, cggsvd3_obj->aref );

    cggsvd3_obj->diff_b =  computeDiff_c( (cggsvd3_obj->p)*(cggsvd3_obj->n), 
                cggsvd3_obj->b, cggsvd3_obj->bref );

    cggsvd3_obj->diff_u =  computeDiff_c( (cggsvd3_obj->ldu)*(cggsvd3_obj->m), 
                cggsvd3_obj->u, cggsvd3_obj->uref );

    cggsvd3_obj->diff_v =  computeDiff_c( (cggsvd3_obj->p)*(cggsvd3_obj->ldv), 
                cggsvd3_obj->v, cggsvd3_obj->vref );

    cggsvd3_obj->diff_q =  computeDiff_c( (cggsvd3_obj->ldq)*(cggsvd3_obj->n), 
                cggsvd3_obj->q, cggsvd3_obj->qref );

    cggsvd3_obj->diff_alpha =  computeDiff_s( cggsvd3_obj->n, 
                cggsvd3_obj->alpha, cggsvd3_obj->alpharef );

    cggsvd3_obj->diff_beta =  computeDiff_s( cggsvd3_obj->n, 
                cggsvd3_obj->beta, cggsvd3_obj->betaref );
				
    cggsvd3_obj->diff_iwork =  computeDiff_i( cggsvd3_obj->n, 
                cggsvd3_obj->iwork, cggsvd3_obj->iworkref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggsvd3 lapack_complex_float: \n diff_a: %f \n diff_b: %f \n \
diff_u: %f \n diff_v: %f \n k: %d \n kref: %d \
diff_q: %f \n l: %d \n lref: %d \n  diff_alpha: %f \n \
diff_beta: %f \n diff_iwork: %d  \n ",
       cggsvd3_obj->diff_a, cggsvd3_obj->diff_b, cggsvd3_obj->diff_u,
       cggsvd3_obj->diff_v, cggsvd3_obj->k, cggsvd3_obj->kref, cggsvd3_obj->diff_q,
       cggsvd3_obj->l, cggsvd3_obj->lref, cggsvd3_obj->diff_alpha,
	   cggsvd3_obj->diff_beta,  cggsvd3_obj->diff_iwork);
#endif
}

TEST_F(cggsvd3_test, cggsvd3_1) {
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_a, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_b, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_q, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_u, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_v, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_alpha, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_beta, cggsvd3_obj->threshold);
    EXPECT_EQ(0, cggsvd3_obj->diff_iwork);
    EXPECT_EQ(cggsvd3_obj->k, cggsvd3_obj->kref);
    EXPECT_EQ(cggsvd3_obj->l, cggsvd3_obj->lref);
}

TEST_F(cggsvd3_test, cggsvd3_2) {
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_a, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_b, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_q, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_u, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_v, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_alpha, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_beta, cggsvd3_obj->threshold);
    EXPECT_EQ(0, cggsvd3_obj->diff_iwork);
    EXPECT_EQ(cggsvd3_obj->k, cggsvd3_obj->kref);
    EXPECT_EQ(cggsvd3_obj->l, cggsvd3_obj->lref);
}

TEST_F(cggsvd3_test, cggsvd3_3) {
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_a, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_b, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_q, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_u, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_v, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_alpha, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_beta, cggsvd3_obj->threshold);
    EXPECT_EQ(0, cggsvd3_obj->diff_iwork);
    EXPECT_EQ(cggsvd3_obj->k, cggsvd3_obj->kref);
    EXPECT_EQ(cggsvd3_obj->l, cggsvd3_obj->lref);
}

TEST_F(cggsvd3_test, cggsvd3_4) {
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_a, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_b, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_q, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_u, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_v, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_alpha, cggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, cggsvd3_obj->diff_beta, cggsvd3_obj->threshold);
    EXPECT_EQ(0, cggsvd3_obj->diff_iwork);
    EXPECT_EQ(cggsvd3_obj->k, cggsvd3_obj->kref);
    EXPECT_EQ(cggsvd3_obj->l, cggsvd3_obj->lref);
}

/* Begin ggsvd3_dcomplex_common_parameters  class definition */
class ggsvd3_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b, diff_u, diff_v, diff_q;
	double diff_alpha, diff_beta;
	int diff_iwork;
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

    lapack_int ldu; // leading dimension of the output array u . ldu≥ max(1, m)
    lapack_int ldv; // The leading dimension of the output array v . ldv≥ max(1, p) 
    lapack_int ldq; // The leading dimension of the output array q . ldq≥ max(1, n)
    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains m-by-n matrix A.
    lapack_complex_double* b, *bref; // contains  p-by-n matrix B.

    /* Output parameters */
    lapack_int k, kref; // the dimension of subblocks.
    lapack_int l, lref; // the dimension of subblocks.
    // Below buffers contain the o/p orthogonal/unitary matrices
    lapack_complex_double* u, *uref;
    lapack_complex_double* v, *vref;
    lapack_complex_double* q, *qref;
	
	/*  generalized singular value pairs of a and b; */
    double* alpha, *alpharef;
    double* beta, *betaref;
	
	/*  stores the sorting information.  */
	lapack_int *iwork, *iworkref;
    /*Return Values */
    int info, inforef;

   public:
      ggsvd3_dcomplex_parameters (int matrix_layout_i, char jobu, char jobv,
                                  char jobq, lapack_int m, lapack_int p, 
                                  lapack_int n);
      ~ggsvd3_dcomplex_parameters ();
};

/* Constructor definition  ggsvd3 lapack_complex_double_common_parameters */
ggsvd3_dcomplex_parameters:: ggsvd3_dcomplex_parameters (int matrix_layout_i,
                    char jobu_i, char jobv_i, char jobq_i, lapack_int m_i,  
            lapack_int p_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    jobu = jobu_i;
    jobv = jobv_i;
    jobq = jobq_i;
    n = n_i;
    m = m_i;
    p = p_i;

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
	diff_alpha =0;
	diff_beta = 0;
	diff_iwork = 0;
    k = kref = 0;
    l = lref = 0;


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggsvd3 lapack_complex_double: matrix_layout: %d   jobu: %c \t \
jobv: %c \t jobq: %c  \n m: %d \t p: %d \t n: %d \n",
matrix_layout, jobu_i, jobv_i, jobq_i, m_i, p_i, n_i);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_int_buffer_pair( &iwork, &iworkref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, m*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, p*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &u, &uref, ldu*m );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &v, &vref, ldv*p );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_double_buffer_pair( &alpha, &alpharef, n );
    lapacke_gtest_alloc_double_buffer_pair( &beta, &betaref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (u==NULL) || (uref==NULL) || \
        (v==NULL) || (vref==NULL) || \
        (iwork==NULL) || (iworkref==NULL) || \
        (alpha==NULL) || (alpharef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (q==NULL) || (qref==NULL)  ){
       EXPECT_FALSE( true) << "ggsvd3_dcomplex_parameters object: malloc error. Exiting ";
       ggsvd3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, m*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, p*n );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( u, uref, ldu*m, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( v, vref, ldv*p, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( q, qref, ldq*n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( alpha, alpharef, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant( iwork, iworkref, n, 0);

   } /* end of Constructor  */

/* Destructor definition  'ggsvd3_dcomplex_common_parameters' */
ggsvd3_dcomplex_parameters :: ~ggsvd3_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggsvd3_free();
} 

//  Test fixture class definition
class zggsvd3_test  : public  ::testing::Test {
public:
   ggsvd3_dcomplex_parameters  *zggsvd3_obj;
   void SetUp();  
   void TearDown () { delete zggsvd3_obj; }
};

void zggsvd3_test::SetUp()
{
    /* LAPACKE_zggsvd3 prototype */
    typedef int (*Fptr_NL_LAPACKE_zggsvd3) ( int matrix_layout, char jobu, 
			char jobv, char jobq, lapack_int m, lapack_int n, lapack_int p, 
			lapack_int * k, lapack_int * l, lapack_complex_double * a, lapack_int lda, 
			lapack_complex_double * b, lapack_int ldb, double * alpha, double * beta, lapack_complex_double * u, 
			lapack_int ldu, lapack_complex_double * v, lapack_int ldv, lapack_complex_double * q, 
			lapack_int ldq, lapack_int * iwork);
			
    Fptr_NL_LAPACKE_zggsvd3 ZGGSVD3;

    zggsvd3_obj = new  ggsvd3_dcomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu,
                                         svd_paramslist[idx].jobv,
                                         svd_paramslist[idx].jobq,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].p,
                                         svd_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);
    zggsvd3_obj->threshold = svd_paramslist[idx].svd_threshold;
    zggsvd3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zggsvd3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zggsvd3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zggsvd3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGSVD3 = (Fptr_NL_LAPACKE_zggsvd3)dlsym(zggsvd3_obj->hModule, "LAPACKE_zggsvd3");
    ASSERT_TRUE(ZGGSVD3 != NULL) << "failed to ppt the Netlib LAPACKE_zggsvd3 symbol";

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zggsvd3_obj->inforef = ZGGSVD3( zggsvd3_obj->matrix_layout,
                                    zggsvd3_obj->jobu,
                                    zggsvd3_obj->jobv,
                                    zggsvd3_obj->jobq,
                                    zggsvd3_obj->m,
                                    zggsvd3_obj->n,
                                    zggsvd3_obj->p,
                                    &zggsvd3_obj->kref,
                                    &zggsvd3_obj->lref,
                                    zggsvd3_obj->aref,
                                    zggsvd3_obj->lda,
                                    zggsvd3_obj->bref,
                                    zggsvd3_obj->ldb,
                                    zggsvd3_obj->alpha,
                                    zggsvd3_obj->beta,
                                    zggsvd3_obj->uref,
                                    zggsvd3_obj->ldu,
                                    zggsvd3_obj->vref,
                                    zggsvd3_obj->ldv,
                                    zggsvd3_obj->qref,
                                    zggsvd3_obj->ldq,
                                    zggsvd3_obj->iworkref
                                    );

    /* Compute libflame's Lapacke o/p  */
    zggsvd3_obj->info = LAPACKE_zggsvd3(  zggsvd3_obj->matrix_layout,
                                    zggsvd3_obj->jobu,
                                    zggsvd3_obj->jobv,
                                    zggsvd3_obj->jobq,
                                    zggsvd3_obj->m,
                                    zggsvd3_obj->n,
                                    zggsvd3_obj->p,
                                    &zggsvd3_obj->k,
                                    &zggsvd3_obj->l,
                                    zggsvd3_obj->a,
                                    zggsvd3_obj->lda,
                                    zggsvd3_obj->b,
                                    zggsvd3_obj->ldb,
                                    zggsvd3_obj->alpha,
                                    zggsvd3_obj->beta,
                                    zggsvd3_obj->u,
                                    zggsvd3_obj->ldu,
                                    zggsvd3_obj->v,
                                    zggsvd3_obj->ldv,
                                    zggsvd3_obj->q,
                                    zggsvd3_obj->ldq,
                                    zggsvd3_obj->iwork
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    zggsvd3_obj->diff_a =  computeDiff_z( (zggsvd3_obj->m)*(zggsvd3_obj->n), 
                zggsvd3_obj->a, zggsvd3_obj->aref );

    zggsvd3_obj->diff_b =  computeDiff_z( (zggsvd3_obj->p)*(zggsvd3_obj->n), 
                zggsvd3_obj->b, zggsvd3_obj->bref );

    zggsvd3_obj->diff_u =  computeDiff_z( (zggsvd3_obj->ldu)*(zggsvd3_obj->m), 
                zggsvd3_obj->u, zggsvd3_obj->uref );

    zggsvd3_obj->diff_v =  computeDiff_z( (zggsvd3_obj->p)*(zggsvd3_obj->ldv), 
                zggsvd3_obj->v, zggsvd3_obj->vref );

    zggsvd3_obj->diff_q =  computeDiff_z( (zggsvd3_obj->ldq)*(zggsvd3_obj->n), 
                zggsvd3_obj->q, zggsvd3_obj->qref );

    zggsvd3_obj->diff_alpha =  computeDiff_d( zggsvd3_obj->n, 
                zggsvd3_obj->alpha, zggsvd3_obj->alpharef );

    zggsvd3_obj->diff_beta =  computeDiff_d( zggsvd3_obj->n, 
                zggsvd3_obj->beta, zggsvd3_obj->betaref );
				
    zggsvd3_obj->diff_iwork =  computeDiff_i( zggsvd3_obj->n, 
                zggsvd3_obj->iwork, zggsvd3_obj->iworkref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggsvd3 lapack_complex_double: \n diff_a: %f \n diff_b: %f \n \
diff_u: %f \n diff_v: %f \n k: %d \n kref: %d \
diff_q: %f \n l: %d \n lref: %d \n  diff_alpha: %f \n \
diff_beta: %f \n diff_iwork: %d  \n ",
       zggsvd3_obj->diff_a, zggsvd3_obj->diff_b, zggsvd3_obj->diff_u,
       zggsvd3_obj->diff_v, zggsvd3_obj->k, zggsvd3_obj->kref, zggsvd3_obj->diff_q,
       zggsvd3_obj->l, zggsvd3_obj->lref, zggsvd3_obj->diff_alpha,
	   zggsvd3_obj->diff_beta,  zggsvd3_obj->diff_iwork);
#endif
}

TEST_F(zggsvd3_test, zggsvd3_1) {
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_a, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_b, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_q, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_u, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_v, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_alpha, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_beta, zggsvd3_obj->threshold);
    EXPECT_EQ(0, zggsvd3_obj->diff_iwork);
    EXPECT_EQ(zggsvd3_obj->k, zggsvd3_obj->kref);
    EXPECT_EQ(zggsvd3_obj->l, zggsvd3_obj->lref);
}

TEST_F(zggsvd3_test, zggsvd3_2) {
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_a, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_b, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_q, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_u, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_v, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_alpha, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_beta, zggsvd3_obj->threshold);
    EXPECT_EQ(0, zggsvd3_obj->diff_iwork);
    EXPECT_EQ(zggsvd3_obj->k, zggsvd3_obj->kref);
    EXPECT_EQ(zggsvd3_obj->l, zggsvd3_obj->lref);
}

TEST_F(zggsvd3_test, zggsvd3_3) {
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_a, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_b, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_q, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_u, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_v, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_alpha, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_beta, zggsvd3_obj->threshold);
    EXPECT_EQ(0, zggsvd3_obj->diff_iwork);
    EXPECT_EQ(zggsvd3_obj->k, zggsvd3_obj->kref);
    EXPECT_EQ(zggsvd3_obj->l, zggsvd3_obj->lref);
}

TEST_F(zggsvd3_test, zggsvd3_4) {
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_a, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_b, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_q, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_u, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_v, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_alpha, zggsvd3_obj->threshold);
    EXPECT_NEAR(0.0, zggsvd3_obj->diff_beta, zggsvd3_obj->threshold);
    EXPECT_EQ(0, zggsvd3_obj->diff_iwork);
    EXPECT_EQ(zggsvd3_obj->k, zggsvd3_obj->kref);
    EXPECT_EQ(zggsvd3_obj->l, zggsvd3_obj->lref);
}


