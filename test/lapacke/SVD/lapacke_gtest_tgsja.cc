#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define tgsja_free() \
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

/* Begin tgsja_float_common_parameters  class definition */
class tgsja_float_parameters{

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
	
	/*  generalized singular value pairs of a and b; */
    float* alpha, *alpharef;
    float* beta, *betaref;
	
	/*  stores the sorting information.  */
	lapack_int ncycle, ncycleref, *iwork, *iworkref;
	
    /*Return Values */
    int info, inforef;

   public:
      tgsja_float_parameters (int matrix_layout_i, char jobu, char jobv,
                                  char jobq, lapack_int m, lapack_int p, 
                                  lapack_int n, float tola, float tolb);
      ~tgsja_float_parameters ();
};

/* Constructor definition  tgsja float_common_parameters */
tgsja_float_parameters:: tgsja_float_parameters (int matrix_layout_i,
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
	diff_alpha =0;
	diff_beta = 0;
	diff_iwork = 0;
    k = kref = 0;
    l = lref = 0;
	ncycle = 0;
	ncycleref = 0;


#if LAPACKE_TEST_VERBOSE
   printf(" \n tgsja float: matrix_layout: %d   jobu: %c \t \
jobv: %c \t jobq: %c  \n m: %d \t p: %d \t n: %d \n tola: %f  \t tolb: %f \n",
matrix_layout, jobu_i, jobv_i, jobq_i, m_i,
p_i, n_i, tola_i, tolb_i);
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
       EXPECT_FALSE( true) << "tgsja_float_parameters object: malloc error. Exiting ";
       tgsja_free();
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

/* Destructor definition  'tgsja_float_common_parameters' */
tgsja_float_parameters :: ~tgsja_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgsja_free();
} 

//  Test fixture class definition
class stgsja_test  : public  ::testing::Test {
public:
   tgsja_float_parameters  *stgsja_obj;
   void SetUp();  
   void TearDown () { delete stgsja_obj; }
};

void stgsja_test::SetUp()
{
    /* LAPACKE stgsja prototype */
    typedef int (*Fptr_NL_LAPACKE_stgsja) (  int matrix_layout, 
		char jobu, char jobv, char jobq, lapack_int m, lapack_int p,
		lapack_int n, lapack_int k, lapack_int l, float* a,
		lapack_int lda, float* b, lapack_int ldb, float tola, float tolb,
		float* alpha, float* beta, float* u, lapack_int ldu, float* v,
		lapack_int ldv, float* q, lapack_int ldq, lapack_int* ncycle);

    Fptr_NL_LAPACKE_stgsja STGSJA;

    /* LAPACKE stgsja prototype */
    typedef int (*Fptr_NL_LAPACKE_sggsvp) ( int matrix_layout, char jobu, 
            char jobv, char jobq, lapack_int m, lapack_int p, lapack_int n, 
            float* a, lapack_int lda, float* b, lapack_int ldb, float tola,
            float tolb, lapack_int* k, lapack_int* l, float* u, lapack_int ldu,
            float* v, lapack_int ldv, float* q, lapack_int ldqe);

    Fptr_NL_LAPACKE_sggsvp SGGSVP;	

    stgsja_obj = new  tgsja_float_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu,
                                         svd_paramslist[idx].jobv,
                                         svd_paramslist[idx].jobq,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].p,
                                         svd_paramslist[idx].n,
                                         svd_paramslist[idx].tola,
                                         svd_paramslist[idx].tolb );

    idx = Circular_Increment_Index(idx);
    stgsja_obj->threshold = svd_paramslist[idx].svd_threshold;
    stgsja_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stgsja_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stgsja_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stgsja_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    STGSJA = (Fptr_NL_LAPACKE_stgsja)dlsym(stgsja_obj->hModule, "LAPACKE_stgsja");
    ASSERT_TRUE(STGSJA != NULL) << "failed to ppt the Netlib LAPACKE_stgsja symbol";

    SGGSVP = (Fptr_NL_LAPACKE_sggsvp)dlsym(stgsja_obj->hModule, "LAPACKE_sggsvp3");
    ASSERT_TRUE(SGGSVP != NULL) << "failed to ppt the Netlib LAPACKE_sggsvp symbol";

    /* Preprocessing for generation of i/ps for reference and libflame apis
       by invoking the below API for both i/p sets 	*/
    stgsja_obj->inforef = SGGSVP(   stgsja_obj->matrix_layout,
                                    stgsja_obj->jobu,
                                    stgsja_obj->jobv,
                                    stgsja_obj->jobq,
                                    stgsja_obj->m,
                                    stgsja_obj->p,
                                    stgsja_obj->n,
                                    stgsja_obj->aref,
                                    stgsja_obj->lda,
                                    stgsja_obj->bref,
                                    stgsja_obj->ldb,
                                    stgsja_obj->tola,
                                    stgsja_obj->tolb,
                                    &stgsja_obj->kref,
                                    &stgsja_obj->lref,
                                    stgsja_obj->uref,
                                    stgsja_obj->ldu,
                                    stgsja_obj->vref,
                                    stgsja_obj->ldv,
                                    stgsja_obj->qref,
                                    stgsja_obj->ldq
                                    );

    stgsja_obj->info = SGGSVP(   stgsja_obj->matrix_layout,
                                    stgsja_obj->jobu,
                                    stgsja_obj->jobv,
                                    stgsja_obj->jobq,
                                    stgsja_obj->m,
                                    stgsja_obj->p,
                                    stgsja_obj->n,
                                    stgsja_obj->a,
                                    stgsja_obj->lda,
                                    stgsja_obj->b,
                                    stgsja_obj->ldb,
                                    stgsja_obj->tola,
                                    stgsja_obj->tolb,
                                    &stgsja_obj->k,
                                    &stgsja_obj->l,
                                    stgsja_obj->u,
                                    stgsja_obj->ldu,
                                    stgsja_obj->v,
                                    stgsja_obj->ldv,
                                    stgsja_obj->q,
                                    stgsja_obj->ldq
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    stgsja_obj->inforef = STGSJA(   stgsja_obj->matrix_layout,
                                    stgsja_obj->jobu,
                                    stgsja_obj->jobv,
                                    stgsja_obj->jobq,
                                    stgsja_obj->m,
                                    stgsja_obj->p,
                                    stgsja_obj->n,
                                    stgsja_obj->kref,
                                    stgsja_obj->lref,
                                    stgsja_obj->aref,
                                    stgsja_obj->lda,
                                    stgsja_obj->bref,
                                    stgsja_obj->ldb,
                                    stgsja_obj->tola,
                                    stgsja_obj->tolb,
                                    stgsja_obj->alpha,
                                    stgsja_obj->beta,
                                    stgsja_obj->uref,
                                    stgsja_obj->ldu,
                                    stgsja_obj->vref,
                                    stgsja_obj->ldv,
                                    stgsja_obj->qref,
                                    stgsja_obj->ldq,
                                    &stgsja_obj->ncycleref
                                    );

    /* Compute libflame's Lapacke o/p  */
    stgsja_obj->info = LAPACKE_stgsja(  stgsja_obj->matrix_layout,
                                    stgsja_obj->jobu,
                                    stgsja_obj->jobv,
                                    stgsja_obj->jobq,
                                    stgsja_obj->m,
                                    stgsja_obj->p,
                                    stgsja_obj->n,
                                    stgsja_obj->k,
                                    stgsja_obj->l,
                                    stgsja_obj->a,
                                    stgsja_obj->lda,
                                    stgsja_obj->b,
                                    stgsja_obj->ldb,
                                    stgsja_obj->tola,
                                    stgsja_obj->tolb,
                                    stgsja_obj->alpha,
                                    stgsja_obj->beta,
                                    stgsja_obj->u,
                                    stgsja_obj->ldu,
                                    stgsja_obj->v,
                                    stgsja_obj->ldv,
                                    stgsja_obj->q,
                                    stgsja_obj->ldq,
                                    &stgsja_obj->ncycle
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    stgsja_obj->diff_a =  computeDiff_s( (stgsja_obj->m)*(stgsja_obj->n), 
                stgsja_obj->a, stgsja_obj->aref );

    stgsja_obj->diff_b =  computeDiff_s( (stgsja_obj->p)*(stgsja_obj->n), 
                stgsja_obj->b, stgsja_obj->bref );

    stgsja_obj->diff_u =  computeDiff_s( (stgsja_obj->ldu)*(stgsja_obj->m), 
                stgsja_obj->u, stgsja_obj->uref );

    stgsja_obj->diff_v =  computeDiff_s( (stgsja_obj->p)*(stgsja_obj->ldv), 
                stgsja_obj->v, stgsja_obj->vref );

    stgsja_obj->diff_q =  computeDiff_s( (stgsja_obj->ldq)*(stgsja_obj->n), 
                stgsja_obj->q, stgsja_obj->qref );

    stgsja_obj->diff_alpha =  computeDiff_s( stgsja_obj->n, 
                stgsja_obj->alpha, stgsja_obj->alpharef );

    stgsja_obj->diff_beta =  computeDiff_s( stgsja_obj->n, 
                stgsja_obj->beta, stgsja_obj->betaref );
				

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgsja float: \n diff_a: %f \n diff_b: %f \n \
diff_u: %f \n diff_v: %f \n k: %d \n kref: %d \
diff_q: %f \n l: %d \n lref: %d \n ",
       stgsja_obj->diff_a, stgsja_obj->diff_b, stgsja_obj->diff_u,
       stgsja_obj->diff_v, stgsja_obj->k, stgsja_obj->kref, stgsja_obj->diff_q,
       stgsja_obj->l, stgsja_obj->lref );
#endif
}

TEST_F(stgsja_test, stgsja1) {
    EXPECT_NEAR(0.0, stgsja_obj->diff_a, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_b, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_q, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_u, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_v, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_alpha, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_beta, stgsja_obj->threshold);
    EXPECT_EQ(stgsja_obj->ncycle, stgsja_obj->ncycleref);
}

TEST_F(stgsja_test, stgsja2) {
    EXPECT_NEAR(0.0, stgsja_obj->diff_a, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_b, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_q, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_u, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_v, stgsja_obj->threshold);
    EXPECT_EQ(stgsja_obj->ncycle, stgsja_obj->ncycleref);
}

TEST_F(stgsja_test, stgsja) {
    EXPECT_NEAR(0.0, stgsja_obj->diff_a, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_b, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_q, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_u, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_v, stgsja_obj->threshold);
    EXPECT_EQ(stgsja_obj->ncycle, stgsja_obj->ncycleref);
}

TEST_F(stgsja_test, stgsja4) {
    EXPECT_NEAR(0.0, stgsja_obj->diff_a, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_b, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_q, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_u, stgsja_obj->threshold);
    EXPECT_NEAR(0.0, stgsja_obj->diff_v, stgsja_obj->threshold);
    EXPECT_EQ(stgsja_obj->ncycle, stgsja_obj->ncycleref);
}


/* Begin tgsja_double_common_parameters  class definition */
class tgsja_double_parameters{

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

    // thresholds to determine the effective numerical rank of matrix B and a subblock of A.
    double tola;
    double tolb;
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
	lapack_int ncycle, ncycleref, *iwork, *iworkref;
	
    /*Return Values */
    int info, inforef;

   public:
      tgsja_double_parameters (int matrix_layout_i, char jobu, char jobv,
                                  char jobq, lapack_int m, lapack_int p, 
                                  lapack_int n, double tola, double tolb);
      ~tgsja_double_parameters ();
};

/* Constructor definition  tgsja double_common_parameters */
tgsja_double_parameters:: tgsja_double_parameters (int matrix_layout_i,
                    char jobu_i, char jobv_i, char jobq_i, lapack_int m_i,  
            lapack_int p_i, lapack_int n_i, double tola_i, double tolb_i)
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
	diff_alpha =0;
	diff_beta = 0;
	diff_iwork = 0;
    k = kref = 0;
    l = lref = 0;
	ncycle = 0;
	ncycleref = 0;


#if LAPACKE_TEST_VERBOSE
   printf(" \n tgsja double: matrix_layout: %d   jobu: %c \t \
jobv: %c \t jobq: %c  \n m: %d \t p: %d \t n: %d \n tola: %f  \t tolb: %f \n",
matrix_layout, jobu_i, jobv_i, jobq_i, m_i,
p_i, n_i, tola_i, tolb_i);
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
       EXPECT_FALSE( true) << "tgsja_double_parameters object: malloc error. Exiting ";
       tgsja_free();
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

/* Destructor definition  'tgsja_double_common_parameters' */
tgsja_double_parameters :: ~tgsja_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgsja_free();
} 

//  Test fixture class definition
class dtgsja_test  : public  ::testing::Test {
public:
   tgsja_double_parameters  *dtgsja_obj;
   void SetUp();  
   void TearDown () { delete dtgsja_obj; }
};

void dtgsja_test::SetUp()
{
    /* LAPACKE dtgsja prototype */
    typedef int (*Fptr_NL_LAPACKE_dtgsja) (  int matrix_layout, 
		char jobu, char jobv, char jobq, lapack_int m, lapack_int p,
		lapack_int n, lapack_int k, lapack_int l, double* a,
		lapack_int lda, double* b, lapack_int ldb, double tola, double tolb,
		double* alpha, double* beta, double* u, lapack_int ldu, double* v,
		lapack_int ldv, double* q, lapack_int ldq, lapack_int* ncycle);

    Fptr_NL_LAPACKE_dtgsja DTGSJA;

    /* LAPACKE dtgsja prototype */
    typedef int (*Fptr_NL_LAPACKE_dggsvp) ( int matrix_layout, char jobu, 
            char jobv, char jobq, lapack_int m, lapack_int p, lapack_int n, 
            double* a, lapack_int lda, double* b, lapack_int ldb, double tola,
            double tolb, lapack_int* k, lapack_int* l, double* u, lapack_int ldu,
            double* v, lapack_int ldv, double* q, lapack_int ldqe);

    Fptr_NL_LAPACKE_dggsvp DGGSVP;	

    dtgsja_obj = new  tgsja_double_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu,
                                         svd_paramslist[idx].jobv,
                                         svd_paramslist[idx].jobq,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].p,
                                         svd_paramslist[idx].n,
                                         svd_paramslist[idx].tola,
                                         svd_paramslist[idx].tolb );

    idx = Circular_Increment_Index(idx);
    dtgsja_obj->threshold = svd_paramslist[idx].svd_threshold;
    dtgsja_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtgsja_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtgsja_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtgsja_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DTGSJA = (Fptr_NL_LAPACKE_dtgsja)dlsym(dtgsja_obj->hModule, "LAPACKE_dtgsja");
    ASSERT_TRUE(DTGSJA != NULL) << "failed to ppt the Netlib LAPACKE_dtgsja symbol";

    DGGSVP = (Fptr_NL_LAPACKE_dggsvp)dlsym(dtgsja_obj->hModule, "LAPACKE_dggsvp3");
    ASSERT_TRUE(DGGSVP != NULL) << "failed to ppt the Netlib LAPACKE_dggsvp symbol";

    /* Preprocessing for generation of i/ps for reference and libflame apis
       by invoking the below API for both i/p sets 	*/
    dtgsja_obj->inforef = DGGSVP(   dtgsja_obj->matrix_layout,
                                    dtgsja_obj->jobu,
                                    dtgsja_obj->jobv,
                                    dtgsja_obj->jobq,
                                    dtgsja_obj->m,
                                    dtgsja_obj->p,
                                    dtgsja_obj->n,
                                    dtgsja_obj->aref,
                                    dtgsja_obj->lda,
                                    dtgsja_obj->bref,
                                    dtgsja_obj->ldb,
                                    dtgsja_obj->tola,
                                    dtgsja_obj->tolb,
                                    &dtgsja_obj->kref,
                                    &dtgsja_obj->lref,
                                    dtgsja_obj->uref,
                                    dtgsja_obj->ldu,
                                    dtgsja_obj->vref,
                                    dtgsja_obj->ldv,
                                    dtgsja_obj->qref,
                                    dtgsja_obj->ldq
                                    );

    dtgsja_obj->info = DGGSVP(   dtgsja_obj->matrix_layout,
                                    dtgsja_obj->jobu,
                                    dtgsja_obj->jobv,
                                    dtgsja_obj->jobq,
                                    dtgsja_obj->m,
                                    dtgsja_obj->p,
                                    dtgsja_obj->n,
                                    dtgsja_obj->a,
                                    dtgsja_obj->lda,
                                    dtgsja_obj->b,
                                    dtgsja_obj->ldb,
                                    dtgsja_obj->tola,
                                    dtgsja_obj->tolb,
                                    &dtgsja_obj->k,
                                    &dtgsja_obj->l,
                                    dtgsja_obj->u,
                                    dtgsja_obj->ldu,
                                    dtgsja_obj->v,
                                    dtgsja_obj->ldv,
                                    dtgsja_obj->q,
                                    dtgsja_obj->ldq
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dtgsja_obj->inforef = DTGSJA(   dtgsja_obj->matrix_layout,
                                    dtgsja_obj->jobu,
                                    dtgsja_obj->jobv,
                                    dtgsja_obj->jobq,
                                    dtgsja_obj->m,
                                    dtgsja_obj->p,
                                    dtgsja_obj->n,
                                    dtgsja_obj->kref,
                                    dtgsja_obj->lref,
                                    dtgsja_obj->aref,
                                    dtgsja_obj->lda,
                                    dtgsja_obj->bref,
                                    dtgsja_obj->ldb,
                                    dtgsja_obj->tola,
                                    dtgsja_obj->tolb,
                                    dtgsja_obj->alpha,
                                    dtgsja_obj->beta,
                                    dtgsja_obj->uref,
                                    dtgsja_obj->ldu,
                                    dtgsja_obj->vref,
                                    dtgsja_obj->ldv,
                                    dtgsja_obj->qref,
                                    dtgsja_obj->ldq,
                                    &dtgsja_obj->ncycleref
                                    );

    /* Compute libflame's Lapacke o/p  */
    dtgsja_obj->info = LAPACKE_dtgsja(  dtgsja_obj->matrix_layout,
                                    dtgsja_obj->jobu,
                                    dtgsja_obj->jobv,
                                    dtgsja_obj->jobq,
                                    dtgsja_obj->m,
                                    dtgsja_obj->p,
                                    dtgsja_obj->n,
                                    dtgsja_obj->k,
                                    dtgsja_obj->l,
                                    dtgsja_obj->a,
                                    dtgsja_obj->lda,
                                    dtgsja_obj->b,
                                    dtgsja_obj->ldb,
                                    dtgsja_obj->tola,
                                    dtgsja_obj->tolb,
                                    dtgsja_obj->alpha,
                                    dtgsja_obj->beta,
                                    dtgsja_obj->u,
                                    dtgsja_obj->ldu,
                                    dtgsja_obj->v,
                                    dtgsja_obj->ldv,
                                    dtgsja_obj->q,
                                    dtgsja_obj->ldq,
                                    &dtgsja_obj->ncycle
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    dtgsja_obj->diff_a =  computeDiff_d( (dtgsja_obj->m)*(dtgsja_obj->n), 
                dtgsja_obj->a, dtgsja_obj->aref );

    dtgsja_obj->diff_b =  computeDiff_d( (dtgsja_obj->p)*(dtgsja_obj->n), 
                dtgsja_obj->b, dtgsja_obj->bref );

    dtgsja_obj->diff_u =  computeDiff_d( (dtgsja_obj->ldu)*(dtgsja_obj->m), 
                dtgsja_obj->u, dtgsja_obj->uref );

    dtgsja_obj->diff_v =  computeDiff_d( (dtgsja_obj->p)*(dtgsja_obj->ldv), 
                dtgsja_obj->v, dtgsja_obj->vref );

    dtgsja_obj->diff_q =  computeDiff_d( (dtgsja_obj->ldq)*(dtgsja_obj->n), 
                dtgsja_obj->q, dtgsja_obj->qref );

    dtgsja_obj->diff_alpha =  computeDiff_d( dtgsja_obj->n, 
                dtgsja_obj->alpha, dtgsja_obj->alpharef );

    dtgsja_obj->diff_beta =  computeDiff_d( dtgsja_obj->n, 
                dtgsja_obj->beta, dtgsja_obj->betaref );
				

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgsja double: \n diff_a: %f \n diff_b: %f \n \
diff_u: %f \n diff_v: %f \n k: %d \n kref: %d \
diff_q: %f \n l: %d \n lref: %d \n ",
       dtgsja_obj->diff_a, dtgsja_obj->diff_b, dtgsja_obj->diff_u,
       dtgsja_obj->diff_v, dtgsja_obj->k, dtgsja_obj->kref, dtgsja_obj->diff_q,
       dtgsja_obj->l, dtgsja_obj->lref );
#endif
}

TEST_F(dtgsja_test, dtgsja1) {
    EXPECT_NEAR(0.0, dtgsja_obj->diff_a, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_b, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_q, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_u, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_v, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_alpha, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_beta, dtgsja_obj->threshold);
    EXPECT_EQ(dtgsja_obj->ncycle, dtgsja_obj->ncycleref);
}

TEST_F(dtgsja_test, dtgsja2) {
    EXPECT_NEAR(0.0, dtgsja_obj->diff_a, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_b, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_q, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_u, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_v, dtgsja_obj->threshold);
    EXPECT_EQ(dtgsja_obj->ncycle, dtgsja_obj->ncycleref);
}

TEST_F(dtgsja_test, dtgsja) {
    EXPECT_NEAR(0.0, dtgsja_obj->diff_a, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_b, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_q, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_u, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_v, dtgsja_obj->threshold);
    EXPECT_EQ(dtgsja_obj->ncycle, dtgsja_obj->ncycleref);
}

TEST_F(dtgsja_test, dtgsja4) {
    EXPECT_NEAR(0.0, dtgsja_obj->diff_a, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_b, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_q, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_u, dtgsja_obj->threshold);
    EXPECT_NEAR(0.0, dtgsja_obj->diff_v, dtgsja_obj->threshold);
    EXPECT_EQ(dtgsja_obj->ncycle, dtgsja_obj->ncycleref);
}

/* Begin tgsja_scomplex_common_parameters  class definition */
class tgsja_scomplex_parameters{

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

    // thresholds to determine the effective numerical rank of matrix B and a subblock of A.
    float tola;
    float tolb;
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
	lapack_int ncycle, ncycleref, *iwork, *iworkref;
	
    /*Return Values */
    int info, inforef;

   public:
      tgsja_scomplex_parameters (int matrix_layout_i, char jobu, char jobv,
                                  char jobq, lapack_int m, lapack_int p, 
                                  lapack_int n, float tola, float tolb);
      ~tgsja_scomplex_parameters ();
};

/* Constructor definition  tgsja lapack_complex_float_common_parameters */
tgsja_scomplex_parameters:: tgsja_scomplex_parameters (int matrix_layout_i,
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
	diff_alpha =0;
	diff_beta = 0;
	diff_iwork = 0;
    k = kref = 0;
    l = lref = 0;
	ncycle = 0;
	ncycleref = 0;


#if LAPACKE_TEST_VERBOSE
   printf(" \n tgsja lapack_complex_float: matrix_layout: %d   jobu: %c \t \
jobv: %c \t jobq: %c  \n m: %d \t p: %d \t n: %d \n tola: %f  \t tolb: %f \n",
matrix_layout, jobu_i, jobv_i, jobq_i, m_i,
p_i, n_i, tola_i, tolb_i);
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
       EXPECT_FALSE( true) << "tgsja_scomplex_parameters object: malloc error. Exiting ";
       tgsja_free();
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

/* Destructor definition  'tgsja_scomplex_common_parameters' */
tgsja_scomplex_parameters :: ~tgsja_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgsja_free();
} 

//  Test fixture class definition
class ctgsja_test  : public  ::testing::Test {
public:
   tgsja_scomplex_parameters  *ctgsja_obj;
   void SetUp();  
   void TearDown () { delete ctgsja_obj; }
};

void ctgsja_test::SetUp()
{
    /* LAPACKE ctgsja prototype */
    typedef int (*Fptr_NL_LAPACKE_ctgsja) (  int matrix_layout, 
		char jobu, char jobv, char jobq, lapack_int m, lapack_int p,
		lapack_int n, lapack_int k, lapack_int l, lapack_complex_float* a,
		lapack_int lda, lapack_complex_float* b, lapack_int ldb, float tola,
		float tolb, float* alpha, float* beta, lapack_complex_float* u,
		lapack_int ldu, lapack_complex_float* v, lapack_int ldv,
		lapack_complex_float* q, lapack_int ldq, lapack_int* ncycle);

    Fptr_NL_LAPACKE_ctgsja CTGSJA;

    /* LAPACKE ctgsja prototype */
    typedef int (*Fptr_NL_LAPACKE_cggsvp) ( int matrix_layout, char jobu,
	char jobv, char jobq, lapack_int m, lapack_int p, lapack_int n,
	lapack_complex_float* a, lapack_int lda, lapack_complex_float* b,
	lapack_int ldb, float tola, float tolb, lapack_int* k, lapack_int* l,
	lapack_complex_float* u, lapack_int ldu, lapack_complex_float* v,
	lapack_int ldv, lapack_complex_float* q, lapack_int ldq);

    Fptr_NL_LAPACKE_cggsvp CGGSVP;	

    ctgsja_obj = new  tgsja_scomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu,
                                         svd_paramslist[idx].jobv,
                                         svd_paramslist[idx].jobq,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].p,
                                         svd_paramslist[idx].n,
                                         svd_paramslist[idx].tola,
                                         svd_paramslist[idx].tolb );

    idx = Circular_Increment_Index(idx);
    ctgsja_obj->threshold = svd_paramslist[idx].svd_threshold;
    ctgsja_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctgsja_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctgsja_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctgsja_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CTGSJA = (Fptr_NL_LAPACKE_ctgsja)dlsym(ctgsja_obj->hModule, "LAPACKE_ctgsja");
    ASSERT_TRUE(CTGSJA != NULL) << "failed to ppt the Netlib LAPACKE_ctgsja symbol";

    CGGSVP = (Fptr_NL_LAPACKE_cggsvp)dlsym(ctgsja_obj->hModule, "LAPACKE_cggsvp3");
    ASSERT_TRUE(CGGSVP != NULL) << "failed to ppt the Netlib LAPACKE_cggsvp symbol";

    /* Preprocessing for generation of i/ps for reference and libflame apis
       by invoking the below API for both i/p sets 	*/
    ctgsja_obj->inforef = CGGSVP(   ctgsja_obj->matrix_layout,
                                    ctgsja_obj->jobu,
                                    ctgsja_obj->jobv,
                                    ctgsja_obj->jobq,
                                    ctgsja_obj->m,
                                    ctgsja_obj->p,
                                    ctgsja_obj->n,
                                    ctgsja_obj->aref,
                                    ctgsja_obj->lda,
                                    ctgsja_obj->bref,
                                    ctgsja_obj->ldb,
                                    ctgsja_obj->tola,
                                    ctgsja_obj->tolb,
                                    &ctgsja_obj->kref,
                                    &ctgsja_obj->lref,
                                    ctgsja_obj->uref,
                                    ctgsja_obj->ldu,
                                    ctgsja_obj->vref,
                                    ctgsja_obj->ldv,
                                    ctgsja_obj->qref,
                                    ctgsja_obj->ldq
                                    );

    ctgsja_obj->info = CGGSVP(   ctgsja_obj->matrix_layout,
                                    ctgsja_obj->jobu,
                                    ctgsja_obj->jobv,
                                    ctgsja_obj->jobq,
                                    ctgsja_obj->m,
                                    ctgsja_obj->p,
                                    ctgsja_obj->n,
                                    ctgsja_obj->a,
                                    ctgsja_obj->lda,
                                    ctgsja_obj->b,
                                    ctgsja_obj->ldb,
                                    ctgsja_obj->tola,
                                    ctgsja_obj->tolb,
                                    &ctgsja_obj->k,
                                    &ctgsja_obj->l,
                                    ctgsja_obj->u,
                                    ctgsja_obj->ldu,
                                    ctgsja_obj->v,
                                    ctgsja_obj->ldv,
                                    ctgsja_obj->q,
                                    ctgsja_obj->ldq
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    ctgsja_obj->inforef = CTGSJA(   ctgsja_obj->matrix_layout,
                                    ctgsja_obj->jobu,
                                    ctgsja_obj->jobv,
                                    ctgsja_obj->jobq,
                                    ctgsja_obj->m,
                                    ctgsja_obj->p,
                                    ctgsja_obj->n,
                                    ctgsja_obj->kref,
                                    ctgsja_obj->lref,
                                    ctgsja_obj->aref,
                                    ctgsja_obj->lda,
                                    ctgsja_obj->bref,
                                    ctgsja_obj->ldb,
                                    ctgsja_obj->tola,
                                    ctgsja_obj->tolb,
                                    ctgsja_obj->alpha,
                                    ctgsja_obj->beta,
                                    ctgsja_obj->uref,
                                    ctgsja_obj->ldu,
                                    ctgsja_obj->vref,
                                    ctgsja_obj->ldv,
                                    ctgsja_obj->qref,
                                    ctgsja_obj->ldq,
                                    &ctgsja_obj->ncycleref
                                    );

    /* Compute libflame's Lapacke o/p  */
    ctgsja_obj->info = LAPACKE_ctgsja(  ctgsja_obj->matrix_layout,
                                    ctgsja_obj->jobu,
                                    ctgsja_obj->jobv,
                                    ctgsja_obj->jobq,
                                    ctgsja_obj->m,
                                    ctgsja_obj->p,
                                    ctgsja_obj->n,
                                    ctgsja_obj->k,
                                    ctgsja_obj->l,
                                    ctgsja_obj->a,
                                    ctgsja_obj->lda,
                                    ctgsja_obj->b,
                                    ctgsja_obj->ldb,
                                    ctgsja_obj->tola,
                                    ctgsja_obj->tolb,
                                    ctgsja_obj->alpha,
                                    ctgsja_obj->beta,
                                    ctgsja_obj->u,
                                    ctgsja_obj->ldu,
                                    ctgsja_obj->v,
                                    ctgsja_obj->ldv,
                                    ctgsja_obj->q,
                                    ctgsja_obj->ldq,
                                    &ctgsja_obj->ncycle
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    ctgsja_obj->diff_a =  computeDiff_c( (ctgsja_obj->m)*(ctgsja_obj->n), 
                ctgsja_obj->a, ctgsja_obj->aref );

    ctgsja_obj->diff_b =  computeDiff_c( (ctgsja_obj->p)*(ctgsja_obj->n), 
                ctgsja_obj->b, ctgsja_obj->bref );

    ctgsja_obj->diff_u =  computeDiff_c( (ctgsja_obj->ldu)*(ctgsja_obj->m), 
                ctgsja_obj->u, ctgsja_obj->uref );

    ctgsja_obj->diff_v =  computeDiff_c( (ctgsja_obj->p)*(ctgsja_obj->ldv), 
                ctgsja_obj->v, ctgsja_obj->vref );

    ctgsja_obj->diff_q =  computeDiff_c( (ctgsja_obj->ldq)*(ctgsja_obj->n), 
                ctgsja_obj->q, ctgsja_obj->qref );

    ctgsja_obj->diff_alpha =  computeDiff_s( ctgsja_obj->n, 
                ctgsja_obj->alpha, ctgsja_obj->alpharef );

    ctgsja_obj->diff_beta =  computeDiff_s( ctgsja_obj->n, 
                ctgsja_obj->beta, ctgsja_obj->betaref );
				

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgsja lapack_complex_float: \n diff_a: %f \n diff_b: %f \n \
diff_u: %f \n diff_v: %f \n k: %d \n kref: %d \
diff_q: %f \n l: %d \n lref: %d \n ",
       ctgsja_obj->diff_a, ctgsja_obj->diff_b, ctgsja_obj->diff_u,
       ctgsja_obj->diff_v, ctgsja_obj->k, ctgsja_obj->kref, ctgsja_obj->diff_q,
       ctgsja_obj->l, ctgsja_obj->lref );
#endif
}

TEST_F(ctgsja_test, ctgsja1) {
    EXPECT_NEAR(0.0, ctgsja_obj->diff_a, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_b, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_q, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_u, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_v, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_alpha, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_beta, ctgsja_obj->threshold);
    EXPECT_EQ(ctgsja_obj->ncycle, ctgsja_obj->ncycleref);
}

TEST_F(ctgsja_test, ctgsja2) {
    EXPECT_NEAR(0.0, ctgsja_obj->diff_a, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_b, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_q, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_u, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_v, ctgsja_obj->threshold);
    EXPECT_EQ(ctgsja_obj->ncycle, ctgsja_obj->ncycleref);
}

TEST_F(ctgsja_test, ctgsja) {
    EXPECT_NEAR(0.0, ctgsja_obj->diff_a, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_b, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_q, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_u, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_v, ctgsja_obj->threshold);
    EXPECT_EQ(ctgsja_obj->ncycle, ctgsja_obj->ncycleref);
}

TEST_F(ctgsja_test, ctgsja4) {
    EXPECT_NEAR(0.0, ctgsja_obj->diff_a, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_b, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_q, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_u, ctgsja_obj->threshold);
    EXPECT_NEAR(0.0, ctgsja_obj->diff_v, ctgsja_obj->threshold);
    EXPECT_EQ(ctgsja_obj->ncycle, ctgsja_obj->ncycleref);
}


/* Begin tgsja_dcomplex_common_parameters  class definition */
class tgsja_dcomplex_parameters{

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

    // thresholds to determine the effective numerical rank of matrix B and a subblock of A.
    double tola;
    double tolb;
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
	lapack_int ncycle, ncycleref, *iwork, *iworkref;
	
    /*Return Values */
    int info, inforef;

   public:
      tgsja_dcomplex_parameters (int matrix_layout_i, char jobu, char jobv,
                                  char jobq, lapack_int m, lapack_int p, 
                                  lapack_int n, double tola, double tolb);
      ~tgsja_dcomplex_parameters ();
};

/* Constructor definition  tgsja lapack_complex_double_common_parameters */
tgsja_dcomplex_parameters:: tgsja_dcomplex_parameters (int matrix_layout_i,
                    char jobu_i, char jobv_i, char jobq_i, lapack_int m_i,  
            lapack_int p_i, lapack_int n_i, double tola_i, double tolb_i)
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
	diff_alpha =0;
	diff_beta = 0;
	diff_iwork = 0;
    k = kref = 0;
    l = lref = 0;
	ncycle = 0;
	ncycleref = 0;


#if LAPACKE_TEST_VERBOSE
   printf(" \n tgsja lapack_complex_double: matrix_layout: %d   jobu: %c \t \
jobv: %c \t jobq: %c  \n m: %d \t p: %d \t n: %d \n tola: %f  \t tolb: %f \n",
matrix_layout, jobu_i, jobv_i, jobq_i, m_i,
p_i, n_i, tola_i, tolb_i);
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
       EXPECT_FALSE( true) << "tgsja_dcomplex_parameters object: malloc error. Exiting ";
       tgsja_free();
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

/* Destructor definition  'tgsja_dcomplex_common_parameters' */
tgsja_dcomplex_parameters :: ~tgsja_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgsja_free();
} 

//  Test fixture class definition
class ztgsja_test  : public  ::testing::Test {
public:
   tgsja_dcomplex_parameters  *ztgsja_obj;
   void SetUp();  
   void TearDown () { delete ztgsja_obj; }
};

void ztgsja_test::SetUp()
{
    /* LAPACKE ztgsja prototype */
    typedef int (*Fptr_NL_LAPACKE_ztgsja) (  int matrix_layout, 
		char jobu, char jobv, char jobq, lapack_int m, lapack_int p,
		lapack_int n, lapack_int k, lapack_int l, lapack_complex_double* a,
		lapack_int lda, lapack_complex_double* b, lapack_int ldb, double tola,
		double tolb, double* alpha, double* beta, lapack_complex_double* u,
		lapack_int ldu, lapack_complex_double* v, lapack_int ldv,
		lapack_complex_double* q, lapack_int ldq, lapack_int* ncycle);

    Fptr_NL_LAPACKE_ztgsja ZTGSJA;

    /* LAPACKE ztgsja prototype */
    typedef int (*Fptr_NL_LAPACKE_zggsvp) ( int matrix_layout, char jobu,
	char jobv, char jobq, lapack_int m, lapack_int p, lapack_int n,
	lapack_complex_double* a, lapack_int lda, lapack_complex_double* b,
	lapack_int ldb, double tola, double tolb, lapack_int* k, lapack_int* l,
	lapack_complex_double* u, lapack_int ldu, lapack_complex_double* v,
	lapack_int ldv, lapack_complex_double* q, lapack_int ldq);

    Fptr_NL_LAPACKE_zggsvp ZGGSVP;	

    ztgsja_obj = new  tgsja_dcomplex_parameters( svd_paramslist[idx].matrix_layout,
                                         svd_paramslist[idx].jobu,
                                         svd_paramslist[idx].jobv,
                                         svd_paramslist[idx].jobq,
                                         svd_paramslist[idx].m,
                                         svd_paramslist[idx].p,
                                         svd_paramslist[idx].n,
                                         svd_paramslist[idx].tola,
                                         svd_paramslist[idx].tolb );

    idx = Circular_Increment_Index(idx);
    ztgsja_obj->threshold = svd_paramslist[idx].svd_threshold;
    ztgsja_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztgsja_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztgsja_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztgsja_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZTGSJA = (Fptr_NL_LAPACKE_ztgsja)dlsym(ztgsja_obj->hModule, "LAPACKE_ztgsja");
    ASSERT_TRUE(ZTGSJA != NULL) << "failed to ppt the Netlib LAPACKE_ztgsja symbol";

    ZGGSVP = (Fptr_NL_LAPACKE_zggsvp)dlsym(ztgsja_obj->hModule, "LAPACKE_zggsvp3");
    ASSERT_TRUE(ZGGSVP != NULL) << "failed to ppt the Netlib LAPACKE_zggsvp symbol";

    /* Preprocessing for generation of i/ps for reference and libflame apis
       by invoking the below API for both i/p sets 	*/
    ztgsja_obj->inforef = ZGGSVP(   ztgsja_obj->matrix_layout,
                                    ztgsja_obj->jobu,
                                    ztgsja_obj->jobv,
                                    ztgsja_obj->jobq,
                                    ztgsja_obj->m,
                                    ztgsja_obj->p,
                                    ztgsja_obj->n,
                                    ztgsja_obj->aref,
                                    ztgsja_obj->lda,
                                    ztgsja_obj->bref,
                                    ztgsja_obj->ldb,
                                    ztgsja_obj->tola,
                                    ztgsja_obj->tolb,
                                    &ztgsja_obj->kref,
                                    &ztgsja_obj->lref,
                                    ztgsja_obj->uref,
                                    ztgsja_obj->ldu,
                                    ztgsja_obj->vref,
                                    ztgsja_obj->ldv,
                                    ztgsja_obj->qref,
                                    ztgsja_obj->ldq
                                    );

    ztgsja_obj->info = ZGGSVP(   ztgsja_obj->matrix_layout,
                                    ztgsja_obj->jobu,
                                    ztgsja_obj->jobv,
                                    ztgsja_obj->jobq,
                                    ztgsja_obj->m,
                                    ztgsja_obj->p,
                                    ztgsja_obj->n,
                                    ztgsja_obj->a,
                                    ztgsja_obj->lda,
                                    ztgsja_obj->b,
                                    ztgsja_obj->ldb,
                                    ztgsja_obj->tola,
                                    ztgsja_obj->tolb,
                                    &ztgsja_obj->k,
                                    &ztgsja_obj->l,
                                    ztgsja_obj->u,
                                    ztgsja_obj->ldu,
                                    ztgsja_obj->v,
                                    ztgsja_obj->ldv,
                                    ztgsja_obj->q,
                                    ztgsja_obj->ldq
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    ztgsja_obj->inforef = ZTGSJA(   ztgsja_obj->matrix_layout,
                                    ztgsja_obj->jobu,
                                    ztgsja_obj->jobv,
                                    ztgsja_obj->jobq,
                                    ztgsja_obj->m,
                                    ztgsja_obj->p,
                                    ztgsja_obj->n,
                                    ztgsja_obj->kref,
                                    ztgsja_obj->lref,
                                    ztgsja_obj->aref,
                                    ztgsja_obj->lda,
                                    ztgsja_obj->bref,
                                    ztgsja_obj->ldb,
                                    ztgsja_obj->tola,
                                    ztgsja_obj->tolb,
                                    ztgsja_obj->alpha,
                                    ztgsja_obj->beta,
                                    ztgsja_obj->uref,
                                    ztgsja_obj->ldu,
                                    ztgsja_obj->vref,
                                    ztgsja_obj->ldv,
                                    ztgsja_obj->qref,
                                    ztgsja_obj->ldq,
                                    &ztgsja_obj->ncycleref
                                    );

    /* Compute libflame's Lapacke o/p  */
    ztgsja_obj->info = LAPACKE_ztgsja(  ztgsja_obj->matrix_layout,
                                    ztgsja_obj->jobu,
                                    ztgsja_obj->jobv,
                                    ztgsja_obj->jobq,
                                    ztgsja_obj->m,
                                    ztgsja_obj->p,
                                    ztgsja_obj->n,
                                    ztgsja_obj->k,
                                    ztgsja_obj->l,
                                    ztgsja_obj->a,
                                    ztgsja_obj->lda,
                                    ztgsja_obj->b,
                                    ztgsja_obj->ldb,
                                    ztgsja_obj->tola,
                                    ztgsja_obj->tolb,
                                    ztgsja_obj->alpha,
                                    ztgsja_obj->beta,
                                    ztgsja_obj->u,
                                    ztgsja_obj->ldu,
                                    ztgsja_obj->v,
                                    ztgsja_obj->ldv,
                                    ztgsja_obj->q,
                                    ztgsja_obj->ldq,
                                    &ztgsja_obj->ncycle
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    ztgsja_obj->diff_a =  computeDiff_z( (ztgsja_obj->m)*(ztgsja_obj->n), 
                ztgsja_obj->a, ztgsja_obj->aref );

    ztgsja_obj->diff_b =  computeDiff_z( (ztgsja_obj->p)*(ztgsja_obj->n), 
                ztgsja_obj->b, ztgsja_obj->bref );

    ztgsja_obj->diff_u =  computeDiff_z( (ztgsja_obj->ldu)*(ztgsja_obj->m), 
                ztgsja_obj->u, ztgsja_obj->uref );

    ztgsja_obj->diff_v =  computeDiff_z( (ztgsja_obj->p)*(ztgsja_obj->ldv), 
                ztgsja_obj->v, ztgsja_obj->vref );

    ztgsja_obj->diff_q =  computeDiff_z( (ztgsja_obj->ldq)*(ztgsja_obj->n), 
                ztgsja_obj->q, ztgsja_obj->qref );

    ztgsja_obj->diff_alpha =  computeDiff_d( ztgsja_obj->n, 
                ztgsja_obj->alpha, ztgsja_obj->alpharef );

    ztgsja_obj->diff_beta =  computeDiff_d( ztgsja_obj->n, 
                ztgsja_obj->beta, ztgsja_obj->betaref );
				

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgsja lapack_complex_double: \n diff_a: %f \n diff_b: %f \n \
diff_u: %f \n diff_v: %f \n k: %d \n kref: %d \
diff_q: %f \n l: %d \n lref: %d \n ",
       ztgsja_obj->diff_a, ztgsja_obj->diff_b, ztgsja_obj->diff_u,
       ztgsja_obj->diff_v, ztgsja_obj->k, ztgsja_obj->kref, ztgsja_obj->diff_q,
       ztgsja_obj->l, ztgsja_obj->lref );
#endif
}

TEST_F(ztgsja_test, ztgsja1) {
    EXPECT_NEAR(0.0, ztgsja_obj->diff_a, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_b, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_q, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_u, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_v, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_alpha, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_beta, ztgsja_obj->threshold);
    EXPECT_EQ(ztgsja_obj->ncycle, ztgsja_obj->ncycleref);
}

TEST_F(ztgsja_test, ztgsja2) {
    EXPECT_NEAR(0.0, ztgsja_obj->diff_a, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_b, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_q, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_u, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_v, ztgsja_obj->threshold);
    EXPECT_EQ(ztgsja_obj->ncycle, ztgsja_obj->ncycleref);
}

TEST_F(ztgsja_test, ztgsja) {
    EXPECT_NEAR(0.0, ztgsja_obj->diff_a, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_b, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_q, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_u, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_v, ztgsja_obj->threshold);
    EXPECT_EQ(ztgsja_obj->ncycle, ztgsja_obj->ncycleref);
}

TEST_F(ztgsja_test, ztgsja4) {
    EXPECT_NEAR(0.0, ztgsja_obj->diff_a, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_b, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_q, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_u, ztgsja_obj->threshold);
    EXPECT_NEAR(0.0, ztgsja_obj->diff_v, ztgsja_obj->threshold);
    EXPECT_EQ(ztgsja_obj->ncycle, ztgsja_obj->ncycleref);
}



