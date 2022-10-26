#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define tgsen_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (select!=NULL)   free(select); \
    if (selectref!=NULL)  free(selectref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (b1Q!=NULL)        free(b1Q); \
    if (b1Z!=NULL)     free(b1Z); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (d!=NULL)     free(d); \
    if (e!=NULL)     free(e); \
    if (alphar!=NULL)        free(alphar); \
    if (alpharref!=NULL)     free(alpharref); \
    if (alphai!=NULL)        free(alphai); \
    if (alphairef!=NULL)     free(alphairef); \
    if (beta!=NULL)        free(beta); \
    if (betaref!=NULL)     free(betaref); \
    if (q!=NULL)        free(q); \
    if (qref!=NULL)     free(qref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin tgsen_float_common_parameters  class definition */
class tgsen_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_q, diff_z;
    float diff_alphai, diff_alphar, diff_beta, diff_dif;
	float threshold;
    void *hModule, *dModule;
    
	// Intermediate buffers to generate q, z i/ps for tgsen API
    float* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    float* d, *e;
    float* tauq, *taup;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    lapack_int ijob; // can be 0 to 5
    char wantq; // Must be 'N', 'I', or 'V'.
    char wantz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int ifst, ifstref; // ifst and ilst mark the rows and columns of A which are to be reduced.
    lapack_int ilst, ilstref;
    
    /* Input / Output parameters */
    
    float* a, *aref; // contains the n-by-n general matrix A.
    lapack_int lda; //  The leading dimension of a
    float* b, *bref; // contains the n-by-n upper triangular matrix B.
    lapack_int ldb; //  The leading dimension of b
    float* q, *qref; // If wantq = 'V', then q is the orthogonal/unitary matrix Q1
    
    lapack_int ldq; // The leading dimension of q;
    float* z, *zref;
    lapack_int ldz; // The leading dimension of z;
    int *select, *selectref;
	float* alphar, *alphai, *beta;
	float* alpharref, *alphairef, *betaref;
	lapack_int m, mref;
	float pl, plref;
	float pr, prref;
	float dif[2], difref[2];
	
    /*Return Values */
    int info, inforef;

   public:
      tgsen_float_parameters (int matrix_layout_i, int ijob, char wantq, char wantz,
                              lapack_int n);

      ~tgsen_float_parameters ();
};

/* Constructor definition  float_common_parameters */
tgsen_float_parameters:: tgsen_float_parameters (int matrix_layout_i,
                            int ijob_i, char wantq_i, char wantz_i, 
							lapack_int n_i)
{
	int randIndex1, randIndex2;

    matrix_layout = matrix_layout_i;
    wantq = wantq_i;
    wantz = wantz_i;
    n = n_i;
	ijob = ijob_i;

    lda = n;
    ldb = n;
    ldq = n;
    ldz = n;
    
    hModule = NULL;
    dModule = NULL;
    diff_a = 0;
    diff_b = 0;
    diff_q = 0;
    diff_z = 0;

	// setting 'ifst', 'ilst' to randomized values within 'n'.
	randIndex1 = (int) (rand() %  n);
	randIndex2 = (int) (rand() %  n);
	ifst = 1;//(randIndex1 < randIndex2)? randIndex1:randIndex2;
	ifstref = ifst;
	ilst = n;//(randIndex1 < randIndex2)? randIndex2:randIndex1;
    ilstref = ilst;
	
#if LAPACKE_TEST_VERBOSE
   printf(" \n tgsen float: matrix_layout: %d n: %d  wantz: %c \
       wantq: %c  \n", matrix_layout, n, wantz, wantq);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_float_buffer_pair( &z, &zref, ldz*n );
    lapacke_gtest_alloc_float_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_float_buffer_pair( &taup, &tauq, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );
    lapacke_gtest_alloc_float_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_float_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_float_buffer_pair( &beta, &betaref, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (b1Q==NULL) || (b1Z==NULL) || \
        (z==NULL) || (zref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (select==NULL) || (selectref==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "tgsen_float_parameters object: malloc error. Exiting ";
       tgsen_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, ldb*n );
    //lapacke_gtest_init_float_buffer_pair_rand( b1Q, b1Z, ldb*n );
	//memcpy (b1Q, b, sizeof(float)*ldb*n );
    lapacke_gtest_init_float_buffer_pair_rand( q, qref, ldq*n );
    lapacke_gtest_init_float_buffer_pair_rand( z, zref, ldz*n );
	lapacke_gtest_init_float_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_float_buffer_pair_with_constant( taup, tauq, n, 0.0);
	lapacke_gtest_init_int_buffer_pair_with_constant( select, selectref, n, 0xFF);
	lapacke_gtest_init_float_buffer_pair_with_constant( alphar, alpharref, n, 0.0);
	lapacke_gtest_init_float_buffer_pair_with_constant( alphai, alphairef, n, 0.0);
	lapacke_gtest_init_float_buffer_pair_with_constant( beta, betaref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'tgsen_float_common_parameters' */
tgsen_float_parameters :: ~tgsen_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgsen_free();
} 

//  Test fixture class definition
class stgsen_test  : public  ::testing::Test {
public:
   tgsen_float_parameters  *stgsen_obj;
   void SetUp();  
   void TearDown () { delete stgsen_obj; }
};

void stgsen_test::SetUp()
{
    /* LAPACKE STGSEN prototype */
    typedef int (*Fptr_NL_LAPACKE_stgsen) ( int matrix_layout,
		lapack_int ijob, lapack_logical wantq, lapack_logical wantz,
		const lapack_logical* select, lapack_int n, float* a, lapack_int lda,
		float* b, lapack_int ldb, float* alphar, float* alphai, float* beta,
		float* q, lapack_int ldq, float* z, lapack_int ldz, lapack_int* m,
		float* pl, float* pr, float* dif );

    Fptr_NL_LAPACKE_stgsen STGSEN;

    stgsen_obj = new  tgsen_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].wantz,
                                         eig_non_sym_paramslist[idx].wantq,
										 eig_non_sym_paramslist[idx].tgsen_ijob,
                                         eig_paramslist[idx].n );
										 
    stgsen_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    stgsen_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stgsen_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stgsen_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stgsen_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    STGSEN = (Fptr_NL_LAPACKE_stgsen)dlsym(stgsen_obj->hModule, "LAPACKE_stgsen");
    ASSERT_TRUE(STGSEN != NULL) << "failed to ppt the Netlib LAPACKE_stgsen symbol";
 
    /* Transform the generalized matrices A, B to real-Schur/Schur canonical
	   form through gges API  */
	
	/* temporary variables for gges API call */
	int sdim;
	float vsl[2], vsr[2];
    stgsen_obj->inforef =  LAPACKE_sgges (   stgsen_obj->matrix_layout,
                                    'N',
                                    'N',
                                    'N',
									(LAPACK_S_SELECT3)NULL,
                                    stgsen_obj->n,
                                    stgsen_obj->a,
                                    stgsen_obj->lda,
                                    stgsen_obj->b,
                                    stgsen_obj->ldb,
                                    &sdim,
                                    stgsen_obj->alphar,
                                    stgsen_obj->alphai,
                                    stgsen_obj->beta,
									stgsen_obj->b1Q,
									stgsen_obj->n,	
									stgsen_obj->b1Z,
									stgsen_obj->n
                                    );
		
    stgsen_obj->inforef =  LAPACKE_sgges (   stgsen_obj->matrix_layout,
                                    'N',
                                    'N',
                                    'N',
									(LAPACK_S_SELECT3)NULL,
                                    stgsen_obj->n,
                                    stgsen_obj->aref,
                                    stgsen_obj->lda,
                                    stgsen_obj->bref,
                                    stgsen_obj->ldb,
                                    &sdim,
                                    stgsen_obj->alpharref,
                                    stgsen_obj->alphairef,
                                    stgsen_obj->betaref,
									NULL,
									stgsen_obj->n,
									NULL,
									stgsen_obj->n
                                    );


    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    stgsen_obj->inforef = /*STGSEN */ LAPACKE_stgsen(   stgsen_obj->matrix_layout,
                                    stgsen_obj->ijob,
                                    (lapack_logical)stgsen_obj->wantq,
                                    (lapack_logical)stgsen_obj->wantz,
                                    stgsen_obj->selectref,
                                    stgsen_obj->n,
                                    stgsen_obj->aref,
                                    stgsen_obj->lda,
                                    stgsen_obj->bref,
                                    stgsen_obj->ldb,
                                    stgsen_obj->alpharref,
                                    stgsen_obj->alphairef,
                                    stgsen_obj->betaref,
                                    stgsen_obj->qref,
                                    stgsen_obj->ldq,
                                    stgsen_obj->zref,
                                    stgsen_obj->ldz,
                                    &stgsen_obj->mref,
                                    &stgsen_obj->plref,
                                    &stgsen_obj->prref,
                                    stgsen_obj->difref
                                    );

    /* Compute libflame's Lapacke o/p  */
    stgsen_obj->info = LAPACKE_stgsen(  stgsen_obj->matrix_layout,
                                    stgsen_obj->ijob,
                                    (lapack_logical)stgsen_obj->wantq,
                                    (lapack_logical)stgsen_obj->wantz,
                                    stgsen_obj->select,
                                    stgsen_obj->n,
                                    stgsen_obj->a,
                                    stgsen_obj->lda,
                                    stgsen_obj->b,
                                    stgsen_obj->ldb,
                                    stgsen_obj->alphar,
                                    stgsen_obj->alphai,
                                    stgsen_obj->beta,
                                    stgsen_obj->q,
                                    stgsen_obj->ldq,
                                    stgsen_obj->z,
                                    stgsen_obj->ldz,
                                    &stgsen_obj->m,
                                    &stgsen_obj->pl,
                                    &stgsen_obj->pr,
                                    stgsen_obj->dif
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    stgsen_obj->diff_a =  computeDiff_s( (stgsen_obj->lda)*(stgsen_obj->n), 
                stgsen_obj->a, stgsen_obj->aref );

    stgsen_obj->diff_b =  computeDiff_s( (stgsen_obj->ldb)*(stgsen_obj->n), 
                stgsen_obj->b, stgsen_obj->bref );

    stgsen_obj->diff_q =  computeDiff_s( (stgsen_obj->ldq)*(stgsen_obj->n), 
                stgsen_obj->q, stgsen_obj->qref );

    stgsen_obj->diff_z =  computeDiff_s( (stgsen_obj->ldz)*(stgsen_obj->n), 
                stgsen_obj->z, stgsen_obj->zref );

    stgsen_obj->diff_alphar =  computeDiff_s( stgsen_obj->n, 
                stgsen_obj->alphar, stgsen_obj->alpharref );

    stgsen_obj->diff_alphai =  computeDiff_s( stgsen_obj->n, 
                stgsen_obj->alphai, stgsen_obj->alphairef );

    stgsen_obj->diff_beta =  computeDiff_s( stgsen_obj->n, 
                stgsen_obj->beta, stgsen_obj->betaref );

    stgsen_obj->diff_dif =  computeDiff_s( 2, 
                stgsen_obj->dif, stgsen_obj->difref );

}

TEST_F(stgsen_test, stgsen1) {
    EXPECT_NEAR(0.0, stgsen_obj->diff_a, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_b, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_q, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_z, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_alphar, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_alphai, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_beta, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_dif, stgsen_obj->threshold);
}

TEST_F(stgsen_test, stgsen2) {
    EXPECT_NEAR(0.0, stgsen_obj->diff_a, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_b, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_q, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_z, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_alphar, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_alphai, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_beta, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_dif, stgsen_obj->threshold);
}

TEST_F(stgsen_test, stgsen3) {
    EXPECT_NEAR(0.0, stgsen_obj->diff_a, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_b, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_q, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_z, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_alphar, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_alphai, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_beta, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_dif, stgsen_obj->threshold);
}

TEST_F(stgsen_test, stgsen4) {
    EXPECT_NEAR(0.0, stgsen_obj->diff_a, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_b, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_q, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_z, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_alphar, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_alphai, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_beta, stgsen_obj->threshold);
    EXPECT_NEAR(0.0, stgsen_obj->diff_dif, stgsen_obj->threshold);
}

