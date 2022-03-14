#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define tgevc_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (b1Q!=NULL)        free(b1Q); \
    if (b1Z!=NULL)     free(b1Z); \
    if (vl!=NULL)       free(vl); \
    if (vlref!=NULL)    free(vlref); \
    if (vr!=NULL)       free(vr); \
    if (vrref!=NULL)    free(vrref); \
    if (lscale!=NULL)   free(lscale); \
    if (lscaleref!=NULL)  free(lscaleref); \
    if (rscale!=NULL)        free(rscale); \
    if (rscaleref!=NULL)     free(rscaleref); \
    if (alphar!=NULL)     free(alphar); \
    if (alpharref!=NULL)     free(alpharref); \
    if (alphai!=NULL)     free(alphai); \
    if (alphairef!=NULL)     free(alphairef); \
    if (beta!=NULL)     free(beta); \
    if (betaref!=NULL)     free(betaref); \
    if (select!=NULL)     free(select); \
    if (selectref!=NULL)     free(selectref); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (d!=NULL)     free(d); \
    if (e!=NULL)     free(e); \
    if (taup!=NULL)     free(taup); \
    if (tauq!=NULL)     free(tauq); \
    if (q!=NULL)        free(q); \
    if (qref!=NULL)     free(qref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin tgevc_float_common_parameters  class definition */
class tgevc_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_h, diff_t, diff_q, diff_z;
    float diff_alphai, diff_alphar, diff_beta;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	//*****************************************************************
    char job; // Must be 'E' or 'S'.
	char compq; // Must be 'N', 'I', or 'V'.
	char compz; // Must be 'N', 'I', or 'V'.
	//******************************************
	char side; // Must be 'R', 'L', or 'B'.
	char howmny; //  Must be 'A', 'B', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvl;
    lapack_int ldvr;
    lapack_int mm; // The number of columns in the arrays vl and/or vr
	// Intermediate buffers to generate i/ps for hgeqz API
    float* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    float* d, *e;
    float* tauq, *taup;

    /* Input / Output parameters */
    float* a, *aref; // contains the n-by-n general matrix A.
    float* b, *bref; // contains the n-by-n upper triangular matrix B.
    float *lscale, *lscaleref;
    float *rscale, *rscaleref;
    float* q, *qref; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    float* z, *zref; //
    lapack_int ldz; // The leading dimension of z;

    /* Input / Output parameters */
    lapack_int ilo, iloref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi, ihiref;

    /* Output parameters */
    float *vl, *vlref; // left eigenvectors.
    float *vr, *vrref; // right eigenvectors
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors

    float* alphar, *alpharref;
	float* alphai, *alphairef;
	float* beta, *betaref;
    /*Return Values */
    int info, inforef;

    public:
      tgevc_float_parameters (int matrix_layout_i, char side_i, char howmny_i,
                              lapack_int n_i );
      ~tgevc_float_parameters ();
};

/* Constructor definition  float_common_parameters */
tgevc_float_parameters:: tgevc_float_parameters (int matrix_layout_i,
                            char side_i, char howmny_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
    howmny = howmny_i;
    n  = n_i;
    mm = n;
	m = 0;
	mref = 0;

    lda = n;//n+1;
    ldb = n; //n+1;
    ldq = n;
    ldz = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;
	
    diff_h = 0;
    diff_t = 0;
    diff_q = 0;
    diff_z = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
    iloref = 0;
    ihiref = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgevc float: matrix_layout: %d n: %d  job: %c \n",
                                 matrix_layout, n, job);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_float_buffer_pair( &rscale, &rscaleref, n );
    lapacke_gtest_alloc_float_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_float_buffer_pair( &z, &zref, ldz*n );
    lapacke_gtest_alloc_float_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_float_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_float_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_float_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_float_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_float_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_float_buffer_pair( &taup, &tauq, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (b1Q==NULL) || (b1Z==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (z==NULL) || (zref==NULL) || \
		(select==NULL) || (selectref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "tgevc_float_parameters object: malloc error. Exiting ";
       tgevc_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    //lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_rand( b1Q, b1Z, ldb*n );
	memcpy (b1Q, b, sizeof(float)*ldb*n );
	lapacke_gtest_init_float_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_float_buffer_pair_with_constant( taup, tauq, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( q, qref, ldq*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( z, zref, ldz*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'tgevc_float_common_parameters' */
tgevc_float_parameters :: ~tgevc_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgevc_free();
} 

//  Test fixture class definition
class stgevc_test  : public  ::testing::Test {
public:
   tgevc_float_parameters  *stgevc_obj;
   void SetUp();  
   void TearDown () { delete stgevc_obj; }
};

void stgevc_test::SetUp()
{
    /* LAPACKE_sggbal prototype */
    typedef int (*Fptr_NL_LAPACKE_sggbal) ( int matrix_layout, char job,
                    lapack_int n, float* a, lapack_int lda, float* b,
                    lapack_int ldb, lapack_int* ilo, lapack_int* ihi,
                    float* lscale, float* rscale);
                 
    Fptr_NL_LAPACKE_sggbal SGGBAL;
	
    /* LAPACKE_sgghrd prototype */	
    typedef int (*Fptr_NL_LAPACKE_sgghrd) ( int matrix_layout, char compq,
                 char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                 float* a, lapack_int lda, float* b, lapack_int ldb, 
                 float* q, lapack_int ldq, float* z, lapack_int ldz);
				 
    Fptr_NL_LAPACKE_sgghrd SGGHRD;

    /* LAPACKE_stgevc prototype */	
    typedef int (*Fptr_NL_LAPACKE_stgevc) (int matrix_layout, char side,
	char howmny, const lapack_logical* select, lapack_int n, const float* s, 
	lapack_int lds, const float* p, lapack_int ldp, float* vl, 
	lapack_int ldvl, float* vr, lapack_int ldvr, lapack_int mm, lapack_int* m);
				 
    Fptr_NL_LAPACKE_stgevc STGEVC;
    /* LAPACKE_shgeqz prototype */	
    typedef int (*Fptr_NL_LAPACKE_shgeqz) ( int matrix_layout, char job,
		char compq, char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
		float* h, lapack_int ldh, float* t, lapack_int ldt, float* alphar,
		float* alphai, float* beta, float* q, lapack_int ldq, float* z,
		lapack_int ldz );
				 
    Fptr_NL_LAPACKE_shgeqz SHGEQZ;
    stgevc_obj = new  tgevc_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].side_tgevc,
										 eig_non_sym_paramslist[idx].howmny,
                                         eig_paramslist[idx].n );
                                         
    stgevc_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    stgevc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stgevc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stgevc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stgevc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGBAL = (Fptr_NL_LAPACKE_sggbal)dlsym(stgevc_obj->hModule, "LAPACKE_sggbal");
    ASSERT_TRUE(SGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_sggbal symbol";

    SGGHRD = (Fptr_NL_LAPACKE_sgghrd)dlsym(stgevc_obj->hModule, "LAPACKE_sgghrd");
    ASSERT_TRUE(SGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_sgghrd symbol";
 
    STGEVC = (Fptr_NL_LAPACKE_stgevc)dlsym(stgevc_obj->hModule, "LAPACKE_stgevc");
    ASSERT_TRUE(STGEVC != NULL) << "failed to ppt the Netlib LAPACKE_stgevc symbol";

    SHGEQZ = (Fptr_NL_LAPACKE_shgeqz)dlsym(stgevc_obj->hModule, "LAPACKE_shgeqz");
    ASSERT_TRUE(SHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_shgeqz symbol";

	/*  Prepare the inputs to the TGEVC API from the 'gebal', 'gghrd' and 'hgeqz' 
	    api sequence  */	
	
	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    stgevc_obj->info = LAPACKE_sgebrd(  stgevc_obj->matrix_layout,
                                    stgevc_obj->n,
                                    stgevc_obj->n,
                                    stgevc_obj->b1Z,
                                    stgevc_obj->ldb,
                                    stgevc_obj->d,
                                    stgevc_obj->e,
                                    stgevc_obj->tauq,
                                    stgevc_obj->taup
                                    );

    stgevc_obj->info = LAPACKE_sorgbr(  stgevc_obj->matrix_layout,
                                    'Q',
                                    stgevc_obj->n,
                                    stgevc_obj->n,
                                    stgevc_obj->n,
                                    stgevc_obj->b1Z,
                                    stgevc_obj->ldb,
                                    stgevc_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (stgevc_obj->z, stgevc_obj->b1Z, sizeof(float)*(stgevc_obj->ldb)*(stgevc_obj->n) );
	memcpy (stgevc_obj->zref, stgevc_obj->b1Z, sizeof(float)*(stgevc_obj->ldb)*(stgevc_obj->n) );

    stgevc_obj->info = LAPACKE_sgebrd(  stgevc_obj->matrix_layout,
                                    stgevc_obj->n,
                                    stgevc_obj->n,
                                    stgevc_obj->b1Q,
                                    stgevc_obj->ldb,
                                    stgevc_obj->d,
                                    stgevc_obj->e,
                                    stgevc_obj->tauq,
                                    stgevc_obj->taup
                                    );

    stgevc_obj->info = LAPACKE_sorgbr(  stgevc_obj->matrix_layout,
                                    'Q',
                                    stgevc_obj->n,
                                    stgevc_obj->n,
                                    stgevc_obj->n,
                                    stgevc_obj->b1Q,
                                    stgevc_obj->ldb,
                                    stgevc_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (stgevc_obj->q, stgevc_obj->b1Q, sizeof(float)*(stgevc_obj->ldb)*(stgevc_obj->n) );
	memcpy (stgevc_obj->qref, stgevc_obj->b1Q, sizeof(float)*(stgevc_obj->ldb)*(stgevc_obj->n) );

    stgevc_obj->inforef = SGGBAL(   stgevc_obj->matrix_layout,
                                    'B',
                                    stgevc_obj->n,
                                    stgevc_obj->aref,
                                    stgevc_obj->lda,
                                    stgevc_obj->bref,
                                    stgevc_obj->ldb,
                                    &stgevc_obj->iloref,
                                    &stgevc_obj->ihiref,
                                    stgevc_obj->lscaleref,
                                    stgevc_obj->rscaleref
                                    );

    stgevc_obj->inforef = SGGHRD(   stgevc_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    stgevc_obj->n,
                                    stgevc_obj->iloref,
                                    stgevc_obj->ihiref, 
                                    stgevc_obj->aref,
                                    stgevc_obj->lda,
                                    stgevc_obj->bref,
                                    stgevc_obj->ldb,
                                    stgevc_obj->qref,
                                    stgevc_obj->ldq,
                                    stgevc_obj->zref,
                                    stgevc_obj->ldz
                                    );

    stgevc_obj->inforef = SHGEQZ(   stgevc_obj->matrix_layout,
                                    'S',
                                    'V',
									'V',
                                    stgevc_obj->n,
                                    stgevc_obj->iloref,
                                    stgevc_obj->ihiref, 
                                    stgevc_obj->aref,
                                    stgevc_obj->lda,
                                    stgevc_obj->bref,
                                    stgevc_obj->ldb,
									stgevc_obj->alpharref,
									stgevc_obj->alphairef,
									stgevc_obj->betaref,
                                    stgevc_obj->qref,
                                    stgevc_obj->ldq,
                                    stgevc_obj->zref,
                                    stgevc_obj->ldz
                                    );
    stgevc_obj->inforef = SGGBAL(   stgevc_obj->matrix_layout,
                                    'B',
                                    stgevc_obj->n,
                                    stgevc_obj->a,
                                    stgevc_obj->lda,
                                    stgevc_obj->b,
                                    stgevc_obj->ldb,
                                    &stgevc_obj->ilo,
                                    &stgevc_obj->ihi,
                                    stgevc_obj->lscale,
                                    stgevc_obj->rscale
                                    );

    stgevc_obj->inforef = SGGHRD(   stgevc_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    stgevc_obj->n,
                                    stgevc_obj->ilo,
                                    stgevc_obj->ihi, 
                                    stgevc_obj->a,
                                    stgevc_obj->lda,
                                    stgevc_obj->b,
                                    stgevc_obj->ldb,
                                    stgevc_obj->q,
                                    stgevc_obj->ldq,
                                    stgevc_obj->z,
                                    stgevc_obj->ldz
                                    );
    stgevc_obj->inforef = SHGEQZ(   stgevc_obj->matrix_layout,
                                    'S',
                                    'V',
									'V',
                                    stgevc_obj->n,
                                    stgevc_obj->ilo,
                                    stgevc_obj->ihi, 
                                    stgevc_obj->a,
                                    stgevc_obj->lda,
                                    stgevc_obj->b,
                                    stgevc_obj->ldb,
									stgevc_obj->alphar,
									stgevc_obj->alphai,
									stgevc_obj->beta,
                                    stgevc_obj->q,
                                    stgevc_obj->ldq,
                                    stgevc_obj->z,
                                    stgevc_obj->ldz
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    stgevc_obj->inforef = STGEVC(   stgevc_obj->matrix_layout,
                                    stgevc_obj->side,
                                    stgevc_obj->howmny,
									stgevc_obj->select,
                                    stgevc_obj->n,
                                    stgevc_obj->aref,
                                    stgevc_obj->lda,
                                    stgevc_obj->bref,
                                    stgevc_obj->ldb,
                                    stgevc_obj->qref,
                                    stgevc_obj->ldq,
                                    stgevc_obj->zref,
                                    stgevc_obj->ldz,
									stgevc_obj->mm,
									&stgevc_obj->mref
                                    );

    /* Compute libflame's Lapacke o/p  */
    stgevc_obj->inforef =  LAPACKE_stgevc (stgevc_obj->matrix_layout,
                                    stgevc_obj->side,
                                    stgevc_obj->howmny,
									stgevc_obj->select,
                                    stgevc_obj->n,
                                    stgevc_obj->a,
                                    stgevc_obj->lda,
                                    stgevc_obj->b,
                                    stgevc_obj->ldb,
                                    stgevc_obj->q,
                                    stgevc_obj->ldq,
                                    stgevc_obj->z,
                                    stgevc_obj->ldz,
									stgevc_obj->mm,
									&stgevc_obj->m
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    stgevc_obj->diff_h =  computeDiff_s( (stgevc_obj->lda)*(stgevc_obj->n), 
                stgevc_obj->a, stgevc_obj->aref );

    stgevc_obj->diff_t =  computeDiff_s( (stgevc_obj->ldb)*(stgevc_obj->n), 
                stgevc_obj->b, stgevc_obj->bref );

    stgevc_obj->diff_q =  computeDiff_s( (stgevc_obj->ldq)*(stgevc_obj->n), 
                stgevc_obj->q, stgevc_obj->qref );

    stgevc_obj->diff_z =  computeDiff_s( (stgevc_obj->ldz)*(stgevc_obj->n), 
                stgevc_obj->z, stgevc_obj->zref );

    stgevc_obj->diff_alphar =  computeDiff_s( stgevc_obj->n, 
                stgevc_obj->alphar, stgevc_obj->alpharref );

    stgevc_obj->diff_alphai =  computeDiff_s( stgevc_obj->n, 
                stgevc_obj->alphai, stgevc_obj->alphairef );

    stgevc_obj->diff_beta =  computeDiff_s( stgevc_obj->n, 
                stgevc_obj->beta, stgevc_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgevc float: \n diff_h: %f \n diff_t: %f \n \
diff_q: %f \n diff_z: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       stgevc_obj->diff_h, stgevc_obj->diff_t, stgevc_obj->diff_q,
	   stgevc_obj->diff_z, stgevc_obj->diff_alphar,
       stgevc_obj->diff_alphai, stgevc_obj->diff_beta );
#endif
}

TEST_F(stgevc_test, stgevc1) {
    EXPECT_NEAR(0.0, stgevc_obj->diff_q, stgevc_obj->threshold);
    EXPECT_NEAR(0.0, stgevc_obj->diff_z, stgevc_obj->threshold);
    EXPECT_EQ(stgevc_obj->m, stgevc_obj->mref);
}

TEST_F(stgevc_test, stgevc2) {
    EXPECT_NEAR(0.0, stgevc_obj->diff_q, stgevc_obj->threshold);
    EXPECT_NEAR(0.0, stgevc_obj->diff_z, stgevc_obj->threshold);
    EXPECT_EQ(stgevc_obj->m, stgevc_obj->mref);
}
TEST_F(stgevc_test, stgevc3) {
    EXPECT_NEAR(0.0, stgevc_obj->diff_q, stgevc_obj->threshold);
    EXPECT_NEAR(0.0, stgevc_obj->diff_z, stgevc_obj->threshold);
    EXPECT_EQ(stgevc_obj->m, stgevc_obj->mref);
}
TEST_F(stgevc_test, stgevc4) {
    EXPECT_NEAR(0.0, stgevc_obj->diff_q, stgevc_obj->threshold);
    EXPECT_NEAR(0.0, stgevc_obj->diff_z, stgevc_obj->threshold);
    EXPECT_EQ(stgevc_obj->m, stgevc_obj->mref);
}


/* Begin tgevc_double_common_parameters  class definition */
class tgevc_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_h, diff_t, diff_q, diff_z;
    double diff_alphai, diff_alphar, diff_beta;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	//*****************************************************************
    char job; // Must be 'E' or 'S'.
	char compq; // Must be 'N', 'I', or 'V'.
	char compz; // Must be 'N', 'I', or 'V'.
	//******************************************
	char side; // Must be 'R', 'L', or 'B'.
	char howmny; //  Must be 'A', 'B', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvl;
    lapack_int ldvr;
    lapack_int mm; // The number of columns in the arrays vl and/or vr
	// Intermediate buffers to generate i/ps for hgeqz API
    double* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    double* d, *e;
    double* tauq, *taup;

    /* Input / Output parameters */
    double* a, *aref; // contains the n-by-n general matrix A.
    double* b, *bref; // contains the n-by-n upper triangular matrix B.
    double *lscale, *lscaleref;
    double *rscale, *rscaleref;
    double* q, *qref; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    double* z, *zref; //
    lapack_int ldz; // The leading dimension of z;

    /* Input / Output parameters */
    lapack_int ilo, iloref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi, ihiref;

    /* Output parameters */
    double *vl, *vlref; // left eigenvectors.
    double *vr, *vrref; // right eigenvectors
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors

    double* alphar, *alpharref;
	double* alphai, *alphairef;
	double* beta, *betaref;
    /*Return Values */
    int info, inforef;

    public:
      tgevc_double_parameters (int matrix_layout_i, char side_i, char howmny_i,
                              lapack_int n_i );
      ~tgevc_double_parameters ();
};

/* Constructor definition  double_common_parameters */
tgevc_double_parameters:: tgevc_double_parameters (int matrix_layout_i,
                            char side_i, char howmny_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
    howmny = howmny_i;
    n  = n_i;
    mm = n;
	m = 0;
	mref = 0;

    lda = n;//n+1;
    ldb = n; //n+1;
    ldq = n;
    ldz = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;
	
    diff_h = 0;
    diff_t = 0;
    diff_q = 0;
    diff_z = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
    iloref = 0;
    ihiref = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgevc double: matrix_layout: %d n: %d  job: %c \n",
                                 matrix_layout, n, job);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_double_buffer_pair( &rscale, &rscaleref, n );
    lapacke_gtest_alloc_double_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_double_buffer_pair( &z, &zref, ldz*n );
    lapacke_gtest_alloc_double_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_double_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_double_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_double_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_double_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_double_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_double_buffer_pair( &taup, &tauq, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (b1Q==NULL) || (b1Z==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (z==NULL) || (zref==NULL) || \
		(select==NULL) || (selectref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "tgevc_double_parameters object: malloc error. Exiting ";
       tgevc_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, lda*n );
    //lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_double_buffer_pair_rand( b1Q, b1Z, ldb*n );
	memcpy (b1Q, b, sizeof(double)*ldb*n );
	lapacke_gtest_init_double_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_double_buffer_pair_with_constant( taup, tauq, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( q, qref, ldq*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( z, zref, ldz*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'tgevc_double_common_parameters' */
tgevc_double_parameters :: ~tgevc_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgevc_free();
} 

//  Test fixture class definition
class dtgevc_test  : public  ::testing::Test {
public:
   tgevc_double_parameters  *dtgevc_obj;
   void SetUp();  
   void TearDown () { delete dtgevc_obj; }
};

void dtgevc_test::SetUp()
{
    /* LAPACKE_dggbal prototype */
    typedef int (*Fptr_NL_LAPACKE_dggbal) ( int matrix_layout, char job,
                    lapack_int n, double* a, lapack_int lda, double* b,
                    lapack_int ldb, lapack_int* ilo, lapack_int* ihi,
                    double* lscale, double* rscale);
                 
    Fptr_NL_LAPACKE_dggbal DGGBAL;
	
    /* LAPACKE_dgghrd prototype */	
    typedef int (*Fptr_NL_LAPACKE_dgghrd) ( int matrix_layout, char compq,
                 char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                 double* a, lapack_int lda, double* b, lapack_int ldb, 
                 double* q, lapack_int ldq, double* z, lapack_int ldz);
				 
    Fptr_NL_LAPACKE_dgghrd DGGHRD;

    /* LAPACKE_dtgevc prototype */	
    typedef int (*Fptr_NL_LAPACKE_dtgevc) (int matrix_layout, char side,
	char howmny, const lapack_logical* select, lapack_int n, const double* s, 
	lapack_int lds, const double* p, lapack_int ldp, double* vl, 
	lapack_int ldvl, double* vr, lapack_int ldvr, lapack_int mm, lapack_int* m);
				 
    Fptr_NL_LAPACKE_dtgevc DTGEVC;
    /* LAPACKE_dhgeqz prototype */	
    typedef int (*Fptr_NL_LAPACKE_dhgeqz) ( int matrix_layout, char job,
		char compq, char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
		double* h, lapack_int ldh, double* t, lapack_int ldt, double* alphar,
		double* alphai, double* beta, double* q, lapack_int ldq, double* z,
		lapack_int ldz );
				 
    Fptr_NL_LAPACKE_dhgeqz DHGEQZ;
    dtgevc_obj = new  tgevc_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].side_tgevc,
										 eig_non_sym_paramslist[idx].howmny,
                                         eig_paramslist[idx].n );
                                         
    dtgevc_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    dtgevc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtgevc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtgevc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtgevc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGBAL = (Fptr_NL_LAPACKE_dggbal)dlsym(dtgevc_obj->hModule, "LAPACKE_dggbal");
    ASSERT_TRUE(DGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_dggbal symbol";

    DGGHRD = (Fptr_NL_LAPACKE_dgghrd)dlsym(dtgevc_obj->hModule, "LAPACKE_dgghrd");
    ASSERT_TRUE(DGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_dgghrd symbol";
 
    DTGEVC = (Fptr_NL_LAPACKE_dtgevc)dlsym(dtgevc_obj->hModule, "LAPACKE_dtgevc");
    ASSERT_TRUE(DTGEVC != NULL) << "failed to ppt the Netlib LAPACKE_dtgevc symbol";

    DHGEQZ = (Fptr_NL_LAPACKE_dhgeqz)dlsym(dtgevc_obj->hModule, "LAPACKE_dhgeqz");
    ASSERT_TRUE(DHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_dhgeqz symbol";

	/*  Prepare the inputs to the TGEVC API from the 'gebal', 'gghrd' and 'hgeqz' 
	    api sequence  */	
	
	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    dtgevc_obj->info = LAPACKE_dgebrd(  dtgevc_obj->matrix_layout,
                                    dtgevc_obj->n,
                                    dtgevc_obj->n,
                                    dtgevc_obj->b1Z,
                                    dtgevc_obj->ldb,
                                    dtgevc_obj->d,
                                    dtgevc_obj->e,
                                    dtgevc_obj->tauq,
                                    dtgevc_obj->taup
                                    );

    dtgevc_obj->info = LAPACKE_dorgbr(  dtgevc_obj->matrix_layout,
                                    'Q',
                                    dtgevc_obj->n,
                                    dtgevc_obj->n,
                                    dtgevc_obj->n,
                                    dtgevc_obj->b1Z,
                                    dtgevc_obj->ldb,
                                    dtgevc_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (dtgevc_obj->z, dtgevc_obj->b1Z, sizeof(double)*(dtgevc_obj->ldb)*(dtgevc_obj->n) );
	memcpy (dtgevc_obj->zref, dtgevc_obj->b1Z, sizeof(double)*(dtgevc_obj->ldb)*(dtgevc_obj->n) );

    dtgevc_obj->info = LAPACKE_dgebrd(  dtgevc_obj->matrix_layout,
                                    dtgevc_obj->n,
                                    dtgevc_obj->n,
                                    dtgevc_obj->b1Q,
                                    dtgevc_obj->ldb,
                                    dtgevc_obj->d,
                                    dtgevc_obj->e,
                                    dtgevc_obj->tauq,
                                    dtgevc_obj->taup
                                    );

    dtgevc_obj->info = LAPACKE_dorgbr(  dtgevc_obj->matrix_layout,
                                    'Q',
                                    dtgevc_obj->n,
                                    dtgevc_obj->n,
                                    dtgevc_obj->n,
                                    dtgevc_obj->b1Q,
                                    dtgevc_obj->ldb,
                                    dtgevc_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (dtgevc_obj->q, dtgevc_obj->b1Q, sizeof(double)*(dtgevc_obj->ldb)*(dtgevc_obj->n) );
	memcpy (dtgevc_obj->qref, dtgevc_obj->b1Q, sizeof(double)*(dtgevc_obj->ldb)*(dtgevc_obj->n) );

    dtgevc_obj->inforef = DGGBAL(   dtgevc_obj->matrix_layout,
                                    'B',
                                    dtgevc_obj->n,
                                    dtgevc_obj->aref,
                                    dtgevc_obj->lda,
                                    dtgevc_obj->bref,
                                    dtgevc_obj->ldb,
                                    &dtgevc_obj->iloref,
                                    &dtgevc_obj->ihiref,
                                    dtgevc_obj->lscaleref,
                                    dtgevc_obj->rscaleref
                                    );

    dtgevc_obj->inforef = DGGHRD(   dtgevc_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    dtgevc_obj->n,
                                    dtgevc_obj->iloref,
                                    dtgevc_obj->ihiref, 
                                    dtgevc_obj->aref,
                                    dtgevc_obj->lda,
                                    dtgevc_obj->bref,
                                    dtgevc_obj->ldb,
                                    dtgevc_obj->qref,
                                    dtgevc_obj->ldq,
                                    dtgevc_obj->zref,
                                    dtgevc_obj->ldz
                                    );

    dtgevc_obj->inforef = DHGEQZ(   dtgevc_obj->matrix_layout,
                                    'S',
                                    'V',
									'V',
                                    dtgevc_obj->n,
                                    dtgevc_obj->iloref,
                                    dtgevc_obj->ihiref, 
                                    dtgevc_obj->aref,
                                    dtgevc_obj->lda,
                                    dtgevc_obj->bref,
                                    dtgevc_obj->ldb,
									dtgevc_obj->alpharref,
									dtgevc_obj->alphairef,
									dtgevc_obj->betaref,
                                    dtgevc_obj->qref,
                                    dtgevc_obj->ldq,
                                    dtgevc_obj->zref,
                                    dtgevc_obj->ldz
                                    );
    dtgevc_obj->inforef = DGGBAL(   dtgevc_obj->matrix_layout,
                                    'B',
                                    dtgevc_obj->n,
                                    dtgevc_obj->a,
                                    dtgevc_obj->lda,
                                    dtgevc_obj->b,
                                    dtgevc_obj->ldb,
                                    &dtgevc_obj->ilo,
                                    &dtgevc_obj->ihi,
                                    dtgevc_obj->lscale,
                                    dtgevc_obj->rscale
                                    );

    dtgevc_obj->inforef = DGGHRD(   dtgevc_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    dtgevc_obj->n,
                                    dtgevc_obj->ilo,
                                    dtgevc_obj->ihi, 
                                    dtgevc_obj->a,
                                    dtgevc_obj->lda,
                                    dtgevc_obj->b,
                                    dtgevc_obj->ldb,
                                    dtgevc_obj->q,
                                    dtgevc_obj->ldq,
                                    dtgevc_obj->z,
                                    dtgevc_obj->ldz
                                    );
    dtgevc_obj->inforef = DHGEQZ(   dtgevc_obj->matrix_layout,
                                    'S',
                                    'V',
									'V',
                                    dtgevc_obj->n,
                                    dtgevc_obj->ilo,
                                    dtgevc_obj->ihi, 
                                    dtgevc_obj->a,
                                    dtgevc_obj->lda,
                                    dtgevc_obj->b,
                                    dtgevc_obj->ldb,
									dtgevc_obj->alphar,
									dtgevc_obj->alphai,
									dtgevc_obj->beta,
                                    dtgevc_obj->q,
                                    dtgevc_obj->ldq,
                                    dtgevc_obj->z,
                                    dtgevc_obj->ldz
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dtgevc_obj->inforef = DTGEVC(   dtgevc_obj->matrix_layout,
                                    dtgevc_obj->side,
                                    dtgevc_obj->howmny,
									dtgevc_obj->select,
                                    dtgevc_obj->n,
                                    dtgevc_obj->aref,
                                    dtgevc_obj->lda,
                                    dtgevc_obj->bref,
                                    dtgevc_obj->ldb,
                                    dtgevc_obj->qref,
                                    dtgevc_obj->ldq,
                                    dtgevc_obj->zref,
                                    dtgevc_obj->ldz,
									dtgevc_obj->mm,
									&dtgevc_obj->mref
                                    );

    /* Compute libflame's Lapacke o/p  */
    dtgevc_obj->inforef =  LAPACKE_dtgevc (dtgevc_obj->matrix_layout,
                                    dtgevc_obj->side,
                                    dtgevc_obj->howmny,
									dtgevc_obj->select,
                                    dtgevc_obj->n,
                                    dtgevc_obj->a,
                                    dtgevc_obj->lda,
                                    dtgevc_obj->b,
                                    dtgevc_obj->ldb,
                                    dtgevc_obj->q,
                                    dtgevc_obj->ldq,
                                    dtgevc_obj->z,
                                    dtgevc_obj->ldz,
									dtgevc_obj->mm,
									&dtgevc_obj->m
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    dtgevc_obj->diff_h =  computeDiff_d( (dtgevc_obj->lda)*(dtgevc_obj->n), 
                dtgevc_obj->a, dtgevc_obj->aref );

    dtgevc_obj->diff_t =  computeDiff_d( (dtgevc_obj->ldb)*(dtgevc_obj->n), 
                dtgevc_obj->b, dtgevc_obj->bref );

    dtgevc_obj->diff_q =  computeDiff_d( (dtgevc_obj->ldq)*(dtgevc_obj->n), 
                dtgevc_obj->q, dtgevc_obj->qref );

    dtgevc_obj->diff_z =  computeDiff_d( (dtgevc_obj->ldz)*(dtgevc_obj->n), 
                dtgevc_obj->z, dtgevc_obj->zref );

    dtgevc_obj->diff_alphar =  computeDiff_d( dtgevc_obj->n, 
                dtgevc_obj->alphar, dtgevc_obj->alpharref );

    dtgevc_obj->diff_alphai =  computeDiff_d( dtgevc_obj->n, 
                dtgevc_obj->alphai, dtgevc_obj->alphairef );

    dtgevc_obj->diff_beta =  computeDiff_d( dtgevc_obj->n, 
                dtgevc_obj->beta, dtgevc_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgevc double: \n diff_h: %f \n diff_t: %f \n \
diff_q: %f \n diff_z: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       dtgevc_obj->diff_h, dtgevc_obj->diff_t, dtgevc_obj->diff_q,
	   dtgevc_obj->diff_z, dtgevc_obj->diff_alphar,
       dtgevc_obj->diff_alphai, dtgevc_obj->diff_beta );
#endif
}

TEST_F(dtgevc_test, dtgevc1) {
    EXPECT_NEAR(0.0, dtgevc_obj->diff_q, dtgevc_obj->threshold);
    EXPECT_NEAR(0.0, dtgevc_obj->diff_z, dtgevc_obj->threshold);
    EXPECT_EQ(dtgevc_obj->m, dtgevc_obj->mref);
}

TEST_F(dtgevc_test, dtgevc2) {
    EXPECT_NEAR(0.0, dtgevc_obj->diff_q, dtgevc_obj->threshold);
    EXPECT_NEAR(0.0, dtgevc_obj->diff_z, dtgevc_obj->threshold);
    EXPECT_EQ(dtgevc_obj->m, dtgevc_obj->mref);
}
TEST_F(dtgevc_test, dtgevc3) {
    EXPECT_NEAR(0.0, dtgevc_obj->diff_q, dtgevc_obj->threshold);
    EXPECT_NEAR(0.0, dtgevc_obj->diff_z, dtgevc_obj->threshold);
    EXPECT_EQ(dtgevc_obj->m, dtgevc_obj->mref);
}
TEST_F(dtgevc_test, dtgevc4) {
    EXPECT_NEAR(0.0, dtgevc_obj->diff_q, dtgevc_obj->threshold);
    EXPECT_NEAR(0.0, dtgevc_obj->diff_z, dtgevc_obj->threshold);
    EXPECT_EQ(dtgevc_obj->m, dtgevc_obj->mref);
}

/* Begin tgevc_scomplex_common_parameters  class definition */
class tgevc_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_h, diff_t, diff_q, diff_z;
    float diff_alphai, diff_alphar, diff_beta;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	//*****************************************************************
    char job; // Must be 'E' or 'S'.
	char compq; // Must be 'N', 'I', or 'V'.
	char compz; // Must be 'N', 'I', or 'V'.
	//******************************************
	char side; // Must be 'R', 'L', or 'B'.
	char howmny; //  Must be 'A', 'B', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvl;
    lapack_int ldvr;
    lapack_int mm; // The number of columns in the arrays vl and/or vr
	// Intermediate buffers to generate i/ps for hgeqz API
    lapack_complex_float* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    float* d, *e;
    lapack_complex_float* tauq, *taup;

    /* Input / Output parameters */
    lapack_complex_float* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_float* b, *bref; // contains the n-by-n upper triangular matrix B.
    float *lscale, *lscaleref;
    float *rscale, *rscaleref;
    lapack_complex_float* q, *qref; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    lapack_complex_float* z, *zref; //
    lapack_int ldz; // The leading dimension of z;

    /* Input / Output parameters */
    lapack_int ilo, iloref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi, ihiref;

    /* Output parameters */
    lapack_complex_float *vl, *vlref; // left eigenvectors.
    lapack_complex_float *vr, *vrref; // right eigenvectors
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors

    lapack_complex_float* alphar, *alpharref;
	lapack_complex_float* alphai, *alphairef;
	lapack_complex_float* beta, *betaref;
    /*Return Values */
    int info, inforef;

    public:
      tgevc_scomplex_parameters (int matrix_layout_i, char side_i, char howmny_i,
                              lapack_int n_i );
      ~tgevc_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
tgevc_scomplex_parameters:: tgevc_scomplex_parameters (int matrix_layout_i,
                            char side_i, char howmny_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
    howmny = howmny_i;
    n  = n_i;
    mm = n;
	m = 0;
	mref = 0;

    lda = n;//n+1;
    ldb = n; //n+1;
    ldq = n;
    ldz = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;
	
    diff_h = 0;
    diff_t = 0;
    diff_q = 0;
    diff_z = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
    iloref = 0;
    ihiref = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgevc lapack_complex_float: matrix_layout: %d n: %d  job: %c \n",
                                 matrix_layout, n, job);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_float_buffer_pair( &rscale, &rscaleref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &z, &zref, ldz*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_float_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &taup, &tauq, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (b1Q==NULL) || (b1Z==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (z==NULL) || (zref==NULL) || \
		(select==NULL) || (selectref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "tgevc_scomplex_parameters object: malloc error. Exiting ";
       tgevc_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    //lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( b1Q, b1Z, ldb*n );
	memcpy (b1Q, b, sizeof(lapack_complex_float)*ldb*n );
	lapacke_gtest_init_float_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_scomplex_buffer_pair_with_constant( taup, tauq, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( q, qref, ldq*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( z, zref, ldz*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'tgevc_scomplex_common_parameters' */
tgevc_scomplex_parameters :: ~tgevc_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgevc_free();
} 

//  Test fixture class definition
class ctgevc_test  : public  ::testing::Test {
public:
   tgevc_scomplex_parameters  *ctgevc_obj;
   void SetUp();  
   void TearDown () { delete ctgevc_obj; }
};

void ctgevc_test::SetUp()
{
    /* LAPACKE_cggbal prototype */
    typedef int (*Fptr_NL_LAPACKE_cggbal) ( int matrix_layout, char job,
                    lapack_int n, lapack_complex_float* a, lapack_int lda, 
					lapack_complex_float* b, lapack_int ldb,
					lapack_int* ilo, lapack_int* ihi,
                    float* lscale, float* rscale);
                 
    Fptr_NL_LAPACKE_cggbal CGGBAL;
	
    /* LAPACKE_cgghrd prototype */	
    typedef int (*Fptr_NL_LAPACKE_cgghrd) ( int matrix_layout, char compq,
                 char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                 lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, 
				 lapack_int ldb, lapack_complex_float* q, lapack_int ldq, 
				 lapack_complex_float* z, lapack_int ldz);
				 
    Fptr_NL_LAPACKE_cgghrd CGGHRD;

    /* LAPACKE_ctgevc prototype */	
    typedef int (*Fptr_NL_LAPACKE_ctgevc) (int matrix_layout, char side,
		char howmny, const lapack_logical* select, lapack_int n,
		const lapack_complex_float* s, lapack_int lds, const lapack_complex_float* p,
		lapack_int ldp, lapack_complex_float* vl, lapack_int ldvl,
		lapack_complex_float* vr, lapack_int ldvr, lapack_int mm, lapack_int* m);
				 
    Fptr_NL_LAPACKE_ctgevc CTGEVC;

    /* LAPACKE_chgeqz prototype */	
    typedef int (*Fptr_NL_LAPACKE_chgeqz) ( int matrix_layout, char job,
		char compq, char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
		lapack_complex_float* h, lapack_int ldh, lapack_complex_float* t, lapack_int ldt, lapack_complex_float* alpha,
		 lapack_complex_float* beta, lapack_complex_float* q, lapack_int ldq, lapack_complex_float* z,
		lapack_int ldz );
				 
    Fptr_NL_LAPACKE_chgeqz CHGEQZ;
    ctgevc_obj = new  tgevc_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].side_tgevc,
										 eig_non_sym_paramslist[idx].howmny,
                                         eig_paramslist[idx].n );
                                         
    ctgevc_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    ctgevc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctgevc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctgevc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctgevc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGBAL = (Fptr_NL_LAPACKE_cggbal)dlsym(ctgevc_obj->hModule, "LAPACKE_cggbal");
    ASSERT_TRUE(CGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_cggbal symbol";

    CGGHRD = (Fptr_NL_LAPACKE_cgghrd)dlsym(ctgevc_obj->hModule, "LAPACKE_cgghrd");
    ASSERT_TRUE(CGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_cgghrd symbol";
 
    CTGEVC = (Fptr_NL_LAPACKE_ctgevc)dlsym(ctgevc_obj->hModule, "LAPACKE_ctgevc");
    ASSERT_TRUE(CTGEVC != NULL) << "failed to ppt the Netlib LAPACKE_ctgevc symbol";

    CHGEQZ = (Fptr_NL_LAPACKE_chgeqz)dlsym(ctgevc_obj->hModule, "LAPACKE_chgeqz");
    ASSERT_TRUE(CHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_chgeqz symbol";

	/*  Prepare the inputs to the TGEVC API from the 'gebal', 'gghrd' and 'hgeqz' 
	    api sequence  */	
	
	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    ctgevc_obj->info = LAPACKE_cgebrd(  ctgevc_obj->matrix_layout,
                                    ctgevc_obj->n,
                                    ctgevc_obj->n,
                                    ctgevc_obj->b1Z,
                                    ctgevc_obj->ldb,
                                    ctgevc_obj->d,
                                    ctgevc_obj->e,
                                    ctgevc_obj->tauq,
                                    ctgevc_obj->taup
                                    );

    ctgevc_obj->info = LAPACKE_cungbr(  ctgevc_obj->matrix_layout,
                                    'Q',
                                    ctgevc_obj->n,
                                    ctgevc_obj->n,
                                    ctgevc_obj->n,
                                    ctgevc_obj->b1Z,
                                    ctgevc_obj->ldb,
                                    ctgevc_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (ctgevc_obj->z, ctgevc_obj->b1Z, sizeof(float)*(ctgevc_obj->ldb)*(ctgevc_obj->n) );
	memcpy (ctgevc_obj->zref, ctgevc_obj->b1Z, sizeof(float)*(ctgevc_obj->ldb)*(ctgevc_obj->n) );

    ctgevc_obj->info = LAPACKE_cgebrd(  ctgevc_obj->matrix_layout,
                                    ctgevc_obj->n,
                                    ctgevc_obj->n,
                                    ctgevc_obj->b1Q,
                                    ctgevc_obj->ldb,
                                    ctgevc_obj->d,
                                    ctgevc_obj->e,
                                    ctgevc_obj->tauq,
                                    ctgevc_obj->taup
                                    );

    ctgevc_obj->info = LAPACKE_cungbr(  ctgevc_obj->matrix_layout,
                                    'Q',
                                    ctgevc_obj->n,
                                    ctgevc_obj->n,
                                    ctgevc_obj->n,
                                    ctgevc_obj->b1Q,
                                    ctgevc_obj->ldb,
                                    ctgevc_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (ctgevc_obj->q, ctgevc_obj->b1Q, sizeof(lapack_complex_float)*(ctgevc_obj->ldb)*(ctgevc_obj->n) );
	memcpy (ctgevc_obj->qref, ctgevc_obj->b1Q, sizeof(lapack_complex_float)*(ctgevc_obj->ldb)*(ctgevc_obj->n) );

    ctgevc_obj->inforef = CGGBAL(   ctgevc_obj->matrix_layout,
                                    'B',
                                    ctgevc_obj->n,
                                    ctgevc_obj->aref,
                                    ctgevc_obj->lda,
                                    ctgevc_obj->bref,
                                    ctgevc_obj->ldb,
                                    &ctgevc_obj->iloref,
                                    &ctgevc_obj->ihiref,
                                    ctgevc_obj->lscaleref,
                                    ctgevc_obj->rscaleref
                                    );

    ctgevc_obj->inforef = CGGHRD(   ctgevc_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    ctgevc_obj->n,
                                    ctgevc_obj->iloref,
                                    ctgevc_obj->ihiref, 
                                    ctgevc_obj->aref,
                                    ctgevc_obj->lda,
                                    ctgevc_obj->bref,
                                    ctgevc_obj->ldb,
                                    ctgevc_obj->qref,
                                    ctgevc_obj->ldq,
                                    ctgevc_obj->zref,
                                    ctgevc_obj->ldz
                                    );

    ctgevc_obj->inforef = CHGEQZ(   ctgevc_obj->matrix_layout,
                                    'S',
                                    'V',
									'V',
                                    ctgevc_obj->n,
                                    ctgevc_obj->iloref,
                                    ctgevc_obj->ihiref, 
                                    ctgevc_obj->aref,
                                    ctgevc_obj->lda,
                                    ctgevc_obj->bref,
                                    ctgevc_obj->ldb,
									ctgevc_obj->alpharref,
									ctgevc_obj->betaref,
                                    ctgevc_obj->qref,
                                    ctgevc_obj->ldq,
                                    ctgevc_obj->zref,
                                    ctgevc_obj->ldz
                                    );
    ctgevc_obj->inforef = CGGBAL(   ctgevc_obj->matrix_layout,
                                    'B',
                                    ctgevc_obj->n,
                                    ctgevc_obj->a,
                                    ctgevc_obj->lda,
                                    ctgevc_obj->b,
                                    ctgevc_obj->ldb,
                                    &ctgevc_obj->ilo,
                                    &ctgevc_obj->ihi,
                                    ctgevc_obj->lscale,
                                    ctgevc_obj->rscale
                                    );

    ctgevc_obj->inforef = CGGHRD(   ctgevc_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    ctgevc_obj->n,
                                    ctgevc_obj->ilo,
                                    ctgevc_obj->ihi, 
                                    ctgevc_obj->a,
                                    ctgevc_obj->lda,
                                    ctgevc_obj->b,
                                    ctgevc_obj->ldb,
                                    ctgevc_obj->q,
                                    ctgevc_obj->ldq,
                                    ctgevc_obj->z,
                                    ctgevc_obj->ldz
                                    );
    ctgevc_obj->inforef = CHGEQZ(   ctgevc_obj->matrix_layout,
                                    'S',
                                    'V',
									'V',
                                    ctgevc_obj->n,
                                    ctgevc_obj->ilo,
                                    ctgevc_obj->ihi, 
                                    ctgevc_obj->a,
                                    ctgevc_obj->lda,
                                    ctgevc_obj->b,
                                    ctgevc_obj->ldb,
									ctgevc_obj->alphar,
									ctgevc_obj->beta,
                                    ctgevc_obj->q,
                                    ctgevc_obj->ldq,
                                    ctgevc_obj->z,
                                    ctgevc_obj->ldz
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    ctgevc_obj->inforef = CTGEVC(   ctgevc_obj->matrix_layout,
                                    ctgevc_obj->side,
                                    ctgevc_obj->howmny,
									ctgevc_obj->select,
                                    ctgevc_obj->n,
                                    ctgevc_obj->aref,
                                    ctgevc_obj->lda,
                                    ctgevc_obj->bref,
                                    ctgevc_obj->ldb,
                                    ctgevc_obj->qref,
                                    ctgevc_obj->ldq,
                                    ctgevc_obj->zref,
                                    ctgevc_obj->ldz,
									ctgevc_obj->mm,
									&ctgevc_obj->mref
                                    );

    /* Compute libflame's Lapacke o/p  */
    ctgevc_obj->inforef =  LAPACKE_ctgevc (ctgevc_obj->matrix_layout,
                                    ctgevc_obj->side,
                                    ctgevc_obj->howmny,
									ctgevc_obj->select,
                                    ctgevc_obj->n,
                                    ctgevc_obj->a,
                                    ctgevc_obj->lda,
                                    ctgevc_obj->b,
                                    ctgevc_obj->ldb,
                                    ctgevc_obj->q,
                                    ctgevc_obj->ldq,
                                    ctgevc_obj->z,
                                    ctgevc_obj->ldz,
									ctgevc_obj->mm,
									&ctgevc_obj->m
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    ctgevc_obj->diff_h =  computeDiff_c( (ctgevc_obj->lda)*(ctgevc_obj->n), 
                ctgevc_obj->a, ctgevc_obj->aref );

    ctgevc_obj->diff_t =  computeDiff_c( (ctgevc_obj->ldb)*(ctgevc_obj->n), 
                ctgevc_obj->b, ctgevc_obj->bref );

    ctgevc_obj->diff_q =  computeDiff_c( (ctgevc_obj->ldq)*(ctgevc_obj->n), 
                ctgevc_obj->q, ctgevc_obj->qref );

    ctgevc_obj->diff_z =  computeDiff_c( (ctgevc_obj->ldz)*(ctgevc_obj->n), 
                ctgevc_obj->z, ctgevc_obj->zref );

    ctgevc_obj->diff_alphar =  computeDiff_c( ctgevc_obj->n, 
                ctgevc_obj->alphar, ctgevc_obj->alpharref );

    ctgevc_obj->diff_alphai =  computeDiff_c( ctgevc_obj->n, 
                ctgevc_obj->alphai, ctgevc_obj->alphairef );

    ctgevc_obj->diff_beta =  computeDiff_c( ctgevc_obj->n, 
                ctgevc_obj->beta, ctgevc_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgevc lapack_complex_float: \n diff_h: %f \n diff_t: %f \n \
diff_q: %f \n diff_z: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       ctgevc_obj->diff_h, ctgevc_obj->diff_t, ctgevc_obj->diff_q,
	   ctgevc_obj->diff_z, ctgevc_obj->diff_alphar,
       ctgevc_obj->diff_alphai, ctgevc_obj->diff_beta );
#endif
}

TEST_F(ctgevc_test, ctgevc1) {
    EXPECT_NEAR(0.0, ctgevc_obj->diff_q, ctgevc_obj->threshold);
    EXPECT_NEAR(0.0, ctgevc_obj->diff_z, ctgevc_obj->threshold);
    EXPECT_EQ(ctgevc_obj->m, ctgevc_obj->mref);
}

TEST_F(ctgevc_test, ctgevc2) {
    EXPECT_NEAR(0.0, ctgevc_obj->diff_q, ctgevc_obj->threshold);
    EXPECT_NEAR(0.0, ctgevc_obj->diff_z, ctgevc_obj->threshold);
    EXPECT_EQ(ctgevc_obj->m, ctgevc_obj->mref);
}
TEST_F(ctgevc_test, ctgevc3) {
    EXPECT_NEAR(0.0, ctgevc_obj->diff_q, ctgevc_obj->threshold);
    EXPECT_NEAR(0.0, ctgevc_obj->diff_z, ctgevc_obj->threshold);
    EXPECT_EQ(ctgevc_obj->m, ctgevc_obj->mref);
}
TEST_F(ctgevc_test, ctgevc4) {
    EXPECT_NEAR(0.0, ctgevc_obj->diff_q, ctgevc_obj->threshold);
    EXPECT_NEAR(0.0, ctgevc_obj->diff_z, ctgevc_obj->threshold);
    EXPECT_EQ(ctgevc_obj->m, ctgevc_obj->mref);
}

/* Begin tgevc_dcomplex_common_parameters  class definition */
class tgevc_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_h, diff_t, diff_q, diff_z;
    double diff_alphai, diff_alphar, diff_beta;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	//*****************************************************************
    char job; // Must be 'E' or 'S'.
	char compq; // Must be 'N', 'I', or 'V'.
	char compz; // Must be 'N', 'I', or 'V'.
	//******************************************
	char side; // Must be 'R', 'L', or 'B'.
	char howmny; //  Must be 'A', 'B', or 'S'.
	int *select, *selectref; // Specifies the eigenvectors to be computed.  

    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
    lapack_int ldvl;
    lapack_int ldvr;
    lapack_int mm; // The number of columns in the arrays vl and/or vr
	// Intermediate buffers to generate i/ps for hgeqz API
    lapack_complex_double* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    double* d, *e;
    lapack_complex_double* tauq, *taup;

    /* Input / Output parameters */
    lapack_complex_double* a, *aref; // contains the n-by-n general matrix A.
    lapack_complex_double* b, *bref; // contains the n-by-n upper triangular matrix B.
    double *lscale, *lscaleref;
    double *rscale, *rscaleref;
    lapack_complex_double* q, *qref; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    lapack_complex_double* z, *zref; //
    lapack_int ldz; // The leading dimension of z;

    /* Input / Output parameters */
    lapack_int ilo, iloref; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi, ihiref;

    /* Output parameters */
    lapack_complex_double *vl, *vlref; // left eigenvectors.
    lapack_complex_double *vr, *vrref; // right eigenvectors
    lapack_int m, mref; // The number of columns in the arrays VL and/or VR required to store the eigenvectors

    lapack_complex_double* alphar, *alpharref;
	lapack_complex_double* alphai, *alphairef;
	lapack_complex_double* beta, *betaref;
    /*Return Values */
    int info, inforef;

    public:
      tgevc_dcomplex_parameters (int matrix_layout_i, char side_i, char howmny_i,
                              lapack_int n_i );
      ~tgevc_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
tgevc_dcomplex_parameters:: tgevc_dcomplex_parameters (int matrix_layout_i,
                            char side_i, char howmny_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    side = side_i;
    howmny = howmny_i;
    n  = n_i;
    mm = n;
	m = 0;
	mref = 0;

    lda = n;//n+1;
    ldb = n; //n+1;
    ldq = n;
    ldz = n;
    ldvl = n;
    ldvr = n;

    hModule = NULL;
    dModule = NULL;
	
    diff_h = 0;
    diff_t = 0;
    diff_q = 0;
    diff_z = 0;
    diff_alphai = 0;
    diff_alphar = 0;
    diff_beta = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
    iloref = 0;
    ihiref = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgevc lapack_complex_double: matrix_layout: %d n: %d  job: %c \n",
                                 matrix_layout, n, job);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &lscale, &lscaleref, n );
    lapacke_gtest_alloc_double_buffer_pair( &rscale, &rscaleref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &z, &zref, ldz*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &alphar, &alpharref, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &alphai, &alphairef, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &beta, &betaref, n );
    lapacke_gtest_alloc_int_buffer_pair( &select, &selectref, n );	
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vl, &vlref, n*ldvl );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vr, &vrref, n*ldvr );
    lapacke_gtest_alloc_double_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &taup, &tauq, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (b1Q==NULL) || (b1Z==NULL) || \
        (lscale==NULL) || (lscaleref==NULL) || \
        (rscale==NULL) || (rscaleref==NULL) || \
        (alphar==NULL) || (alpharref==NULL) || \
        (alphai==NULL) || (alphairef==NULL) || \
        (beta==NULL) || (betaref==NULL) || \
        (vl==NULL) || (vlref==NULL) || \
        (vr==NULL) || (vrref==NULL) || \
        (z==NULL) || (zref==NULL) || \
		(select==NULL) || (selectref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "tgevc_dcomplex_parameters object: malloc error. Exiting ";
       tgevc_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    //lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( b, bref, n,ldb, 'U');
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b1Q, b1Z, ldb*n );
	memcpy (b1Q, b, sizeof(lapack_complex_double)*ldb*n );
	lapacke_gtest_init_double_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_dcomplex_buffer_pair_with_constant( taup, tauq, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( lscale, lscaleref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rscale, rscaleref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( q, qref, ldq*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( z, zref, ldz*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( alphar, alpharref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( alphai, alphairef, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( beta, betaref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(select, selectref, n, -1);

   } /* end of Constructor  */
    

/* Destructor definition  'tgevc_dcomplex_common_parameters' */
tgevc_dcomplex_parameters :: ~tgevc_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgevc_free();
} 

//  Test fixture class definition
class ztgevc_test  : public  ::testing::Test {
public:
   tgevc_dcomplex_parameters  *ztgevc_obj;
   void SetUp();  
   void TearDown () { delete ztgevc_obj; }
};

void ztgevc_test::SetUp()
{
    /* LAPACKE_zggbal prototype */
    typedef int (*Fptr_NL_LAPACKE_zggbal) ( int matrix_layout, char job,
                    lapack_int n, lapack_complex_double* a, lapack_int lda, 
					lapack_complex_double* b, lapack_int ldb,
					lapack_int* ilo, lapack_int* ihi,
                    double* lscale, double* rscale);
                 
    Fptr_NL_LAPACKE_zggbal ZGGBAL;
	
    /* LAPACKE_zgghrd prototype */	
    typedef int (*Fptr_NL_LAPACKE_zgghrd) ( int matrix_layout, char compq,
                 char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                 lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, 
				 lapack_int ldb, lapack_complex_double* q, lapack_int ldq, 
				 lapack_complex_double* z, lapack_int ldz);
				 
    Fptr_NL_LAPACKE_zgghrd ZGGHRD;

    /* LAPACKE_ztgevc prototype */	
    typedef int (*Fptr_NL_LAPACKE_ztgevc) (int matrix_layout, char side,
		char howmny, const lapack_logical* select, lapack_int n,
		const lapack_complex_double* s, lapack_int lds, const lapack_complex_double* p,
		lapack_int ldp, lapack_complex_double* vl, lapack_int ldvl,
		lapack_complex_double* vr, lapack_int ldvr, lapack_int mm, lapack_int* m);
				 
    Fptr_NL_LAPACKE_ztgevc ZTGEVC;

    /* LAPACKE_zhgeqz prototype */	
    typedef int (*Fptr_NL_LAPACKE_zhgeqz) ( int matrix_layout, char job,
		char compq, char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
		lapack_complex_double* h, lapack_int ldh, lapack_complex_double* t, lapack_int ldt, lapack_complex_double* alpha,
		 lapack_complex_double* beta, lapack_complex_double* q, lapack_int ldq, lapack_complex_double* z,
		lapack_int ldz );
				 
    Fptr_NL_LAPACKE_zhgeqz ZHGEQZ;
    ztgevc_obj = new  tgevc_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].side_tgevc,
										 eig_non_sym_paramslist[idx].howmny,
                                         eig_paramslist[idx].n );
                                         
    ztgevc_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    ztgevc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztgevc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztgevc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztgevc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGBAL = (Fptr_NL_LAPACKE_zggbal)dlsym(ztgevc_obj->hModule, "LAPACKE_zggbal");
    ASSERT_TRUE(ZGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_zggbal symbol";

    ZGGHRD = (Fptr_NL_LAPACKE_zgghrd)dlsym(ztgevc_obj->hModule, "LAPACKE_zgghrd");
    ASSERT_TRUE(ZGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_zgghrd symbol";
 
    ZTGEVC = (Fptr_NL_LAPACKE_ztgevc)dlsym(ztgevc_obj->hModule, "LAPACKE_ztgevc");
    ASSERT_TRUE(ZTGEVC != NULL) << "failed to ppt the Netlib LAPACKE_ztgevc symbol";

    ZHGEQZ = (Fptr_NL_LAPACKE_zhgeqz)dlsym(ztgevc_obj->hModule, "LAPACKE_zhgeqz");
    ASSERT_TRUE(ZHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_zhgeqz symbol";

	/*  Prepare the inputs to the TGEVC API from the 'gebal', 'gghrd' and 'hgeqz' 
	    api sequence  */	
	
	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    ztgevc_obj->info = LAPACKE_zgebrd(  ztgevc_obj->matrix_layout,
                                    ztgevc_obj->n,
                                    ztgevc_obj->n,
                                    ztgevc_obj->b1Z,
                                    ztgevc_obj->ldb,
                                    ztgevc_obj->d,
                                    ztgevc_obj->e,
                                    ztgevc_obj->tauq,
                                    ztgevc_obj->taup
                                    );

    ztgevc_obj->info = LAPACKE_zungbr(  ztgevc_obj->matrix_layout,
                                    'Q',
                                    ztgevc_obj->n,
                                    ztgevc_obj->n,
                                    ztgevc_obj->n,
                                    ztgevc_obj->b1Z,
                                    ztgevc_obj->ldb,
                                    ztgevc_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (ztgevc_obj->z, ztgevc_obj->b1Z, sizeof(double)*(ztgevc_obj->ldb)*(ztgevc_obj->n) );
	memcpy (ztgevc_obj->zref, ztgevc_obj->b1Z, sizeof(double)*(ztgevc_obj->ldb)*(ztgevc_obj->n) );

    ztgevc_obj->info = LAPACKE_zgebrd(  ztgevc_obj->matrix_layout,
                                    ztgevc_obj->n,
                                    ztgevc_obj->n,
                                    ztgevc_obj->b1Q,
                                    ztgevc_obj->ldb,
                                    ztgevc_obj->d,
                                    ztgevc_obj->e,
                                    ztgevc_obj->tauq,
                                    ztgevc_obj->taup
                                    );

    ztgevc_obj->info = LAPACKE_zungbr(  ztgevc_obj->matrix_layout,
                                    'Q',
                                    ztgevc_obj->n,
                                    ztgevc_obj->n,
                                    ztgevc_obj->n,
                                    ztgevc_obj->b1Q,
                                    ztgevc_obj->ldb,
                                    ztgevc_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (ztgevc_obj->q, ztgevc_obj->b1Q, sizeof(lapack_complex_double)*(ztgevc_obj->ldb)*(ztgevc_obj->n) );
	memcpy (ztgevc_obj->qref, ztgevc_obj->b1Q, sizeof(lapack_complex_double)*(ztgevc_obj->ldb)*(ztgevc_obj->n) );

    ztgevc_obj->inforef = ZGGBAL(   ztgevc_obj->matrix_layout,
                                    'B',
                                    ztgevc_obj->n,
                                    ztgevc_obj->aref,
                                    ztgevc_obj->lda,
                                    ztgevc_obj->bref,
                                    ztgevc_obj->ldb,
                                    &ztgevc_obj->iloref,
                                    &ztgevc_obj->ihiref,
                                    ztgevc_obj->lscaleref,
                                    ztgevc_obj->rscaleref
                                    );

    ztgevc_obj->inforef = ZGGHRD(   ztgevc_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    ztgevc_obj->n,
                                    ztgevc_obj->iloref,
                                    ztgevc_obj->ihiref, 
                                    ztgevc_obj->aref,
                                    ztgevc_obj->lda,
                                    ztgevc_obj->bref,
                                    ztgevc_obj->ldb,
                                    ztgevc_obj->qref,
                                    ztgevc_obj->ldq,
                                    ztgevc_obj->zref,
                                    ztgevc_obj->ldz
                                    );

    ztgevc_obj->inforef = ZHGEQZ(   ztgevc_obj->matrix_layout,
                                    'S',
                                    'V',
									'V',
                                    ztgevc_obj->n,
                                    ztgevc_obj->iloref,
                                    ztgevc_obj->ihiref, 
                                    ztgevc_obj->aref,
                                    ztgevc_obj->lda,
                                    ztgevc_obj->bref,
                                    ztgevc_obj->ldb,
									ztgevc_obj->alpharref,
									ztgevc_obj->betaref,
                                    ztgevc_obj->qref,
                                    ztgevc_obj->ldq,
                                    ztgevc_obj->zref,
                                    ztgevc_obj->ldz
                                    );
    ztgevc_obj->inforef = ZGGBAL(   ztgevc_obj->matrix_layout,
                                    'B',
                                    ztgevc_obj->n,
                                    ztgevc_obj->a,
                                    ztgevc_obj->lda,
                                    ztgevc_obj->b,
                                    ztgevc_obj->ldb,
                                    &ztgevc_obj->ilo,
                                    &ztgevc_obj->ihi,
                                    ztgevc_obj->lscale,
                                    ztgevc_obj->rscale
                                    );

    ztgevc_obj->inforef = ZGGHRD(   ztgevc_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    ztgevc_obj->n,
                                    ztgevc_obj->ilo,
                                    ztgevc_obj->ihi, 
                                    ztgevc_obj->a,
                                    ztgevc_obj->lda,
                                    ztgevc_obj->b,
                                    ztgevc_obj->ldb,
                                    ztgevc_obj->q,
                                    ztgevc_obj->ldq,
                                    ztgevc_obj->z,
                                    ztgevc_obj->ldz
                                    );
    ztgevc_obj->inforef = ZHGEQZ(   ztgevc_obj->matrix_layout,
                                    'S',
                                    'V',
									'V',
                                    ztgevc_obj->n,
                                    ztgevc_obj->ilo,
                                    ztgevc_obj->ihi, 
                                    ztgevc_obj->a,
                                    ztgevc_obj->lda,
                                    ztgevc_obj->b,
                                    ztgevc_obj->ldb,
									ztgevc_obj->alphar,
									ztgevc_obj->beta,
                                    ztgevc_obj->q,
                                    ztgevc_obj->ldq,
                                    ztgevc_obj->z,
                                    ztgevc_obj->ldz
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    ztgevc_obj->inforef = ZTGEVC(   ztgevc_obj->matrix_layout,
                                    ztgevc_obj->side,
                                    ztgevc_obj->howmny,
									ztgevc_obj->select,
                                    ztgevc_obj->n,
                                    ztgevc_obj->aref,
                                    ztgevc_obj->lda,
                                    ztgevc_obj->bref,
                                    ztgevc_obj->ldb,
                                    ztgevc_obj->qref,
                                    ztgevc_obj->ldq,
                                    ztgevc_obj->zref,
                                    ztgevc_obj->ldz,
									ztgevc_obj->mm,
									&ztgevc_obj->mref
                                    );

    /* Compute libflame's Lapacke o/p  */
    ztgevc_obj->inforef =  LAPACKE_ztgevc (ztgevc_obj->matrix_layout,
                                    ztgevc_obj->side,
                                    ztgevc_obj->howmny,
									ztgevc_obj->select,
                                    ztgevc_obj->n,
                                    ztgevc_obj->a,
                                    ztgevc_obj->lda,
                                    ztgevc_obj->b,
                                    ztgevc_obj->ldb,
                                    ztgevc_obj->q,
                                    ztgevc_obj->ldq,
                                    ztgevc_obj->z,
                                    ztgevc_obj->ldz,
									ztgevc_obj->mm,
									&ztgevc_obj->m
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    ztgevc_obj->diff_h =  computeDiff_z( (ztgevc_obj->lda)*(ztgevc_obj->n), 
                ztgevc_obj->a, ztgevc_obj->aref );

    ztgevc_obj->diff_t =  computeDiff_z( (ztgevc_obj->ldb)*(ztgevc_obj->n), 
                ztgevc_obj->b, ztgevc_obj->bref );

    ztgevc_obj->diff_q =  computeDiff_z( (ztgevc_obj->ldq)*(ztgevc_obj->n), 
                ztgevc_obj->q, ztgevc_obj->qref );

    ztgevc_obj->diff_z =  computeDiff_z( (ztgevc_obj->ldz)*(ztgevc_obj->n), 
                ztgevc_obj->z, ztgevc_obj->zref );

    ztgevc_obj->diff_alphar =  computeDiff_z( ztgevc_obj->n, 
                ztgevc_obj->alphar, ztgevc_obj->alpharref );

    ztgevc_obj->diff_alphai =  computeDiff_z( ztgevc_obj->n, 
                ztgevc_obj->alphai, ztgevc_obj->alphairef );

    ztgevc_obj->diff_beta =  computeDiff_z( ztgevc_obj->n, 
                ztgevc_obj->beta, ztgevc_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n tgevc lapack_complex_double: \n diff_h: %f \n diff_t: %f \n \
diff_q: %f \n diff_z: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       ztgevc_obj->diff_h, ztgevc_obj->diff_t, ztgevc_obj->diff_q,
	   ztgevc_obj->diff_z, ztgevc_obj->diff_alphar,
       ztgevc_obj->diff_alphai, ztgevc_obj->diff_beta );
#endif
}

TEST_F(ztgevc_test, ztgevc1) {
    EXPECT_NEAR(0.0, ztgevc_obj->diff_q, ztgevc_obj->threshold);
    EXPECT_NEAR(0.0, ztgevc_obj->diff_z, ztgevc_obj->threshold);
    EXPECT_EQ(ztgevc_obj->m, ztgevc_obj->mref);
}

TEST_F(ztgevc_test, ztgevc2) {
    EXPECT_NEAR(0.0, ztgevc_obj->diff_q, ztgevc_obj->threshold);
    EXPECT_NEAR(0.0, ztgevc_obj->diff_z, ztgevc_obj->threshold);
    EXPECT_EQ(ztgevc_obj->m, ztgevc_obj->mref);
}
TEST_F(ztgevc_test, ztgevc3) {
    EXPECT_NEAR(0.0, ztgevc_obj->diff_q, ztgevc_obj->threshold);
    EXPECT_NEAR(0.0, ztgevc_obj->diff_z, ztgevc_obj->threshold);
    EXPECT_EQ(ztgevc_obj->m, ztgevc_obj->mref);
}
TEST_F(ztgevc_test, ztgevc4) {
    EXPECT_NEAR(0.0, ztgevc_obj->diff_q, ztgevc_obj->threshold);
    EXPECT_NEAR(0.0, ztgevc_obj->diff_z, ztgevc_obj->threshold);
    EXPECT_EQ(ztgevc_obj->m, ztgevc_obj->mref);
}
