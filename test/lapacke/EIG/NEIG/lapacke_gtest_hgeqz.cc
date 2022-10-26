#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define LAPACKE_TEST_VERBOSE (1)

#define hgeqz_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (b1Q!=NULL)        free(b1Q); \
    if (b1Z!=NULL)     free(b1Z); \
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

/* Begin hgeqz_float_common_parameters  class definition */
class hgeqz_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_h, diff_t, diff_q, diff_z;
    float diff_alphai, diff_alphar, diff_beta;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'E' or 'S'.
	char compq; // Must be 'N', 'I', or 'V'.
	char compz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
	
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
    float* alphar, *alpharref;
	float* alphai, *alphairef;
	float* beta, *betaref;
    /*Return Values */
    int info, inforef;

    public:
      hgeqz_float_parameters (int matrix_layout_i, char job, char compq,
                         	  char compz, lapack_int n);
      ~hgeqz_float_parameters ();
};

/* Constructor definition  float_common_parameters */
hgeqz_float_parameters:: hgeqz_float_parameters (int matrix_layout_i, char job_i, 
                                      char compq_i, char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
	compq = compq_i;
	compz = compz_i;
    n = n_i;

    lda = n;
    ldb = n;
    ldq = n;
    ldz = n;

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
   printf(" \n hgeqz float: matrix_layout: %d n: %d  job: %c \n",
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
        (z==NULL) || (zref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "hgeqz_float_parameters object: malloc error. Exiting ";
       hgeqz_free();
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
    lapacke_gtest_init_float_buffer_pair_with_constant( lscale, lscaleref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( rscale, rscaleref, n, 0.0);
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( q, qref, n,ldq, 'I');
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( z, zref, n,ldq, 'I');
    lapacke_gtest_init_float_buffer_pair_with_constant( alphar, alpharref, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( alphai, alphairef, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( beta, betaref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'hgeqz_float_common_parameters' */
hgeqz_float_parameters :: ~hgeqz_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hgeqz_free();
} 

//  Test fixture class definition
class shgeqz_test  : public  ::testing::Test {
public:
   hgeqz_float_parameters  *shgeqz_obj;
   void SetUp();  
   void TearDown () { delete shgeqz_obj; }
};

void shgeqz_test::SetUp()
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

    /* LAPACKE_shgeqz prototype */	
    typedef int (*Fptr_NL_LAPACKE_shgeqz) ( int matrix_layout, char job,
		char compq, char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
		float* h, lapack_int ldh, float* t, lapack_int ldt, float* alphar,
		float* alphai, float* beta, float* q, lapack_int ldq, float* z,
		lapack_int ldz );
				 
    Fptr_NL_LAPACKE_shgeqz SHGEQZ;

    shgeqz_obj = new  hgeqz_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].job_seqr,
										 eig_non_sym_paramslist[idx].compq_hgeqz,
										 eig_non_sym_paramslist[idx].compz_hgeqz,
                                         eig_paramslist[idx].n );
                                         
    shgeqz_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    shgeqz_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    shgeqz_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(shgeqz_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(shgeqz_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGBAL = (Fptr_NL_LAPACKE_sggbal)dlsym(shgeqz_obj->hModule, "LAPACKE_sggbal");
    ASSERT_TRUE(SGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_sggbal symbol";

    SGGHRD = (Fptr_NL_LAPACKE_sgghrd)dlsym(shgeqz_obj->hModule, "LAPACKE_sgghrd");
    ASSERT_TRUE(SGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_sgghrd symbol";
 
    SHGEQZ = (Fptr_NL_LAPACKE_shgeqz)dlsym(shgeqz_obj->hModule, "LAPACKE_shgeqz");
    ASSERT_TRUE(SHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_shgeqz symbol";

	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    shgeqz_obj->info = LAPACKE_sgebrd(  shgeqz_obj->matrix_layout,
                                    shgeqz_obj->n,
                                    shgeqz_obj->n,
                                    shgeqz_obj->b1Z,
                                    shgeqz_obj->ldb,
                                    shgeqz_obj->d,
                                    shgeqz_obj->e,
                                    shgeqz_obj->tauq,
                                    shgeqz_obj->taup
                                    );

    shgeqz_obj->info = LAPACKE_sorgbr(  shgeqz_obj->matrix_layout,
                                    'Q',
                                    shgeqz_obj->n,
                                    shgeqz_obj->n,
                                    shgeqz_obj->n,
                                    shgeqz_obj->b1Z,
                                    shgeqz_obj->ldb,
                                    shgeqz_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (shgeqz_obj->z, shgeqz_obj->b1Z, sizeof(float)*(shgeqz_obj->ldb)*(shgeqz_obj->n) );
	memcpy (shgeqz_obj->zref, shgeqz_obj->b1Z, sizeof(float)*(shgeqz_obj->ldb)*(shgeqz_obj->n) );

    shgeqz_obj->info = LAPACKE_sgebrd(  shgeqz_obj->matrix_layout,
                                    shgeqz_obj->n,
                                    shgeqz_obj->n,
                                    shgeqz_obj->b1Q,
                                    shgeqz_obj->ldb,
                                    shgeqz_obj->d,
                                    shgeqz_obj->e,
                                    shgeqz_obj->tauq,
                                    shgeqz_obj->taup
                                    );

    shgeqz_obj->info = LAPACKE_sorgbr(  shgeqz_obj->matrix_layout,
                                    'Q',
                                    shgeqz_obj->n,
                                    shgeqz_obj->n,
                                    shgeqz_obj->n,
                                    shgeqz_obj->b1Q,
                                    shgeqz_obj->ldb,
                                    shgeqz_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (shgeqz_obj->q, shgeqz_obj->b1Q, sizeof(float)*(shgeqz_obj->ldb)*(shgeqz_obj->n) );
	memcpy (shgeqz_obj->qref, shgeqz_obj->b1Q, sizeof(float)*(shgeqz_obj->ldb)*(shgeqz_obj->n) );


    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    shgeqz_obj->inforef = SGGBAL(   shgeqz_obj->matrix_layout,
                                    'B',
                                    shgeqz_obj->n,
                                    shgeqz_obj->aref,
                                    shgeqz_obj->lda,
                                    shgeqz_obj->bref,
                                    shgeqz_obj->ldb,
                                    &shgeqz_obj->iloref,
                                    &shgeqz_obj->ihiref,
                                    shgeqz_obj->lscaleref,
                                    shgeqz_obj->rscaleref
                                    );


    shgeqz_obj->inforef = SGGHRD(   shgeqz_obj->matrix_layout,
                                    'I',
                                    'I', 
                                    shgeqz_obj->n,
                                    shgeqz_obj->iloref,
                                    shgeqz_obj->ihiref, 
                                    shgeqz_obj->aref,
                                    shgeqz_obj->lda,
                                    shgeqz_obj->bref,
                                    shgeqz_obj->ldb,
                                    shgeqz_obj->qref,
                                    shgeqz_obj->ldq,
                                    shgeqz_obj->zref,
                                    shgeqz_obj->ldz
                                    );

    shgeqz_obj->inforef = SHGEQZ(   shgeqz_obj->matrix_layout,
                                    shgeqz_obj->job,
                                    shgeqz_obj->compq,
									shgeqz_obj->compz,
                                    shgeqz_obj->n,
                                    shgeqz_obj->iloref,
                                    shgeqz_obj->ihiref, 
                                    shgeqz_obj->aref,
                                    shgeqz_obj->lda,
                                    shgeqz_obj->bref,
                                    shgeqz_obj->ldb,
									shgeqz_obj->alpharref,
									shgeqz_obj->alphairef,
									shgeqz_obj->betaref,
                                    shgeqz_obj->qref,
                                    shgeqz_obj->ldq,
                                    shgeqz_obj->zref,
                                    shgeqz_obj->ldz
                                    );

    /* Compute libflame's Lapacke o/p  */
    shgeqz_obj->inforef = SGGBAL(   shgeqz_obj->matrix_layout,
                                    'B',
                                    shgeqz_obj->n,
                                    shgeqz_obj->a,
                                    shgeqz_obj->lda,
                                    shgeqz_obj->b,
                                    shgeqz_obj->ldb,
                                    &shgeqz_obj->ilo,
                                    &shgeqz_obj->ihi,
                                    shgeqz_obj->lscale,
                                    shgeqz_obj->rscale
                                    );

    shgeqz_obj->inforef = SGGHRD(   shgeqz_obj->matrix_layout,
                                    'I',
                                    'I', 
                                    shgeqz_obj->n,
                                    shgeqz_obj->ilo,
                                    shgeqz_obj->ihi, 
                                    shgeqz_obj->a,
                                    shgeqz_obj->lda,
                                    shgeqz_obj->b,
                                    shgeqz_obj->ldb,
                                    shgeqz_obj->q,
                                    shgeqz_obj->ldq,
                                    shgeqz_obj->z,
                                    shgeqz_obj->ldz
                                    );

    shgeqz_obj->inforef = LAPACKE_shgeqz(   shgeqz_obj->matrix_layout,
                                    shgeqz_obj->job,
                                    shgeqz_obj->compq,
									shgeqz_obj->compz,
                                    shgeqz_obj->n,
                                    shgeqz_obj->ilo,
                                    shgeqz_obj->ihi, 
                                    shgeqz_obj->a,
                                    shgeqz_obj->lda,
                                    shgeqz_obj->b,
                                    shgeqz_obj->ldb,
									shgeqz_obj->alphar,
									shgeqz_obj->alphai,
									shgeqz_obj->beta,
                                    shgeqz_obj->q,
                                    shgeqz_obj->ldq,
                                    shgeqz_obj->z,
                                    shgeqz_obj->ldz
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    shgeqz_obj->diff_h =  computeDiff_s( (shgeqz_obj->lda)*(shgeqz_obj->n), 
                shgeqz_obj->a, shgeqz_obj->aref );

    shgeqz_obj->diff_t =  computeDiff_s( (shgeqz_obj->ldb)*(shgeqz_obj->n), 
                shgeqz_obj->b, shgeqz_obj->bref );

    shgeqz_obj->diff_q =  computeDiff_s( (shgeqz_obj->ldq)*(shgeqz_obj->n), 
                shgeqz_obj->q, shgeqz_obj->qref );

    shgeqz_obj->diff_z =  computeDiff_s( (shgeqz_obj->ldz)*(shgeqz_obj->n), 
                shgeqz_obj->z, shgeqz_obj->zref );

    shgeqz_obj->diff_alphar =  computeDiff_s( shgeqz_obj->n, 
                shgeqz_obj->alphar, shgeqz_obj->alpharref );

    shgeqz_obj->diff_alphai =  computeDiff_s( shgeqz_obj->n, 
                shgeqz_obj->alphai, shgeqz_obj->alphairef );

    shgeqz_obj->diff_beta =  computeDiff_s( shgeqz_obj->n, 
                shgeqz_obj->beta, shgeqz_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n hgeqz float: \n diff_h: %f \n diff_t: %f \n \
diff_q: %f \n diff_z: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       shgeqz_obj->diff_h, shgeqz_obj->diff_t, shgeqz_obj->diff_q,
	   shgeqz_obj->diff_z, shgeqz_obj->diff_alphar,
       shgeqz_obj->diff_alphai, shgeqz_obj->diff_beta );
#endif
}

TEST_F(shgeqz_test, shgeqz1) {
    EXPECT_NEAR(0.0, shgeqz_obj->diff_h, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_t, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_q, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_z, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_alphai, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_alphar, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_beta, shgeqz_obj->threshold);
}

TEST_F(shgeqz_test, shgeqz2) {
    EXPECT_NEAR(0.0, shgeqz_obj->diff_h, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_t, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_q, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_z, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_alphai, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_alphar, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_beta, shgeqz_obj->threshold);
}

TEST_F(shgeqz_test, shgeqz3) {
    EXPECT_NEAR(0.0, shgeqz_obj->diff_h, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_t, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_q, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_z, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_alphai, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_alphar, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_beta, shgeqz_obj->threshold);
}

TEST_F(shgeqz_test, shgeqz4) {
    EXPECT_NEAR(0.0, shgeqz_obj->diff_h, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_t, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_q, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_z, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_alphai, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_alphar, shgeqz_obj->threshold);
    EXPECT_NEAR(0.0, shgeqz_obj->diff_beta, shgeqz_obj->threshold);
}


/* Begin hgeqz_double_common_parameters  class definition */
class hgeqz_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_h, diff_t, diff_q, diff_z;
    double diff_alphai, diff_alphar, diff_beta;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'E' or 'S'.
	char compq; // Must be 'N', 'I', or 'V'.
	char compz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
	
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
    double* alphar, *alpharref;
	double* alphai, *alphairef;
	double* beta, *betaref;
    /*Return Values */
    int info, inforef;

    public:
      hgeqz_double_parameters (int matrix_layout_i, char job, char compq,
                         	  char compz, lapack_int n);
      ~hgeqz_double_parameters ();
};

/* Constructor definition  double_common_parameters */
hgeqz_double_parameters:: hgeqz_double_parameters (int matrix_layout_i, char job_i, 
                                      char compq_i, char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
	compq = compq_i;
	compz = compz_i;
    n = n_i;

    lda = n;
    ldb = n;
    ldq = n;
    ldz = n;

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
   printf(" \n hgeqz double: matrix_layout: %d n: %d  job: %c \n",
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
        (z==NULL) || (zref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "hgeqz_double_parameters object: malloc error. Exiting ";
       hgeqz_free();
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
    lapacke_gtest_init_double_buffer_pair_with_constant( lscale, lscaleref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( rscale, rscaleref, n, 0.0);
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( q, qref, n,ldq, 'I');
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( z, zref, n,ldq, 'I');
    lapacke_gtest_init_double_buffer_pair_with_constant( alphar, alpharref, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( alphai, alphairef, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( beta, betaref, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'hgeqz_double_common_parameters' */
hgeqz_double_parameters :: ~hgeqz_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hgeqz_free();
} 

//  Test fixture class definition
class dhgeqz_test  : public  ::testing::Test {
public:
   hgeqz_double_parameters  *dhgeqz_obj;
   void SetUp();  
   void TearDown () { delete dhgeqz_obj; }
};

void dhgeqz_test::SetUp()
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

    /* LAPACKE_dhgeqz prototype */	
    typedef int (*Fptr_NL_LAPACKE_dhgeqz) ( int matrix_layout, char job,
		char compq, char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
		double* h, lapack_int ldh, double* t, lapack_int ldt, double* alphar,
		double* alphai, double* beta, double* q, lapack_int ldq, double* z,
		lapack_int ldz );
				 
    Fptr_NL_LAPACKE_dhgeqz DHGEQZ;

    dhgeqz_obj = new  hgeqz_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].job_seqr,
										 eig_non_sym_paramslist[idx].compq_hgeqz,
										 eig_non_sym_paramslist[idx].compz_hgeqz,
                                         eig_paramslist[idx].n );
                                         
    dhgeqz_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
	
    dhgeqz_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dhgeqz_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dhgeqz_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dhgeqz_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGBAL = (Fptr_NL_LAPACKE_dggbal)dlsym(dhgeqz_obj->hModule, "LAPACKE_dggbal");
    ASSERT_TRUE(DGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_dggbal symbol";

    DGGHRD = (Fptr_NL_LAPACKE_dgghrd)dlsym(dhgeqz_obj->hModule, "LAPACKE_dgghrd");
    ASSERT_TRUE(DGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_dgghrd symbol";
 
    DHGEQZ = (Fptr_NL_LAPACKE_dhgeqz)dlsym(dhgeqz_obj->hModule, "LAPACKE_dhgeqz");
    ASSERT_TRUE(DHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_dhgeqz symbol";

	/* Prepare the i/p data to 'dhgeqz' API from the sequence of gebrd, dorgbr,
	   'dggbal', 'dgghrd' APIs for both NETLIB LAPACKE, AOCL-lapacke APIs    */
	   
	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    dhgeqz_obj->info = LAPACKE_dgebrd(  dhgeqz_obj->matrix_layout,
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->b1Z,
                                    dhgeqz_obj->ldb,
                                    dhgeqz_obj->d,
                                    dhgeqz_obj->e,
                                    dhgeqz_obj->tauq,
                                    dhgeqz_obj->taup
                                    );

    dhgeqz_obj->info = LAPACKE_dorgbr(  dhgeqz_obj->matrix_layout,
                                    'Q',
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->b1Z,
                                    dhgeqz_obj->ldb,
                                    dhgeqz_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (dhgeqz_obj->z, dhgeqz_obj->b1Z, sizeof(double)*(dhgeqz_obj->ldb)*(dhgeqz_obj->n) );
	memcpy (dhgeqz_obj->zref, dhgeqz_obj->b1Z, sizeof(double)*(dhgeqz_obj->ldb)*(dhgeqz_obj->n) );

	// Generate the orthogonal matrix 'Q' for i/p to gghrd call.
    dhgeqz_obj->info = LAPACKE_dgebrd(  dhgeqz_obj->matrix_layout,
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->b1Q,
                                    dhgeqz_obj->ldb,
                                    dhgeqz_obj->d,
                                    dhgeqz_obj->e,
                                    dhgeqz_obj->tauq,
                                    dhgeqz_obj->taup
                                    );

    dhgeqz_obj->info = LAPACKE_dorgbr(  dhgeqz_obj->matrix_layout,
                                    'Q',
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->b1Q,
                                    dhgeqz_obj->ldb,
                                    dhgeqz_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (dhgeqz_obj->q, dhgeqz_obj->b1Q, sizeof(double)*(dhgeqz_obj->ldb)*(dhgeqz_obj->n) );
	memcpy (dhgeqz_obj->qref, dhgeqz_obj->b1Q, sizeof(double)*(dhgeqz_obj->ldb)*(dhgeqz_obj->n) );


    dhgeqz_obj->inforef = DGGBAL(   dhgeqz_obj->matrix_layout,
                                    'B',
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->a,
                                    dhgeqz_obj->lda,
                                    dhgeqz_obj->b,
                                    dhgeqz_obj->ldb,
                                    &dhgeqz_obj->ilo,
                                    &dhgeqz_obj->ihi,
                                    dhgeqz_obj->lscale,
                                    dhgeqz_obj->rscale
                                    );

    dhgeqz_obj->inforef = DGGHRD(   dhgeqz_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->ilo,
                                    dhgeqz_obj->ihi, 
                                    dhgeqz_obj->a,
                                    dhgeqz_obj->lda,
                                    dhgeqz_obj->b,
                                    dhgeqz_obj->ldb,
                                    dhgeqz_obj->q,
                                    dhgeqz_obj->ldq,
                                    dhgeqz_obj->z,
                                    dhgeqz_obj->ldz
                                    );

    dhgeqz_obj->inforef = DGGBAL(   dhgeqz_obj->matrix_layout,
                                    'B',
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->aref,
                                    dhgeqz_obj->lda,
                                    dhgeqz_obj->bref,
                                    dhgeqz_obj->ldb,
                                    &dhgeqz_obj->iloref,
                                    &dhgeqz_obj->ihiref,
                                    dhgeqz_obj->lscaleref,
                                    dhgeqz_obj->rscaleref
                                    );

    dhgeqz_obj->inforef = DGGHRD(   dhgeqz_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->iloref,
                                    dhgeqz_obj->ihiref, 
                                    dhgeqz_obj->aref,
                                    dhgeqz_obj->lda,
                                    dhgeqz_obj->bref,
                                    dhgeqz_obj->ldb,
                                    dhgeqz_obj->qref,
                                    dhgeqz_obj->ldq,
                                    dhgeqz_obj->zref,
                                    dhgeqz_obj->ldz
                                    );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dhgeqz_obj->inforef = DHGEQZ(   dhgeqz_obj->matrix_layout,
                                    dhgeqz_obj->job,
                                    dhgeqz_obj->compq,
									dhgeqz_obj->compz,
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->iloref,
                                    dhgeqz_obj->ihiref, 
                                    dhgeqz_obj->aref,
                                    dhgeqz_obj->lda,
                                    dhgeqz_obj->bref,
                                    dhgeqz_obj->ldb,
									dhgeqz_obj->alpharref,
									dhgeqz_obj->alphairef,
									dhgeqz_obj->betaref,
                                    dhgeqz_obj->qref,
                                    dhgeqz_obj->ldq,
                                    dhgeqz_obj->zref,
                                    dhgeqz_obj->ldz
                                    );

    /* Compute libflame's Lapacke o/p  */
    dhgeqz_obj->inforef = LAPACKE_dhgeqz(   dhgeqz_obj->matrix_layout,
                                    dhgeqz_obj->job,
                                    dhgeqz_obj->compq,
									dhgeqz_obj->compz,
                                    dhgeqz_obj->n,
                                    dhgeqz_obj->ilo,
                                    dhgeqz_obj->ihi, 
                                    dhgeqz_obj->a,
                                    dhgeqz_obj->lda,
                                    dhgeqz_obj->b,
                                    dhgeqz_obj->ldb,
									dhgeqz_obj->alphar,
									dhgeqz_obj->alphai,
									dhgeqz_obj->beta,
                                    dhgeqz_obj->q,
                                    dhgeqz_obj->ldq,
                                    dhgeqz_obj->z,
                                    dhgeqz_obj->ldz
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    dhgeqz_obj->diff_h =  computeDiff_d( (dhgeqz_obj->lda)*(dhgeqz_obj->n), 
                dhgeqz_obj->a, dhgeqz_obj->aref );

    dhgeqz_obj->diff_t =  computeDiff_d( (dhgeqz_obj->ldb)*(dhgeqz_obj->n), 
                dhgeqz_obj->b, dhgeqz_obj->bref );

    dhgeqz_obj->diff_q =  computeDiff_d( (dhgeqz_obj->ldq)*(dhgeqz_obj->n), 
                dhgeqz_obj->q, dhgeqz_obj->qref );

    dhgeqz_obj->diff_z =  computeDiff_d( (dhgeqz_obj->ldz)*(dhgeqz_obj->n), 
                dhgeqz_obj->z, dhgeqz_obj->zref );

    dhgeqz_obj->diff_alphar =  computeDiff_d( dhgeqz_obj->n, 
                dhgeqz_obj->alphar, dhgeqz_obj->alpharref );

    dhgeqz_obj->diff_alphai =  computeDiff_d( dhgeqz_obj->n, 
                dhgeqz_obj->alphai, dhgeqz_obj->alphairef );

    dhgeqz_obj->diff_beta =  computeDiff_d( dhgeqz_obj->n, 
                dhgeqz_obj->beta, dhgeqz_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n hgeqz double: \n diff_h: %f \n diff_t: %f \n \
diff_q: %f \n diff_z: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       dhgeqz_obj->diff_h, dhgeqz_obj->diff_t, dhgeqz_obj->diff_q,
	   dhgeqz_obj->diff_z, dhgeqz_obj->diff_alphar,
       dhgeqz_obj->diff_alphai, dhgeqz_obj->diff_beta );
#endif
}

TEST_F(dhgeqz_test, dhgeqz1) {
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_h, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_t, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_q, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_z, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_alphai, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_alphar, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_beta, dhgeqz_obj->threshold);
}

TEST_F(dhgeqz_test, dhgeqz2) {
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_h, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_t, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_q, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_z, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_alphai, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_alphar, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_beta, dhgeqz_obj->threshold);
}

TEST_F(dhgeqz_test, dhgeqz3) {
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_h, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_t, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_q, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_z, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_alphai, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_alphar, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_beta, dhgeqz_obj->threshold);
}

TEST_F(dhgeqz_test, dhgeqz4) {
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_h, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_t, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_q, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_z, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_alphai, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_alphar, dhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, dhgeqz_obj->diff_beta, dhgeqz_obj->threshold);
}

/* Begin hgeqz_scomplex_common_parameters  class definition */
class hgeqz_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_h, diff_t, diff_q, diff_z;
    float diff_alphai, diff_alphar, diff_beta;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'E' or 'S'.
	char compq; // Must be 'N', 'I', or 'V'.
	char compz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
	
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
    lapack_complex_float* alphar, *alpharref;
	lapack_complex_float* alphai, *alphairef;
	lapack_complex_float* beta, *betaref;
    /*Return Values */
    int info, inforef;

    public:
      hgeqz_scomplex_parameters (int matrix_layout_i, char job, char compq,
                         	  char compz, lapack_int n);
      ~hgeqz_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
hgeqz_scomplex_parameters:: hgeqz_scomplex_parameters (int matrix_layout_i, char job_i, 
                                      char compq_i, char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
	compq = compq_i;
	compz = compz_i;
    n = n_i;

    lda = n;
    ldb = n;
    ldq = n;
    ldz = n;

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
   printf(" \n hgeqz lapack_complex_float: matrix_layout: %d n: %d  job: %c \n",
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
        (z==NULL) || (zref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "hgeqz_scomplex_parameters object: malloc error. Exiting ";
       hgeqz_free();
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
    
   } /* end of Constructor  */
    

/* Destructor definition  'hgeqz_scomplex_common_parameters' */
hgeqz_scomplex_parameters :: ~hgeqz_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hgeqz_free();
} 

//  Test fixture class definition
class chgeqz_test  : public  ::testing::Test {
public:
   hgeqz_scomplex_parameters  *chgeqz_obj;
   void SetUp();  
   void TearDown () { delete chgeqz_obj; }
};

void chgeqz_test::SetUp()
{
    /* LAPACKE_cggbal prototype */
    typedef int (*Fptr_NL_LAPACKE_cggbal) ( int matrix_layout, char job,
                    lapack_int n, lapack_complex_float* a, lapack_int lda, lapack_complex_float* b,
                    lapack_int ldb, lapack_int* ilo, lapack_int* ihi,
                    float* lscale, float* rscale);
                 
    Fptr_NL_LAPACKE_cggbal CGGBAL;
	
    /* LAPACKE_cgghrd prototype */	
    typedef int (*Fptr_NL_LAPACKE_cgghrd) ( int matrix_layout, char compq,
                 char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                 lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, lapack_int ldb, 
                 lapack_complex_float* q, lapack_int ldq, lapack_complex_float* z, lapack_int ldz);
				 
    Fptr_NL_LAPACKE_cgghrd CGGHRD;

    /* LAPACKE_chgeqz prototype */	
    typedef int (*Fptr_NL_LAPACKE_chgeqz) ( int matrix_layout, char job,
		char compq, char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
		lapack_complex_float* h, lapack_int ldh, lapack_complex_float* t, lapack_int ldt, lapack_complex_float* alpha,
		 lapack_complex_float* beta, lapack_complex_float* q, lapack_int ldq, lapack_complex_float* z,
		lapack_int ldz );
				 
    Fptr_NL_LAPACKE_chgeqz CHGEQZ;

    chgeqz_obj = new  hgeqz_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].job_seqr,
										 eig_non_sym_paramslist[idx].compq_hgeqz,
										 eig_non_sym_paramslist[idx].compz_hgeqz,
                                         eig_paramslist[idx].n );
                                         
    chgeqz_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    chgeqz_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chgeqz_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chgeqz_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chgeqz_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGBAL = (Fptr_NL_LAPACKE_cggbal)dlsym(chgeqz_obj->hModule, "LAPACKE_cggbal");
    ASSERT_TRUE(CGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_cggbal symbol";

    CGGHRD = (Fptr_NL_LAPACKE_cgghrd)dlsym(chgeqz_obj->hModule, "LAPACKE_cgghrd");
    ASSERT_TRUE(CGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_cgghrd symbol";
 
    CHGEQZ = (Fptr_NL_LAPACKE_chgeqz)dlsym(chgeqz_obj->hModule, "LAPACKE_chgeqz");
    ASSERT_TRUE(CHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_chgeqz symbol";

	/* Generate Q matrix for gghrd through gebrd & ungbr APIs */
    chgeqz_obj->info = LAPACKE_cgebrd(  chgeqz_obj->matrix_layout,
                                    chgeqz_obj->n,
                                    chgeqz_obj->n,
                                    chgeqz_obj->b1Z,
                                    chgeqz_obj->ldb,
                                    chgeqz_obj->d,
                                    chgeqz_obj->e,
                                    chgeqz_obj->tauq,
                                    chgeqz_obj->taup
                                    );

    chgeqz_obj->info = LAPACKE_cungbr(  chgeqz_obj->matrix_layout,
                                    'Q',
                                    chgeqz_obj->n,
                                    chgeqz_obj->n,
                                    chgeqz_obj->n,
                                    chgeqz_obj->b1Z,
                                    chgeqz_obj->ldb,
                                    chgeqz_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (chgeqz_obj->z, chgeqz_obj->b1Z, sizeof(float)*(chgeqz_obj->ldb)*(chgeqz_obj->n) );
	memcpy (chgeqz_obj->zref, chgeqz_obj->b1Z, sizeof(float)*(chgeqz_obj->ldb)*(chgeqz_obj->n) );

    chgeqz_obj->info = LAPACKE_cgebrd(  chgeqz_obj->matrix_layout,
                                    chgeqz_obj->n,
                                    chgeqz_obj->n,
                                    chgeqz_obj->b1Q,
                                    chgeqz_obj->ldb,
                                    chgeqz_obj->d,
                                    chgeqz_obj->e,
                                    chgeqz_obj->tauq,
                                    chgeqz_obj->taup
                                    );

    chgeqz_obj->info = LAPACKE_cungbr(  chgeqz_obj->matrix_layout,
                                    'Q',
                                    chgeqz_obj->n,
                                    chgeqz_obj->n,
                                    chgeqz_obj->n,
                                    chgeqz_obj->b1Q,
                                    chgeqz_obj->ldb,
                                    chgeqz_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (chgeqz_obj->q, chgeqz_obj->b1Q, sizeof(float)*(chgeqz_obj->ldb)*(chgeqz_obj->n) );
	memcpy (chgeqz_obj->qref, chgeqz_obj->b1Q, sizeof(float)*(chgeqz_obj->ldb)*(chgeqz_obj->n) );


    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    chgeqz_obj->inforef = CGGBAL(   chgeqz_obj->matrix_layout,
                                    'B',
                                    chgeqz_obj->n,
                                    chgeqz_obj->aref,
                                    chgeqz_obj->lda,
                                    chgeqz_obj->bref,
                                    chgeqz_obj->ldb,
                                    &chgeqz_obj->iloref,
                                    &chgeqz_obj->ihiref,
                                    chgeqz_obj->lscaleref,
                                    chgeqz_obj->rscaleref
                                    );

    chgeqz_obj->inforef = CGGHRD(   chgeqz_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    chgeqz_obj->n,
                                    chgeqz_obj->iloref,
                                    chgeqz_obj->ihiref, 
                                    chgeqz_obj->aref,
                                    chgeqz_obj->lda,
                                    chgeqz_obj->bref,
                                    chgeqz_obj->ldb,
                                    chgeqz_obj->qref,
                                    chgeqz_obj->ldq,
                                    chgeqz_obj->zref,
                                    chgeqz_obj->ldz
                                    );

    chgeqz_obj->inforef = CHGEQZ(   chgeqz_obj->matrix_layout,
                                    chgeqz_obj->job,
                                    chgeqz_obj->compq,
									chgeqz_obj->compz,
                                    chgeqz_obj->n,
                                    chgeqz_obj->iloref,
                                    chgeqz_obj->ihiref, 
                                    chgeqz_obj->aref,
                                    chgeqz_obj->lda,
                                    chgeqz_obj->bref,
                                    chgeqz_obj->ldb,
									chgeqz_obj->alpharref,
									chgeqz_obj->betaref,
                                    chgeqz_obj->qref,
                                    chgeqz_obj->ldq,
                                    chgeqz_obj->zref,
                                    chgeqz_obj->ldz
                                    );

    /* Compute libflame's Lapacke o/p  */
    chgeqz_obj->info = LAPACKE_cggbal(  chgeqz_obj->matrix_layout,
                                    'B',
                                    chgeqz_obj->n,
                                    chgeqz_obj->a,
                                    chgeqz_obj->lda,
                                    chgeqz_obj->b,
                                    chgeqz_obj->ldb,
                                    &chgeqz_obj->ilo,
                                    &chgeqz_obj->ihi,
                                    chgeqz_obj->lscale,
                                    chgeqz_obj->rscale
                                    );

    chgeqz_obj->info = LAPACKE_cgghrd(  chgeqz_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    chgeqz_obj->n,
                                    chgeqz_obj->ilo,
                                    chgeqz_obj->ihi, 
                                    chgeqz_obj->a,
                                    chgeqz_obj->lda,
                                    chgeqz_obj->b,
                                    chgeqz_obj->ldb,
                                    chgeqz_obj->q,
                                    chgeqz_obj->ldq,
                                    chgeqz_obj->z,
                                    chgeqz_obj->ldz
                                    );

    chgeqz_obj->inforef = LAPACKE_chgeqz(   chgeqz_obj->matrix_layout,
                                    chgeqz_obj->job,
                                    chgeqz_obj->compq,
									chgeqz_obj->compz,
                                    chgeqz_obj->n,
                                    chgeqz_obj->ilo,
                                    chgeqz_obj->ihi, 
                                    chgeqz_obj->a,
                                    chgeqz_obj->lda,
                                    chgeqz_obj->b,
                                    chgeqz_obj->ldb,
									chgeqz_obj->alphar,
									chgeqz_obj->beta,
                                    chgeqz_obj->q,
                                    chgeqz_obj->ldq,
                                    chgeqz_obj->z,
                                    chgeqz_obj->ldz
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    chgeqz_obj->diff_h =  computeDiff_c( (chgeqz_obj->lda)*(chgeqz_obj->n), 
                chgeqz_obj->a, chgeqz_obj->aref );

    chgeqz_obj->diff_t =  computeDiff_c( (chgeqz_obj->ldb)*(chgeqz_obj->n), 
                chgeqz_obj->b, chgeqz_obj->bref );

    chgeqz_obj->diff_q =  computeDiff_c( (chgeqz_obj->ldq)*(chgeqz_obj->n), 
                chgeqz_obj->q, chgeqz_obj->qref );

    chgeqz_obj->diff_z =  computeDiff_c( (chgeqz_obj->ldz)*(chgeqz_obj->n), 
                chgeqz_obj->z, chgeqz_obj->zref );

    chgeqz_obj->diff_alphar =  computeDiff_c( chgeqz_obj->n, 
                chgeqz_obj->alphar, chgeqz_obj->alpharref );

    chgeqz_obj->diff_alphai =  computeDiff_c( chgeqz_obj->n, 
                chgeqz_obj->alphai, chgeqz_obj->alphairef );

    chgeqz_obj->diff_beta =  computeDiff_c( chgeqz_obj->n, 
                chgeqz_obj->beta, chgeqz_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n hgeqz lapack_complex_float: \n diff_h: %f \n diff_t: %f \n \
diff_q: %f \n diff_z: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       chgeqz_obj->diff_h, chgeqz_obj->diff_t, chgeqz_obj->diff_q,
	   chgeqz_obj->diff_z, chgeqz_obj->diff_alphar,
       chgeqz_obj->diff_alphai, chgeqz_obj->diff_beta );
#endif
}

TEST_F(chgeqz_test, chgeqz1) {
    EXPECT_NEAR(0.0, chgeqz_obj->diff_h, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_t, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_q, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_z, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_alphai, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_alphar, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_beta, chgeqz_obj->threshold);
}

TEST_F(chgeqz_test, chgeqz2) {
    EXPECT_NEAR(0.0, chgeqz_obj->diff_h, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_t, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_q, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_z, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_alphai, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_alphar, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_beta, chgeqz_obj->threshold);
}

TEST_F(chgeqz_test, chgeqz3) {
    EXPECT_NEAR(0.0, chgeqz_obj->diff_h, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_t, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_q, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_z, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_alphai, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_alphar, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_beta, chgeqz_obj->threshold);
}

TEST_F(chgeqz_test, chgeqz4) {
    EXPECT_NEAR(0.0, chgeqz_obj->diff_h, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_t, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_q, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_z, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_alphai, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_alphar, chgeqz_obj->threshold);
    EXPECT_NEAR(0.0, chgeqz_obj->diff_beta, chgeqz_obj->threshold);
}

/* Begin hgeqz_dcomplex_common_parameters  class definition */
class hgeqz_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_h, diff_t, diff_q, diff_z;
    double diff_alphai, diff_alphar, diff_beta;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'E' or 'S'.
	char compq; // Must be 'N', 'I', or 'V'.
	char compz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
	
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
    lapack_complex_double* alphar, *alpharref;
	lapack_complex_double* alphai, *alphairef;
	lapack_complex_double* beta, *betaref;
    /*Return Values */
    int info, inforef;

    public:
      hgeqz_dcomplex_parameters (int matrix_layout_i, char job, char compq,
                         	  char compz, lapack_int n);
      ~hgeqz_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
hgeqz_dcomplex_parameters:: hgeqz_dcomplex_parameters (int matrix_layout_i, char job_i, 
                                      char compq_i, char compz_i, lapack_int n_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
	compq = compq_i;
	compz = compz_i;
    n = n_i;

    lda = n;
    ldb = n;
    ldq = n;
    ldz = n;

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
   printf(" \n hgeqz lapack_complex_double: matrix_layout: %d n: %d  job: %c \n",
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
        (z==NULL) || (zref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "hgeqz_dcomplex_parameters object: malloc error. Exiting ";
       hgeqz_free();
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
    
   } /* end of Constructor  */
    

/* Destructor definition  'hgeqz_dcomplex_common_parameters' */
hgeqz_dcomplex_parameters :: ~hgeqz_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   hgeqz_free();
} 

//  Test fixture class definition
class zhgeqz_test  : public  ::testing::Test {
public:
   hgeqz_dcomplex_parameters  *zhgeqz_obj;
   void SetUp();  
   void TearDown () { delete zhgeqz_obj; }
};

void zhgeqz_test::SetUp()
{
    /* LAPACKE_zggbal prototype */
    typedef int (*Fptr_NL_LAPACKE_zggbal) ( int matrix_layout, char job,
                    lapack_int n, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b,
                    lapack_int ldb, lapack_int* ilo, lapack_int* ihi,
                    double* lscale, double* rscale);
                 
    Fptr_NL_LAPACKE_zggbal ZGGBAL;
	
    /* LAPACKE_zgghrd prototype */	
    typedef int (*Fptr_NL_LAPACKE_zgghrd) ( int matrix_layout, char compq,
                 char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
                 lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, 
                 lapack_complex_double* q, lapack_int ldq, lapack_complex_double* z, lapack_int ldz);
				 
    Fptr_NL_LAPACKE_zgghrd ZGGHRD;

    /* LAPACKE_zhgeqz prototype */	
    typedef int (*Fptr_NL_LAPACKE_zhgeqz) ( int matrix_layout, char job,
		char compq, char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
		lapack_complex_double* h, lapack_int ldh, lapack_complex_double* t, lapack_int ldt, lapack_complex_double* alpha,
		 lapack_complex_double* beta, lapack_complex_double* q, lapack_int ldq, lapack_complex_double* z,
		lapack_int ldz );
				 
    Fptr_NL_LAPACKE_zhgeqz ZHGEQZ;

    zhgeqz_obj = new  hgeqz_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].job_seqr,
										 eig_non_sym_paramslist[idx].compq_hgeqz,
										 eig_non_sym_paramslist[idx].compz_hgeqz,
                                         eig_paramslist[idx].n );
                                         
    zhgeqz_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    zhgeqz_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhgeqz_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhgeqz_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhgeqz_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGBAL = (Fptr_NL_LAPACKE_zggbal)dlsym(zhgeqz_obj->hModule, "LAPACKE_zggbal");
    ASSERT_TRUE(ZGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_zggbal symbol";

    ZGGHRD = (Fptr_NL_LAPACKE_zgghrd)dlsym(zhgeqz_obj->hModule, "LAPACKE_zgghrd");
    ASSERT_TRUE(ZGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_zgghrd symbol";
 
    ZHGEQZ = (Fptr_NL_LAPACKE_zhgeqz)dlsym(zhgeqz_obj->hModule, "LAPACKE_zhgeqz");
    ASSERT_TRUE(ZHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_zhgeqz symbol";

	/* Generate Q matrix for gghrd through gebrd & ungbr APIs */
    zhgeqz_obj->info = LAPACKE_zgebrd(  zhgeqz_obj->matrix_layout,
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->b1Z,
                                    zhgeqz_obj->ldb,
                                    zhgeqz_obj->d,
                                    zhgeqz_obj->e,
                                    zhgeqz_obj->tauq,
                                    zhgeqz_obj->taup
                                    );

    zhgeqz_obj->info = LAPACKE_zungbr(  zhgeqz_obj->matrix_layout,
                                    'Q',
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->b1Z,
                                    zhgeqz_obj->ldb,
                                    zhgeqz_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (zhgeqz_obj->z, zhgeqz_obj->b1Z, sizeof(double)*(zhgeqz_obj->ldb)*(zhgeqz_obj->n) );
	memcpy (zhgeqz_obj->zref, zhgeqz_obj->b1Z, sizeof(double)*(zhgeqz_obj->ldb)*(zhgeqz_obj->n) );

    zhgeqz_obj->info = LAPACKE_zgebrd(  zhgeqz_obj->matrix_layout,
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->b1Q,
                                    zhgeqz_obj->ldb,
                                    zhgeqz_obj->d,
                                    zhgeqz_obj->e,
                                    zhgeqz_obj->tauq,
                                    zhgeqz_obj->taup
                                    );

    zhgeqz_obj->info = LAPACKE_zungbr(  zhgeqz_obj->matrix_layout,
                                    'Q',
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->b1Q,
                                    zhgeqz_obj->ldb,
                                    zhgeqz_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (zhgeqz_obj->q, zhgeqz_obj->b1Q, sizeof(double)*(zhgeqz_obj->ldb)*(zhgeqz_obj->n) );
	memcpy (zhgeqz_obj->qref, zhgeqz_obj->b1Q, sizeof(double)*(zhgeqz_obj->ldb)*(zhgeqz_obj->n) );


    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zhgeqz_obj->inforef = ZGGBAL(   zhgeqz_obj->matrix_layout,
                                    'B',
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->aref,
                                    zhgeqz_obj->lda,
                                    zhgeqz_obj->bref,
                                    zhgeqz_obj->ldb,
                                    &zhgeqz_obj->iloref,
                                    &zhgeqz_obj->ihiref,
                                    zhgeqz_obj->lscaleref,
                                    zhgeqz_obj->rscaleref
                                    );

    zhgeqz_obj->inforef = ZGGHRD(   zhgeqz_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->iloref,
                                    zhgeqz_obj->ihiref, 
                                    zhgeqz_obj->aref,
                                    zhgeqz_obj->lda,
                                    zhgeqz_obj->bref,
                                    zhgeqz_obj->ldb,
                                    zhgeqz_obj->qref,
                                    zhgeqz_obj->ldq,
                                    zhgeqz_obj->zref,
                                    zhgeqz_obj->ldz
                                    );

    zhgeqz_obj->inforef = ZHGEQZ(   zhgeqz_obj->matrix_layout,
                                    zhgeqz_obj->job,
                                    zhgeqz_obj->compq,
									zhgeqz_obj->compz,
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->iloref,
                                    zhgeqz_obj->ihiref, 
                                    zhgeqz_obj->aref,
                                    zhgeqz_obj->lda,
                                    zhgeqz_obj->bref,
                                    zhgeqz_obj->ldb,
									zhgeqz_obj->alpharref,
									zhgeqz_obj->betaref,
                                    zhgeqz_obj->qref,
                                    zhgeqz_obj->ldq,
                                    zhgeqz_obj->zref,
                                    zhgeqz_obj->ldz
                                    );

    /* Compute libflame's Lapacke o/p  */
    zhgeqz_obj->info = LAPACKE_zggbal(  zhgeqz_obj->matrix_layout,
                                    'B',
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->a,
                                    zhgeqz_obj->lda,
                                    zhgeqz_obj->b,
                                    zhgeqz_obj->ldb,
                                    &zhgeqz_obj->ilo,
                                    &zhgeqz_obj->ihi,
                                    zhgeqz_obj->lscale,
                                    zhgeqz_obj->rscale
                                    );

    zhgeqz_obj->info = LAPACKE_zgghrd(  zhgeqz_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->ilo,
                                    zhgeqz_obj->ihi, 
                                    zhgeqz_obj->a,
                                    zhgeqz_obj->lda,
                                    zhgeqz_obj->b,
                                    zhgeqz_obj->ldb,
                                    zhgeqz_obj->q,
                                    zhgeqz_obj->ldq,
                                    zhgeqz_obj->z,
                                    zhgeqz_obj->ldz
                                    );

    zhgeqz_obj->inforef = LAPACKE_zhgeqz(   zhgeqz_obj->matrix_layout,
                                    zhgeqz_obj->job,
                                    zhgeqz_obj->compq,
									zhgeqz_obj->compz,
                                    zhgeqz_obj->n,
                                    zhgeqz_obj->ilo,
                                    zhgeqz_obj->ihi, 
                                    zhgeqz_obj->a,
                                    zhgeqz_obj->lda,
                                    zhgeqz_obj->b,
                                    zhgeqz_obj->ldb,
									zhgeqz_obj->alphar,
									zhgeqz_obj->beta,
                                    zhgeqz_obj->q,
                                    zhgeqz_obj->ldq,
                                    zhgeqz_obj->z,
                                    zhgeqz_obj->ldz
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zhgeqz_obj->diff_h =  computeDiff_z( (zhgeqz_obj->lda)*(zhgeqz_obj->n), 
                zhgeqz_obj->a, zhgeqz_obj->aref );

    zhgeqz_obj->diff_t =  computeDiff_z( (zhgeqz_obj->ldb)*(zhgeqz_obj->n), 
                zhgeqz_obj->b, zhgeqz_obj->bref );

    zhgeqz_obj->diff_q =  computeDiff_z( (zhgeqz_obj->ldq)*(zhgeqz_obj->n), 
                zhgeqz_obj->q, zhgeqz_obj->qref );

    zhgeqz_obj->diff_z =  computeDiff_z( (zhgeqz_obj->ldz)*(zhgeqz_obj->n), 
                zhgeqz_obj->z, zhgeqz_obj->zref );

    zhgeqz_obj->diff_alphar =  computeDiff_z( zhgeqz_obj->n, 
                zhgeqz_obj->alphar, zhgeqz_obj->alpharref );

    zhgeqz_obj->diff_alphai =  computeDiff_z( zhgeqz_obj->n, 
                zhgeqz_obj->alphai, zhgeqz_obj->alphairef );

    zhgeqz_obj->diff_beta =  computeDiff_z( zhgeqz_obj->n, 
                zhgeqz_obj->beta, zhgeqz_obj->betaref );

#if LAPACKE_TEST_VERBOSE
   printf(" \n hgeqz lapack_complex_double: \n diff_h: %f \n diff_t: %f \n \
diff_q: %f \n diff_z: %f \n diff_alphar: %f \n\
diff_alphai: %f \n diff_beta: %f \n ",
       zhgeqz_obj->diff_h, zhgeqz_obj->diff_t, zhgeqz_obj->diff_q,
	   zhgeqz_obj->diff_z, zhgeqz_obj->diff_alphar,
       zhgeqz_obj->diff_alphai, zhgeqz_obj->diff_beta );
#endif
}

TEST_F(zhgeqz_test, zhgeqz1) {
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_h, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_t, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_q, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_z, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_alphai, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_alphar, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_beta, zhgeqz_obj->threshold);
}

TEST_F(zhgeqz_test, zhgeqz2) {
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_h, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_t, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_q, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_z, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_alphai, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_alphar, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_beta, zhgeqz_obj->threshold);
}

TEST_F(zhgeqz_test, zhgeqz3) {
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_h, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_t, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_q, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_z, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_alphai, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_alphar, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_beta, zhgeqz_obj->threshold);
}

TEST_F(zhgeqz_test, zhgeqz4) {
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_h, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_t, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_q, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_z, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_alphai, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_alphar, zhgeqz_obj->threshold);
    EXPECT_NEAR(0.0, zhgeqz_obj->diff_beta, zhgeqz_obj->threshold);
}
