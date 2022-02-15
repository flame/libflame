#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define tgexc_free() \
    if (a!=NULL)        free(a); \
    if (aref!=NULL)     free(aref); \
    if (b!=NULL)        free(b); \
    if (bref!=NULL)     free(bref); \
    if (b1Q!=NULL)        free(b1Q); \
    if (b1Z!=NULL)     free(b1Z); \
    if (z!=NULL)        free(z); \
    if (zref!=NULL)     free(zref); \
    if (d!=NULL)     free(d); \
    if (e!=NULL)     free(e); \
    if (q!=NULL)        free(q); \
    if (qref!=NULL)     free(qref)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin tgexc_float_common_parameters  class definition */
class tgexc_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_q, diff_z;
	float threshold;
    void *hModule, *dModule;
    
	// Intermediate buffers to generate q, z i/ps for tgexc API
    float* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    float* d, *e;
    float* tauq, *taup;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    
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
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    float* z, *zref; //
    lapack_int ldz; // The leading dimension of z;
    
    /*Return Values */
    int info, inforef;

   public:
      tgexc_float_parameters (int matrix_layout_i, char wantq, char wantz,
                              lapack_int n);

      ~tgexc_float_parameters ();
};

/* Constructor definition  float_common_parameters */
tgexc_float_parameters:: tgexc_float_parameters (int matrix_layout_i,
                            char wantq_i, char wantz_i, lapack_int n_i)
{
	int randIndex1, randIndex2;

    matrix_layout = matrix_layout_i;
    wantq = wantq_i;
    wantz = wantz_i;
    n = n_i;

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
   printf(" \n tgexc float: matrix_layout: %d n: %d  wantz: %c \
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

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (b1Q==NULL) || (b1Z==NULL) || \
        (z==NULL) || (zref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "tgexc_float_parameters object: malloc error. Exiting ";
       tgexc_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_float_buffer_pair_rand( b1Q, b1Z, ldb*n );
	memcpy (b1Q, b, sizeof(float)*ldb*n );
    lapacke_gtest_init_float_buffer_pair_rand( q, qref, ldq*n );
    lapacke_gtest_init_float_buffer_pair_rand( z, zref, ldz*n );
	lapacke_gtest_init_float_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_float_buffer_pair_with_constant( taup, tauq, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'tgexc_float_common_parameters' */
tgexc_float_parameters :: ~tgexc_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgexc_free();
} 

//  Test fixture class definition
class stgexc_test  : public  ::testing::Test {
public:
   tgexc_float_parameters  *stgexc_obj;
   void SetUp();  
   void TearDown () { delete stgexc_obj; }
};

void stgexc_test::SetUp()
{
    /* LAPACKE STGEXC prototype */
    typedef int (*Fptr_NL_LAPACKE_stgexc) ( int matrix_layout,
		lapack_logical wantq, lapack_logical wantz, lapack_int n,
		float* a, lapack_int lda, float* b, lapack_int ldb, float* q,
		lapack_int ldq, float* z, lapack_int ldz, lapack_int* ifst,
		lapack_int* ilst);

    Fptr_NL_LAPACKE_stgexc STGEXC;

    stgexc_obj = new  tgexc_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].wantz,
                                         eig_non_sym_paramslist[idx].wantq,
                                         eig_paramslist[idx].n );
										 
    stgexc_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    stgexc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stgexc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stgexc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stgexc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    STGEXC = (Fptr_NL_LAPACKE_stgexc)dlsym(stgexc_obj->hModule, "LAPACKE_stgexc");
    ASSERT_TRUE(STGEXC != NULL) << "failed to ppt the Netlib LAPACKE_stgexc symbol";
 
 	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    stgexc_obj->info = LAPACKE_sgebrd(  stgexc_obj->matrix_layout,
                                    stgexc_obj->n,
                                    stgexc_obj->n,
                                    stgexc_obj->b1Z,
                                    stgexc_obj->ldb,
                                    stgexc_obj->d,
                                    stgexc_obj->e,
                                    stgexc_obj->tauq,
                                    stgexc_obj->taup
                                    );

    stgexc_obj->info = LAPACKE_sorgbr(  stgexc_obj->matrix_layout,
                                    'Q',
                                    stgexc_obj->n,
                                    stgexc_obj->n,
                                    stgexc_obj->n,
                                    stgexc_obj->b1Z,
                                    stgexc_obj->ldb,
                                    stgexc_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (stgexc_obj->z, stgexc_obj->b1Z, sizeof(float)*(stgexc_obj->ldb)*(stgexc_obj->n) );
	memcpy (stgexc_obj->zref, stgexc_obj->b1Z, sizeof(float)*(stgexc_obj->ldb)*(stgexc_obj->n) );

    stgexc_obj->info = LAPACKE_sgebrd(  stgexc_obj->matrix_layout,
                                    stgexc_obj->n,
                                    stgexc_obj->n,
                                    stgexc_obj->b1Q,
                                    stgexc_obj->ldb,
                                    stgexc_obj->d,
                                    stgexc_obj->e,
                                    stgexc_obj->tauq,
                                    stgexc_obj->taup
                                    );

    stgexc_obj->info = LAPACKE_sorgbr(  stgexc_obj->matrix_layout,
                                    'Q',
                                    stgexc_obj->n,
                                    stgexc_obj->n,
                                    stgexc_obj->n,
                                    stgexc_obj->b1Q,
                                    stgexc_obj->ldb,
                                    stgexc_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (stgexc_obj->q, stgexc_obj->b1Q, sizeof(float)*(stgexc_obj->ldb)*(stgexc_obj->n) );
	memcpy (stgexc_obj->qref, stgexc_obj->b1Q, sizeof(float)*(stgexc_obj->ldb)*(stgexc_obj->n) );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    stgexc_obj->inforef = STGEXC(   stgexc_obj->matrix_layout,
                                    (lapack_logical)stgexc_obj->wantq,
                                    (lapack_logical)stgexc_obj->wantz, 
                                    stgexc_obj->n,
                                    stgexc_obj->aref,
                                    stgexc_obj->lda,
                                    stgexc_obj->bref,
                                    stgexc_obj->ldb,
                                    stgexc_obj->qref,
                                    stgexc_obj->ldq,
                                    stgexc_obj->zref,
                                    stgexc_obj->ldz,
                                    &stgexc_obj->ifstref,
                                    &stgexc_obj->ilstref
                                    );

    /* Compute libflame's Lapacke o/p  */
    stgexc_obj->info = LAPACKE_stgexc(  stgexc_obj->matrix_layout,
                                    (lapack_logical)stgexc_obj->wantq,
                                    (lapack_logical)stgexc_obj->wantz, 
                                    stgexc_obj->n,
                                    stgexc_obj->a,
                                    stgexc_obj->lda,
                                    stgexc_obj->b,
                                    stgexc_obj->ldb,
                                    stgexc_obj->q,
                                    stgexc_obj->ldq,
                                    stgexc_obj->z,
                                    stgexc_obj->ldz,
                                    &stgexc_obj->ifst,
                                    &stgexc_obj->ilst
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    stgexc_obj->diff_a =  computeDiff_s( (stgexc_obj->lda)*(stgexc_obj->n), 
                stgexc_obj->a, stgexc_obj->aref );

    stgexc_obj->diff_b =  computeDiff_s( (stgexc_obj->ldb)*(stgexc_obj->n), 
                stgexc_obj->b, stgexc_obj->bref );

    stgexc_obj->diff_q =  computeDiff_s( (stgexc_obj->ldq)*(stgexc_obj->n), 
                stgexc_obj->q, stgexc_obj->qref );

    stgexc_obj->diff_z =  computeDiff_s( (stgexc_obj->ldz)*(stgexc_obj->n), 
                stgexc_obj->z, stgexc_obj->zref );
}

TEST_F(stgexc_test, stgexc1) {
    EXPECT_NEAR(0.0, stgexc_obj->diff_a, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_b, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_q, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_z, stgexc_obj->threshold);
}

TEST_F(stgexc_test, stgexc2) {
    EXPECT_NEAR(0.0, stgexc_obj->diff_a, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_b, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_q, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_z, stgexc_obj->threshold);
}

TEST_F(stgexc_test, stgexc3) {
    EXPECT_NEAR(0.0, stgexc_obj->diff_a, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_b, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_q, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_z, stgexc_obj->threshold);
}

TEST_F(stgexc_test, stgexc4) {
    EXPECT_NEAR(0.0, stgexc_obj->diff_a, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_b, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_q, stgexc_obj->threshold);
    EXPECT_NEAR(0.0, stgexc_obj->diff_z, stgexc_obj->threshold);
}

/* Begin tgexc_double_common_parameters  class definition */
class tgexc_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b, diff_q, diff_z;
	float threshold;
    void *hModule, *dModule;
    
	// Intermediate buffers to generate q, z i/ps for tgexc API
    double* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    double* d, *e;
    double* tauq, *taup;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    
    char wantq; // Must be 'N', 'I', or 'V'.
    char wantz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int ifst, ifstref; // ifst and ilst mark the rows and columns of A which are to be reduced.
    lapack_int ilst, ilstref;
    
    /* Input / Output parameters */
    
    double* a, *aref; // contains the n-by-n general matrix A.
    lapack_int lda; //  The leading dimension of a
    double* b, *bref; // contains the n-by-n upper triangular matrix B.
    lapack_int ldb; //  The leading dimension of b
    double* q, *qref; // If wantq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    double* z, *zref; //
    lapack_int ldz; // The leading dimension of z;
    
    /*Return Values */
    int info, inforef;

   public:
      tgexc_double_parameters (int matrix_layout_i, char wantq, char wantz,
                              lapack_int n);

      ~tgexc_double_parameters ();
};

/* Constructor definition  double_common_parameters */
tgexc_double_parameters:: tgexc_double_parameters (int matrix_layout_i,
                            char wantq_i, char wantz_i, lapack_int n_i)
{
	int randIndex1, randIndex2;

    matrix_layout = matrix_layout_i;
    wantq = wantq_i;
    wantz = wantz_i;
    n = n_i;

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
   printf(" \n tgexc double: matrix_layout: %d n: %d  wantz: %c \
       wantq: %c  \n", matrix_layout, n, wantz, wantq);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_double_buffer_pair( &z, &zref, ldz*n );
    lapacke_gtest_alloc_double_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_double_buffer_pair( &taup, &tauq, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (b1Q==NULL) || (b1Z==NULL) || \
        (z==NULL) || (zref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "tgexc_double_parameters object: malloc error. Exiting ";
       tgexc_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_double_buffer_pair_rand( b1Q, b1Z, ldb*n );
	memcpy (b1Q, b, sizeof(double)*ldb*n );
    lapacke_gtest_init_double_buffer_pair_rand( q, qref, ldq*n );
    lapacke_gtest_init_double_buffer_pair_rand( z, zref, ldz*n );
	lapacke_gtest_init_double_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_double_buffer_pair_with_constant( taup, tauq, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'tgexc_double_common_parameters' */
tgexc_double_parameters :: ~tgexc_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgexc_free();
} 

//  Test fixture class definition
class dtgexc_test  : public  ::testing::Test {
public:
   tgexc_double_parameters  *dtgexc_obj;
   void SetUp();  
   void TearDown () { delete dtgexc_obj; }
};

void dtgexc_test::SetUp()
{
    /* LAPACKE DTGEXC prototype */
    typedef int (*Fptr_NL_LAPACKE_dtgexc) ( int matrix_layout,
		lapack_logical wantq, lapack_logical wantz, lapack_int n,
		double* a, lapack_int lda, double* b, lapack_int ldb, double* q,
		lapack_int ldq, double* z, lapack_int ldz, lapack_int* ifst,
		lapack_int* ilst);

    Fptr_NL_LAPACKE_dtgexc DTGEXC;

    dtgexc_obj = new  tgexc_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].wantz,
                                         eig_non_sym_paramslist[idx].wantq,
                                         eig_paramslist[idx].n );
										 
    dtgexc_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    dtgexc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtgexc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtgexc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtgexc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DTGEXC = (Fptr_NL_LAPACKE_dtgexc)dlsym(dtgexc_obj->hModule, "LAPACKE_dtgexc");
    ASSERT_TRUE(DTGEXC != NULL) << "failed to ppt the Netlib LAPACKE_dtgexc symbol";
 
 	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    dtgexc_obj->info = LAPACKE_dgebrd(  dtgexc_obj->matrix_layout,
                                    dtgexc_obj->n,
                                    dtgexc_obj->n,
                                    dtgexc_obj->b1Z,
                                    dtgexc_obj->ldb,
                                    dtgexc_obj->d,
                                    dtgexc_obj->e,
                                    dtgexc_obj->tauq,
                                    dtgexc_obj->taup
                                    );

    dtgexc_obj->info = LAPACKE_dorgbr(  dtgexc_obj->matrix_layout,
                                    'Q',
                                    dtgexc_obj->n,
                                    dtgexc_obj->n,
                                    dtgexc_obj->n,
                                    dtgexc_obj->b1Z,
                                    dtgexc_obj->ldb,
                                    dtgexc_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (dtgexc_obj->z, dtgexc_obj->b1Z, sizeof(double)*(dtgexc_obj->ldb)*(dtgexc_obj->n) );
	memcpy (dtgexc_obj->zref, dtgexc_obj->b1Z, sizeof(double)*(dtgexc_obj->ldb)*(dtgexc_obj->n) );

    dtgexc_obj->info = LAPACKE_dgebrd(  dtgexc_obj->matrix_layout,
                                    dtgexc_obj->n,
                                    dtgexc_obj->n,
                                    dtgexc_obj->b1Q,
                                    dtgexc_obj->ldb,
                                    dtgexc_obj->d,
                                    dtgexc_obj->e,
                                    dtgexc_obj->tauq,
                                    dtgexc_obj->taup
                                    );

    dtgexc_obj->info = LAPACKE_dorgbr(  dtgexc_obj->matrix_layout,
                                    'Q',
                                    dtgexc_obj->n,
                                    dtgexc_obj->n,
                                    dtgexc_obj->n,
                                    dtgexc_obj->b1Q,
                                    dtgexc_obj->ldb,
                                    dtgexc_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (dtgexc_obj->q, dtgexc_obj->b1Q, sizeof(double)*(dtgexc_obj->ldb)*(dtgexc_obj->n) );
	memcpy (dtgexc_obj->qref, dtgexc_obj->b1Q, sizeof(double)*(dtgexc_obj->ldb)*(dtgexc_obj->n) );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dtgexc_obj->inforef = DTGEXC(   dtgexc_obj->matrix_layout,
                                    (lapack_logical)dtgexc_obj->wantq,
                                    (lapack_logical)dtgexc_obj->wantz, 
                                    dtgexc_obj->n,
                                    dtgexc_obj->aref,
                                    dtgexc_obj->lda,
                                    dtgexc_obj->bref,
                                    dtgexc_obj->ldb,
                                    dtgexc_obj->qref,
                                    dtgexc_obj->ldq,
                                    dtgexc_obj->zref,
                                    dtgexc_obj->ldz,
                                    &dtgexc_obj->ifstref,
                                    &dtgexc_obj->ilstref
                                    );

    /* Compute libflame's Lapacke o/p  */
    dtgexc_obj->info = LAPACKE_dtgexc(  dtgexc_obj->matrix_layout,
                                    (lapack_logical)dtgexc_obj->wantq,
                                    (lapack_logical)dtgexc_obj->wantz, 
                                    dtgexc_obj->n,
                                    dtgexc_obj->a,
                                    dtgexc_obj->lda,
                                    dtgexc_obj->b,
                                    dtgexc_obj->ldb,
                                    dtgexc_obj->q,
                                    dtgexc_obj->ldq,
                                    dtgexc_obj->z,
                                    dtgexc_obj->ldz,
                                    &dtgexc_obj->ifst,
                                    &dtgexc_obj->ilst
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    dtgexc_obj->diff_a =  computeDiff_d( (dtgexc_obj->lda)*(dtgexc_obj->n), 
                dtgexc_obj->a, dtgexc_obj->aref );

    dtgexc_obj->diff_b =  computeDiff_d( (dtgexc_obj->ldb)*(dtgexc_obj->n), 
                dtgexc_obj->b, dtgexc_obj->bref );

    dtgexc_obj->diff_q =  computeDiff_d( (dtgexc_obj->ldq)*(dtgexc_obj->n), 
                dtgexc_obj->q, dtgexc_obj->qref );

    dtgexc_obj->diff_z =  computeDiff_d( (dtgexc_obj->ldz)*(dtgexc_obj->n), 
                dtgexc_obj->z, dtgexc_obj->zref );
}

TEST_F(dtgexc_test, dtgexc1) {
    EXPECT_NEAR(0.0, dtgexc_obj->diff_a, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_b, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_q, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_z, dtgexc_obj->threshold);
}

TEST_F(dtgexc_test, dtgexc2) {
    EXPECT_NEAR(0.0, dtgexc_obj->diff_a, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_b, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_q, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_z, dtgexc_obj->threshold);
}

TEST_F(dtgexc_test, dtgexc3) {
    EXPECT_NEAR(0.0, dtgexc_obj->diff_a, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_b, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_q, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_z, dtgexc_obj->threshold);
}

TEST_F(dtgexc_test, dtgexc4) {
    EXPECT_NEAR(0.0, dtgexc_obj->diff_a, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_b, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_q, dtgexc_obj->threshold);
    EXPECT_NEAR(0.0, dtgexc_obj->diff_z, dtgexc_obj->threshold);
}

/* Begin tgexc_scomplex_common_parameters  class definition */
class tgexc_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_a, diff_b, diff_q, diff_z;
	float threshold;
    void *hModule, *dModule;
    
	// Intermediate buffers to generate q, z i/ps for tgexc API
    lapack_complex_float* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    float* d, *e;
    lapack_complex_float* tauq, *taup;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    
    char wantq; // Must be 'N', 'I', or 'V'.
    char wantz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int ifst, ifstref; // ifst and ilst mark the rows and columns of A which are to be reduced.
    lapack_int ilst, ilstref;
    
    /* Input / Output parameters */
    
    lapack_complex_float* a, *aref; // contains the n-by-n general matrix A.
    lapack_int lda; //  The leading dimension of a
    lapack_complex_float* b, *bref; // contains the n-by-n upper triangular matrix B.
    lapack_int ldb; //  The leading dimension of b
    lapack_complex_float* q, *qref; // If wantq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    lapack_complex_float* z, *zref; //
    lapack_int ldz; // The leading dimension of z;
    
    /*Return Values */
    int info, inforef;

   public:
      tgexc_scomplex_parameters (int matrix_layout_i, char wantq, char wantz,
                              lapack_int n);

      ~tgexc_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
tgexc_scomplex_parameters:: tgexc_scomplex_parameters (int matrix_layout_i,
                            char wantq_i, char wantz_i, lapack_int n_i)
{
	int randIndex1, randIndex2;

    matrix_layout = matrix_layout_i;
    wantq = wantq_i;
    wantz = wantz_i;
    n = n_i;

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
   printf(" \n tgexc lapack_complex_float: matrix_layout: %d n: %d  wantz: %c \
       wantq: %c  \n", matrix_layout, n, wantz, wantq);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &z, &zref, ldz*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &taup, &tauq, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (b1Q==NULL) || (b1Z==NULL) || \
        (z==NULL) || (zref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "tgexc_scomplex_parameters object: malloc error. Exiting ";
       tgexc_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( b1Q, b1Z, ldb*n );
	memcpy (b1Q, b, sizeof(lapack_complex_float)*ldb*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( q, qref, ldq*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( z, zref, ldz*n );
	lapacke_gtest_init_float_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_scomplex_buffer_pair_with_constant( taup, tauq, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'tgexc_scomplex_common_parameters' */
tgexc_scomplex_parameters :: ~tgexc_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgexc_free();
} 

//  Test fixture class definition
class ctgexc_test  : public  ::testing::Test {
public:
   tgexc_scomplex_parameters  *ctgexc_obj;
   void SetUp();  
   void TearDown () { delete ctgexc_obj; }
};

void ctgexc_test::SetUp()
{
    /* LAPACKE CTGEXC prototype */
    typedef int (*Fptr_NL_LAPACKE_ctgexc) ( int matrix_layout, 
		lapack_logical wantq, lapack_logical wantz, lapack_int n,
		lapack_complex_float* a, lapack_int lda, lapack_complex_float* b,
		lapack_int ldb, lapack_complex_float* q, lapack_int ldq,
		lapack_complex_float* z, lapack_int ldz, lapack_int ifst,
		lapack_int ilst);

    Fptr_NL_LAPACKE_ctgexc CTGEXC;

    ctgexc_obj = new  tgexc_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].wantz,
                                         eig_non_sym_paramslist[idx].wantq,
                                         eig_paramslist[idx].n );
										 
    ctgexc_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    ctgexc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctgexc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctgexc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctgexc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CTGEXC = (Fptr_NL_LAPACKE_ctgexc)dlsym(ctgexc_obj->hModule, "LAPACKE_ctgexc");
    ASSERT_TRUE(CTGEXC != NULL) << "failed to ppt the Netlib LAPACKE_ctgexc symbol";
 
 	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    ctgexc_obj->info = LAPACKE_cgebrd(  ctgexc_obj->matrix_layout,
                                    ctgexc_obj->n,
                                    ctgexc_obj->n,
                                    ctgexc_obj->b1Z,
                                    ctgexc_obj->ldb,
                                    ctgexc_obj->d,
                                    ctgexc_obj->e,
                                    ctgexc_obj->tauq,
                                    ctgexc_obj->taup
                                    );

    ctgexc_obj->info = LAPACKE_cungbr(  ctgexc_obj->matrix_layout,
                                    'Q',
                                    ctgexc_obj->n,
                                    ctgexc_obj->n,
                                    ctgexc_obj->n,
                                    ctgexc_obj->b1Z,
                                    ctgexc_obj->ldb,
                                    ctgexc_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (ctgexc_obj->z, ctgexc_obj->b1Z, sizeof(lapack_complex_float)*(ctgexc_obj->ldb)*(ctgexc_obj->n) );
	memcpy (ctgexc_obj->zref, ctgexc_obj->b1Z, sizeof(lapack_complex_float)*(ctgexc_obj->ldb)*(ctgexc_obj->n) );

    ctgexc_obj->info = LAPACKE_cgebrd(  ctgexc_obj->matrix_layout,
                                    ctgexc_obj->n,
                                    ctgexc_obj->n,
                                    ctgexc_obj->b1Q,
                                    ctgexc_obj->ldb,
                                    ctgexc_obj->d,
                                    ctgexc_obj->e,
                                    ctgexc_obj->tauq,
                                    ctgexc_obj->taup
                                    );

    ctgexc_obj->info = LAPACKE_cungbr(  ctgexc_obj->matrix_layout,
                                    'Q',
                                    ctgexc_obj->n,
                                    ctgexc_obj->n,
                                    ctgexc_obj->n,
                                    ctgexc_obj->b1Q,
                                    ctgexc_obj->ldb,
                                    ctgexc_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (ctgexc_obj->q, ctgexc_obj->b1Q, sizeof(lapack_complex_float)*(ctgexc_obj->ldb)*(ctgexc_obj->n) );
	memcpy (ctgexc_obj->qref, ctgexc_obj->b1Q, sizeof(lapack_complex_float)*(ctgexc_obj->ldb)*(ctgexc_obj->n) );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    ctgexc_obj->inforef = CTGEXC(   ctgexc_obj->matrix_layout,
                                    (lapack_logical)ctgexc_obj->wantq,
                                    (lapack_logical)ctgexc_obj->wantz, 
                                    ctgexc_obj->n,
                                    ctgexc_obj->aref,
                                    ctgexc_obj->lda,
                                    ctgexc_obj->bref,
                                    ctgexc_obj->ldb,
                                    ctgexc_obj->qref,
                                    ctgexc_obj->ldq,
                                    ctgexc_obj->zref,
                                    ctgexc_obj->ldz,
                                    ctgexc_obj->ifstref,
                                    ctgexc_obj->ilstref
                                    );

    /* Compute libflame's Lapacke o/p  */
    ctgexc_obj->info = LAPACKE_ctgexc(  ctgexc_obj->matrix_layout,
                                    (lapack_logical)ctgexc_obj->wantq,
                                    (lapack_logical)ctgexc_obj->wantz, 
                                    ctgexc_obj->n,
                                    ctgexc_obj->a,
                                    ctgexc_obj->lda,
                                    ctgexc_obj->b,
                                    ctgexc_obj->ldb,
                                    ctgexc_obj->q,
                                    ctgexc_obj->ldq,
                                    ctgexc_obj->z,
                                    ctgexc_obj->ldz,
                                    ctgexc_obj->ifst,
                                    ctgexc_obj->ilst
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    ctgexc_obj->diff_a =  computeDiff_c( (ctgexc_obj->lda)*(ctgexc_obj->n), 
                ctgexc_obj->a, ctgexc_obj->aref );

    ctgexc_obj->diff_b =  computeDiff_c( (ctgexc_obj->ldb)*(ctgexc_obj->n), 
                ctgexc_obj->b, ctgexc_obj->bref );

    ctgexc_obj->diff_q =  computeDiff_c( (ctgexc_obj->ldq)*(ctgexc_obj->n), 
                ctgexc_obj->q, ctgexc_obj->qref );

    ctgexc_obj->diff_z =  computeDiff_c( (ctgexc_obj->ldz)*(ctgexc_obj->n), 
                ctgexc_obj->z, ctgexc_obj->zref );
}

TEST_F(ctgexc_test, ctgexc1) {
    EXPECT_NEAR(0.0, ctgexc_obj->diff_a, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_b, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_q, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_z, ctgexc_obj->threshold);
}

TEST_F(ctgexc_test, ctgexc2) {
    EXPECT_NEAR(0.0, ctgexc_obj->diff_a, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_b, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_q, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_z, ctgexc_obj->threshold);
}

TEST_F(ctgexc_test, ctgexc3) {
    EXPECT_NEAR(0.0, ctgexc_obj->diff_a, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_b, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_q, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_z, ctgexc_obj->threshold);
}

TEST_F(ctgexc_test, ctgexc4) {
    EXPECT_NEAR(0.0, ctgexc_obj->diff_a, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_b, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_q, ctgexc_obj->threshold);
    EXPECT_NEAR(0.0, ctgexc_obj->diff_z, ctgexc_obj->threshold);
}


/* Begin tgexc_dcomplex_common_parameters  class definition */
class tgexc_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_a, diff_b, diff_q, diff_z;
	float threshold;
    void *hModule, *dModule;
    
	// Intermediate buffers to generate q, z i/ps for tgexc API
    lapack_complex_double* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    double* d, *e;
    lapack_complex_double* tauq, *taup;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    
    char wantq; // Must be 'N', 'I', or 'V'.
    char wantz; // Must be 'N', 'I', or 'V'.
    lapack_int n; // The order of the matrices A and B (n≥ 0).
    lapack_int ifst, ifstref; // ifst and ilst mark the rows and columns of A which are to be reduced.
    lapack_int ilst, ilstref;
    
    /* Input / Output parameters */
    
    lapack_complex_double* a, *aref; // contains the n-by-n general matrix A.
    lapack_int lda; //  The leading dimension of a
    lapack_complex_double* b, *bref; // contains the n-by-n upper triangular matrix B.
    lapack_int ldb; //  The leading dimension of b
    lapack_complex_double* q, *qref; // If wantq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    lapack_complex_double* z, *zref; //
    lapack_int ldz; // The leading dimension of z;
    
    /*Return Values */
    int info, inforef;

   public:
      tgexc_dcomplex_parameters (int matrix_layout_i, char wantq, char wantz,
                              lapack_int n);

      ~tgexc_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
tgexc_dcomplex_parameters:: tgexc_dcomplex_parameters (int matrix_layout_i,
                            char wantq_i, char wantz_i, lapack_int n_i)
{
	int randIndex1, randIndex2;

    matrix_layout = matrix_layout_i;
    wantq = wantq_i;
    wantz = wantz_i;
    n = n_i;

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
   printf(" \n tgexc lapack_complex_double: matrix_layout: %d n: %d  wantz: %c \
       wantq: %c  \n", matrix_layout, n, wantz, wantq);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, ldb*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &q, &qref, ldq*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &z, &zref, ldz*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &taup, &tauq, n );

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) || \
        (b1Q==NULL) || (b1Z==NULL) || \
        (z==NULL) || (zref==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) || \
        (q==NULL) || (qref==NULL) ){
       EXPECT_FALSE( true) << "tgexc_dcomplex_parameters object: malloc error. Exiting ";
       tgexc_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, lda*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, ldb*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b1Q, b1Z, ldb*n );
	memcpy (b1Q, b, sizeof(lapack_complex_double)*ldb*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( q, qref, ldq*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( z, zref, ldz*n );
	lapacke_gtest_init_double_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_dcomplex_buffer_pair_with_constant( taup, tauq, n, 0.0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'tgexc_dcomplex_common_parameters' */
tgexc_dcomplex_parameters :: ~tgexc_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   tgexc_free();
} 

//  Test fixture class definition
class ztgexc_test  : public  ::testing::Test {
public:
   tgexc_dcomplex_parameters  *ztgexc_obj;
   void SetUp();  
   void TearDown () { delete ztgexc_obj; }
};

void ztgexc_test::SetUp()
{
    /* LAPACKE ZTGEXC prototype */
    typedef int (*Fptr_NL_LAPACKE_ztgexc) ( int matrix_layout, 
		lapack_logical wantq, lapack_logical wantz, lapack_int n,
		lapack_complex_double* a, lapack_int lda, lapack_complex_double* b,
		lapack_int ldb, lapack_complex_double* q, lapack_int ldq,
		lapack_complex_double* z, lapack_int ldz, lapack_int ifst,
		lapack_int ilst);

    Fptr_NL_LAPACKE_ztgexc ZTGEXC;

    ztgexc_obj = new  tgexc_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_non_sym_paramslist[idx].wantz,
                                         eig_non_sym_paramslist[idx].wantq,
                                         eig_paramslist[idx].n );
										 
    ztgexc_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    idx = Circular_Increment_Index(idx);
    ztgexc_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztgexc_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztgexc_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztgexc_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZTGEXC = (Fptr_NL_LAPACKE_ztgexc)dlsym(ztgexc_obj->hModule, "LAPACKE_ztgexc");
    ASSERT_TRUE(ZTGEXC != NULL) << "failed to ppt the Netlib LAPACKE_ztgexc symbol";
 
 	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    ztgexc_obj->info = LAPACKE_zgebrd(  ztgexc_obj->matrix_layout,
                                    ztgexc_obj->n,
                                    ztgexc_obj->n,
                                    ztgexc_obj->b1Z,
                                    ztgexc_obj->ldb,
                                    ztgexc_obj->d,
                                    ztgexc_obj->e,
                                    ztgexc_obj->tauq,
                                    ztgexc_obj->taup
                                    );

    ztgexc_obj->info = LAPACKE_zungbr(  ztgexc_obj->matrix_layout,
                                    'Q',
                                    ztgexc_obj->n,
                                    ztgexc_obj->n,
                                    ztgexc_obj->n,
                                    ztgexc_obj->b1Z,
                                    ztgexc_obj->ldb,
                                    ztgexc_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (ztgexc_obj->z, ztgexc_obj->b1Z, sizeof(lapack_complex_double)*(ztgexc_obj->ldb)*(ztgexc_obj->n) );
	memcpy (ztgexc_obj->zref, ztgexc_obj->b1Z, sizeof(lapack_complex_double)*(ztgexc_obj->ldb)*(ztgexc_obj->n) );

    ztgexc_obj->info = LAPACKE_zgebrd(  ztgexc_obj->matrix_layout,
                                    ztgexc_obj->n,
                                    ztgexc_obj->n,
                                    ztgexc_obj->b1Q,
                                    ztgexc_obj->ldb,
                                    ztgexc_obj->d,
                                    ztgexc_obj->e,
                                    ztgexc_obj->tauq,
                                    ztgexc_obj->taup
                                    );

    ztgexc_obj->info = LAPACKE_zungbr(  ztgexc_obj->matrix_layout,
                                    'Q',
                                    ztgexc_obj->n,
                                    ztgexc_obj->n,
                                    ztgexc_obj->n,
                                    ztgexc_obj->b1Q,
                                    ztgexc_obj->ldb,
                                    ztgexc_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (ztgexc_obj->q, ztgexc_obj->b1Q, sizeof(lapack_complex_double)*(ztgexc_obj->ldb)*(ztgexc_obj->n) );
	memcpy (ztgexc_obj->qref, ztgexc_obj->b1Q, sizeof(lapack_complex_double)*(ztgexc_obj->ldb)*(ztgexc_obj->n) );

    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    ztgexc_obj->inforef = ZTGEXC(   ztgexc_obj->matrix_layout,
                                    (lapack_logical)ztgexc_obj->wantq,
                                    (lapack_logical)ztgexc_obj->wantz, 
                                    ztgexc_obj->n,
                                    ztgexc_obj->aref,
                                    ztgexc_obj->lda,
                                    ztgexc_obj->bref,
                                    ztgexc_obj->ldb,
                                    ztgexc_obj->qref,
                                    ztgexc_obj->ldq,
                                    ztgexc_obj->zref,
                                    ztgexc_obj->ldz,
                                    ztgexc_obj->ifstref,
                                    ztgexc_obj->ilstref
                                    );

    /* Compute libflame's Lapacke o/p  */
    ztgexc_obj->info = LAPACKE_ztgexc(  ztgexc_obj->matrix_layout,
                                    (lapack_logical)ztgexc_obj->wantq,
                                    (lapack_logical)ztgexc_obj->wantz, 
                                    ztgexc_obj->n,
                                    ztgexc_obj->a,
                                    ztgexc_obj->lda,
                                    ztgexc_obj->b,
                                    ztgexc_obj->ldb,
                                    ztgexc_obj->q,
                                    ztgexc_obj->ldq,
                                    ztgexc_obj->z,
                                    ztgexc_obj->ldz,
                                    ztgexc_obj->ifst,
                                    ztgexc_obj->ilst
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    /* If howmny = 'A' or 'S', then vl need not be set. */
    ztgexc_obj->diff_a =  computeDiff_z( (ztgexc_obj->lda)*(ztgexc_obj->n), 
                ztgexc_obj->a, ztgexc_obj->aref );

    ztgexc_obj->diff_b =  computeDiff_z( (ztgexc_obj->ldb)*(ztgexc_obj->n), 
                ztgexc_obj->b, ztgexc_obj->bref );

    ztgexc_obj->diff_q =  computeDiff_z( (ztgexc_obj->ldq)*(ztgexc_obj->n), 
                ztgexc_obj->q, ztgexc_obj->qref );

    ztgexc_obj->diff_z =  computeDiff_z( (ztgexc_obj->ldz)*(ztgexc_obj->n), 
                ztgexc_obj->z, ztgexc_obj->zref );
}

TEST_F(ztgexc_test, ztgexc1) {
    EXPECT_NEAR(0.0, ztgexc_obj->diff_a, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_b, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_q, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_z, ztgexc_obj->threshold);
}

TEST_F(ztgexc_test, ztgexc2) {
    EXPECT_NEAR(0.0, ztgexc_obj->diff_a, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_b, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_q, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_z, ztgexc_obj->threshold);
}

TEST_F(ztgexc_test, ztgexc3) {
    EXPECT_NEAR(0.0, ztgexc_obj->diff_a, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_b, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_q, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_z, ztgexc_obj->threshold);
}

TEST_F(ztgexc_test, ztgexc4) {
    EXPECT_NEAR(0.0, ztgexc_obj->diff_a, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_b, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_q, ztgexc_obj->threshold);
    EXPECT_NEAR(0.0, ztgexc_obj->diff_z, ztgexc_obj->threshold);
}
