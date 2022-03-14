#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define ggbak_free() \
    if (a!=NULL)        free(a); \
    if (b!=NULL)        free(b); \
    if (b1Q!=NULL)        free(b1Q); \
    if (b1Z!=NULL)     free(b1Z); \
    if (v!=NULL)        free(v); \
    if (vref!=NULL)     free(vref); \
    if (lscale!=NULL)   free(lscale); \
    if (rscale!=NULL)   free(rscale); \
    if (alphar!=NULL)     free(alphar); \
    if (alphai!=NULL)     free(alphai); \
    if (beta!=NULL)     free(beta); \
    if (z!=NULL)        free(z); \
    if (d!=NULL)        free(d); \
    if (taup!=NULL)     free(taup); \
    if (tauq!=NULL)     free(tauq); \
    if (e!=NULL)        free(e); \
    if (q!=NULL)        free(q);

#define LAPACKE_TEST_VERBOSE  (1)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ggbak_float_common_parameters  class definition */
class ggbak_float_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_v;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'N', 'P', 'S', or 'B'.
	char side; // Must be 'L', or 'R'.
     
	int mm, m_trevc; // The number of columns in the arrays vl and/or vr (mm≥m)
    lapack_int n, m; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b

	// Intermediate buffers to generate i/ps for hgeqz API
    float* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    float* d, *e;
    float* tauq, *taup;
	//int *select; // Specifies the eigenvectors to be computed.  


    /* Input / Output parameters */
    float* a; // contains the n-by-n general matrix A.
    float* b; // contains the n-by-n upper triangular matrix B.
    float *lscale;
    float *rscale;
    float* q; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    lapack_int ldv; // The leading dimension of vl, vr, v matrices;
    float* z; //
    lapack_int ldz; // The leading dimension of z;

    /* Input / Output parameters */
    lapack_int ilo; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi;
    float* alphar;
	float* alphai;
	float* beta;

    /* Output parameters */
    float* v, *vref; // contain the o/p eigen vectors from ggbak o/p.
    /*Return Values */
    int info, inforef;

    public:
      ggbak_float_parameters (int matrix_layout_i, char job, char side,
                         	  lapack_int n, lapack_int m);
      ~ggbak_float_parameters ();
};

/* Constructor definition  float_common_parameters */
ggbak_float_parameters:: ggbak_float_parameters (int matrix_layout_i, char job_i, 
                                      char side_i, lapack_int n_i, lapack_int m_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
	side = side_i;
    n = n_i;
    m = m_i;
	mm = m_i;

    lda = m;
    ldb = m;
    ldq = m;
    ldz = m;
	ldv = m;

    hModule = NULL;
    dModule = NULL;
	
    diff_v = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
	m_trevc = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbak float: matrix_layout: %d  job: %c \
   side: %c  n: %d  m: %d \n", matrix_layout, job, side, n, m);
#endif

    /* Memory allocation of the buffers */
	
    lapacke_gtest_alloc_float_buffer_pair( &a, &b, lda*n );
    lapacke_gtest_alloc_float_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &lscale, &rscale, n );
    lapacke_gtest_alloc_float_buffer_pair( &q, &z, ldq*n );
    lapacke_gtest_alloc_float_buffer_pair( &alphar, &alphai, n );
    lapacke_gtest_alloc_float_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_float_buffer_pair( &taup, &tauq, n );
    lapacke_gtest_alloc_float_buffer_pair( &v, &vref, n*ldv );
    beta = (float *)malloc (n * sizeof(float));

    if( (a==NULL) || (b==NULL) ||  \
        (b1Q==NULL) || (b1Z==NULL) || \
        (lscale==NULL) || (rscale==NULL) || \
        (alphar==NULL) || (alphai==NULL) || \
        (beta==NULL) ||  (q==NULL) || \
        (v==NULL)  || (vref==NULL) || \
        (z==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) )
    {
       EXPECT_FALSE( true) << "ggbak_float_parameters object: malloc error. Exiting ";
       ggbak_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_rand( a, lda*n );
    lapacke_gtest_init_float_buffer_pair_rand( b, b1Q, ldb*n );
    lapacke_gtest_init_float_buffer_rand( b1Z, ldb*n );
	lapacke_gtest_init_float_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_float_buffer_pair_with_constant( taup, tauq, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( lscale, rscale, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( q, z, ldq*n, 0 );
    lapacke_gtest_init_float_buffer_pair_with_constant( alphar, alphai, n, 0);
    lapacke_gtest_init_float_buffer_with_constant( beta, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant( v, vref, n*ldv, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'ggbak_float_common_parameters' */
ggbak_float_parameters :: ~ggbak_float_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggbak_free();
} 

//  Test fixture class definition
class sggbak_test  : public  ::testing::Test {
public:
   ggbak_float_parameters  *sggbak_obj;
   void SetUp();  
   void TearDown () { delete sggbak_obj; }
};

void sggbak_test::SetUp()
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

    /* LAPACKE_stgevc prototype */	
    typedef int (*Fptr_NL_LAPACKE_stgevc) (int matrix_layout, char side,
	char howmny, const lapack_logical* select, lapack_int n, const float* s, 
	lapack_int lds, const float* p, lapack_int ldp, float* vl, 
	lapack_int ldvl, float* vr, lapack_int ldvr, lapack_int mm, lapack_int* m);
				 
    Fptr_NL_LAPACKE_stgevc STGEVC;

    /* LAPACKE_sggbak prototype */	
    typedef int (*Fptr_NL_LAPACKE_sggbak) ( int matrix_layout, char job,
				char side, lapack_int n, lapack_int ilo, lapack_int ihi,
				const float* lscale, const float* rscale, lapack_int m,
				float* v, lapack_int ldv );
				 
    Fptr_NL_LAPACKE_sggbak SGGBAK;

    sggbak_obj = new  ggbak_float_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].n );
                                         
    sggbak_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    sggbak_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    idx = Circular_Increment_Index(idx);
    sggbak_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sggbak_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sggbak_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SGGBAL = (Fptr_NL_LAPACKE_sggbal)dlsym(sggbak_obj->hModule, "LAPACKE_sggbal");
    ASSERT_TRUE(SGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_sggbal symbol";

    SGGHRD = (Fptr_NL_LAPACKE_sgghrd)dlsym(sggbak_obj->hModule, "LAPACKE_sgghrd");
    ASSERT_TRUE(SGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_sgghrd symbol";
 
    SHGEQZ = (Fptr_NL_LAPACKE_shgeqz)dlsym(sggbak_obj->hModule, "LAPACKE_shgeqz");
    ASSERT_TRUE(SHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_sggbak symbol";

    STGEVC = (Fptr_NL_LAPACKE_stgevc)dlsym(sggbak_obj->hModule, "LAPACKE_stgevc");
    ASSERT_TRUE(STGEVC != NULL) << "failed to ppt the Netlib LAPACKE_stgevc symbol";

    SGGBAK = (Fptr_NL_LAPACKE_sggbak)dlsym(sggbak_obj->hModule, "LAPACKE_sggbak");
    ASSERT_TRUE(SGGBAK != NULL) << "failed to ppt the Netlib LAPACKE_sggbak symbol";

	/*  Prepare the inputs to the TGEVC API from the 'gebal', 'gghrd' and 'hgeqz' 
	    api sequence  */	
	
	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    sggbak_obj->info = LAPACKE_sgebrd(  sggbak_obj->matrix_layout,
                                    sggbak_obj->n,
                                    sggbak_obj->n,
                                    sggbak_obj->b1Z,
                                    sggbak_obj->ldb,
                                    sggbak_obj->d,
                                    sggbak_obj->e,
                                    sggbak_obj->tauq,
                                    sggbak_obj->taup
                                    );

    sggbak_obj->info = LAPACKE_sorgbr(  sggbak_obj->matrix_layout,
                                    'Q',
                                    sggbak_obj->n,
                                    sggbak_obj->n,
                                    sggbak_obj->n,
                                    sggbak_obj->b1Z,
                                    sggbak_obj->ldb,
                                    sggbak_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (sggbak_obj->z, sggbak_obj->b1Z, sizeof(float)*(sggbak_obj->ldb)*(sggbak_obj->n) );

    sggbak_obj->info = LAPACKE_sgebrd(  sggbak_obj->matrix_layout,
                                    sggbak_obj->n,
                                    sggbak_obj->n,
                                    sggbak_obj->b1Q,
                                    sggbak_obj->ldb,
                                    sggbak_obj->d,
                                    sggbak_obj->e,
                                    sggbak_obj->tauq,
                                    sggbak_obj->taup
                                    );

    sggbak_obj->info = LAPACKE_sorgbr(  sggbak_obj->matrix_layout,
                                    'Q',
                                    sggbak_obj->n,
                                    sggbak_obj->n,
                                    sggbak_obj->n,
                                    sggbak_obj->b1Q,
                                    sggbak_obj->ldb,
                                    sggbak_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (sggbak_obj->q, sggbak_obj->b1Q, sizeof(float)*(sggbak_obj->ldb)*(sggbak_obj->n) );

    /* Compute the common i/p by invoking Netlib-Lapack's API  */
    sggbak_obj->inforef = SGGBAL(   sggbak_obj->matrix_layout,
                                    sggbak_obj->job, //'B',
                                    sggbak_obj->n,
                                    sggbak_obj->a,
                                    sggbak_obj->lda,
                                    sggbak_obj->b,
                                    sggbak_obj->ldb,
                                    &sggbak_obj->ilo,
                                    &sggbak_obj->ihi,
                                    sggbak_obj->lscale,
                                    sggbak_obj->rscale
                                    );

    sggbak_obj->inforef = SGGHRD(   sggbak_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    sggbak_obj->n,
                                    sggbak_obj->ilo,
                                    sggbak_obj->ihi, 
                                    sggbak_obj->a,
                                    sggbak_obj->lda,
                                    sggbak_obj->b,
                                    sggbak_obj->ldb,
                                    sggbak_obj->q,
                                    sggbak_obj->ldq,
                                    sggbak_obj->z,
                                    sggbak_obj->ldz
                                    );

    sggbak_obj->inforef = SHGEQZ(   sggbak_obj->matrix_layout,
                                    'S',
                                    'V',
                                    'V',
                                    sggbak_obj->n,
                                    sggbak_obj->ilo,
                                    sggbak_obj->ihi, 
                                    sggbak_obj->a,
                                    sggbak_obj->lda,
                                    sggbak_obj->b,
                                    sggbak_obj->ldb,
									sggbak_obj->alphar,
									sggbak_obj->alphai,
									sggbak_obj->beta,
                                    sggbak_obj->q,
                                    sggbak_obj->ldq,
                                    sggbak_obj->z,
                                    sggbak_obj->ldz
                                    );
									
    sggbak_obj->inforef = STGEVC(   sggbak_obj->matrix_layout,
                                    'B',
                                    'B',
				                    NULL,
                                    sggbak_obj->n,
                                    sggbak_obj->a,
                                    sggbak_obj->lda,
                                    sggbak_obj->b,
                                    sggbak_obj->ldb,
                                    sggbak_obj->q,
                                    sggbak_obj->ldq,
                                    sggbak_obj->z,
                                    sggbak_obj->ldz,
									sggbak_obj->mm,
									&sggbak_obj->m_trevc
                                    );

    if( sggbak_obj->side == 'L')
	{
		memcpy ( sggbak_obj->v, sggbak_obj->q, sizeof(float)*sggbak_obj->ldq*sggbak_obj->n);
		memcpy ( sggbak_obj->vref, sggbak_obj->q, sizeof(float)*sggbak_obj->ldq*sggbak_obj->n);
	}
	else
	{
		memcpy ( sggbak_obj->v, sggbak_obj->z, sizeof(float)*sggbak_obj->ldz*sggbak_obj->n);
		memcpy ( sggbak_obj->vref, sggbak_obj->z, sizeof(float)*sggbak_obj->ldz*sggbak_obj->n);
	}
	
    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    sggbak_obj->inforef = SGGBAK(   sggbak_obj->matrix_layout,
                                    sggbak_obj->job,
                                    sggbak_obj->side,
                                    sggbak_obj->n,
                                    sggbak_obj->ilo,
                                    sggbak_obj->ihi,
                                    sggbak_obj->lscale,
                                    sggbak_obj->rscale,
                                    sggbak_obj->m,
                                    sggbak_obj->vref,
                                    sggbak_obj->ldv
                                    );

    /* Compute libflame's Lapacke o/p  */
    sggbak_obj->inforef = LAPACKE_sggbak(   sggbak_obj->matrix_layout,
                                    sggbak_obj->job,
                                    sggbak_obj->side,
                                    sggbak_obj->n,
                                    sggbak_obj->ilo,
                                    sggbak_obj->ihi,
                                    sggbak_obj->lscale,
                                    sggbak_obj->rscale,
                                    sggbak_obj->m,
                                    sggbak_obj->v,
                                    sggbak_obj->ldv
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    sggbak_obj->diff_v =  computeDiff_s( (sggbak_obj->ldv)*(sggbak_obj->n), 
                sggbak_obj->v, sggbak_obj->vref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbak float: \n diff_v: %f  ",sggbak_obj->diff_v );
#endif
}

TEST_F(sggbak_test, sggbak1) {
    EXPECT_NEAR(0.0, sggbak_obj->diff_v, sggbak_obj->threshold);
}

TEST_F(sggbak_test, sggbak2) {
    EXPECT_NEAR(0.0, sggbak_obj->diff_v, sggbak_obj->threshold);
}

TEST_F(sggbak_test, sggbak3) {
    EXPECT_NEAR(0.0, sggbak_obj->diff_v, sggbak_obj->threshold);
}

TEST_F(sggbak_test, sggbak4) {
    EXPECT_NEAR(0.0, sggbak_obj->diff_v, sggbak_obj->threshold);
}


/* Begin ggbak_double_common_parameters  class definition */
class ggbak_double_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_v;

    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'N', 'P', 'S', or 'B'.
	char side; // Must be 'L', or 'R'.
     
	int mm, m_trevc; // The number of columns in the arrays vl and/or vr (mm≥m)
    lapack_int n, m; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
	// Intermediate buffers to generate i/ps for hgeqz API
    double* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    double* d, *e;
    double* tauq, *taup;
	//int *select; // Specifies the eigenvectors to be computed.  


    /* Input / Output parameters */
    double* a; // contains the n-by-n general matrix A.
    double* b; // contains the n-by-n upper triangular matrix B.
    double *lscale;
    double *rscale;
    double* q; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    
    lapack_int ldq; // The leading dimension of q;
    lapack_int ldv; // The leading dimension of vl, vr, v matrices;
    double* z; //
    lapack_int ldz; // The leading dimension of z;

    /* Input / Output parameters */
    lapack_int ilo; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi;
    double* alphar;
	double* alphai;
	double* beta;

    /* Output parameters */
    double* v, *vref; // contain the o/p eigen vectors from ggbak o/p.
    /*Return Values */
    int info, inforef;

    public:
      ggbak_double_parameters (int matrix_layout_i, char job, char side,
                         	  lapack_int n, lapack_int m);
      ~ggbak_double_parameters ();
};

/* Constructor definition  double_common_parameters */
ggbak_double_parameters:: ggbak_double_parameters (int matrix_layout_i, char job_i, 
                                      char side_i, lapack_int n_i, lapack_int m_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
	side = side_i;
    n = n_i;
    m = m_i;
	mm = m_i;

    lda = m;
    ldb = m;
    ldq = m;
    ldz = m;
	ldv = m;

    hModule = NULL;
    dModule = NULL;
	
    diff_v = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
	m_trevc = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbak double: matrix_layout: %d  job: %c \
   side: %c  n: %d  m: %d \n", matrix_layout, job, side, n, m);
#endif

    /* Memory allocation of the buffers */
	
    lapacke_gtest_alloc_double_buffer_pair( &a, &b, lda*n );
    lapacke_gtest_alloc_double_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &lscale, &rscale, n );
    lapacke_gtest_alloc_double_buffer_pair( &q, &z, ldq*n );
    lapacke_gtest_alloc_double_buffer_pair( &alphar, &alphai, n );
    lapacke_gtest_alloc_double_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_double_buffer_pair( &taup, &tauq, n );
    lapacke_gtest_alloc_double_buffer_pair( &v, &vref, n*ldv );
    beta = (double *)malloc (n * sizeof(double));

    if( (a==NULL) || (b==NULL) ||  \
        (b1Q==NULL) || (b1Z==NULL) || \
        (lscale==NULL) || (rscale==NULL) || \
        (alphar==NULL) || (alphai==NULL) || \
        (beta==NULL) ||  (q==NULL) || \
        (v==NULL)  || (vref==NULL) || \
        (z==NULL) ||  \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) )
    {
       EXPECT_FALSE( true) << "ggbak_double_parameters object: malloc error. Exiting ";
       ggbak_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_rand( a, lda*n );
    lapacke_gtest_init_double_buffer_pair_rand( b, b1Q, ldb*n );
    lapacke_gtest_init_double_buffer_rand( b1Z, ldb*n );
	lapacke_gtest_init_double_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_double_buffer_pair_with_constant( taup, tauq, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( lscale, rscale, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( q, z, ldq*n, 0 );
    lapacke_gtest_init_double_buffer_pair_with_constant( alphar, alphai, n, 0);
    lapacke_gtest_init_double_buffer_with_constant( beta, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant( v, vref, n*ldv, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'ggbak_double_common_parameters' */
ggbak_double_parameters :: ~ggbak_double_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggbak_free();
} 

//  Test fixture class definition
class dggbak_test  : public  ::testing::Test {
public:
   ggbak_double_parameters  *dggbak_obj;
   void SetUp();  
   void TearDown () { delete dggbak_obj; }
};

void dggbak_test::SetUp()
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

    /* LAPACKE_dtgevc prototype */	
    typedef int (*Fptr_NL_LAPACKE_dtgevc) (int matrix_layout, char side,
	char howmny, const lapack_logical* select, lapack_int n, const double* s, 
	lapack_int lds, const double* p, lapack_int ldp, double* vl, 
	lapack_int ldvl, double* vr, lapack_int ldvr, lapack_int mm, lapack_int* m);
				 
    Fptr_NL_LAPACKE_dtgevc DTGEVC;

    /* LAPACKE_dggbak prototype */	
    typedef int (*Fptr_NL_LAPACKE_dggbak) ( int matrix_layout, char job,
				char side, lapack_int n, lapack_int ilo, lapack_int ihi,
				const double* lscale, const double* rscale, lapack_int m,
				double* v, lapack_int ldv );
				 
    Fptr_NL_LAPACKE_dggbak DGGBAK;

    dggbak_obj = new  ggbak_double_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].n );
                                         
    dggbak_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    dggbak_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    idx = Circular_Increment_Index(idx);
    dggbak_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dggbak_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dggbak_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DGGBAL = (Fptr_NL_LAPACKE_dggbal)dlsym(dggbak_obj->hModule, "LAPACKE_dggbal");
    ASSERT_TRUE(DGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_dggbal symbol";

    DGGHRD = (Fptr_NL_LAPACKE_dgghrd)dlsym(dggbak_obj->hModule, "LAPACKE_dgghrd");
    ASSERT_TRUE(DGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_dgghrd symbol";
 
    DHGEQZ = (Fptr_NL_LAPACKE_dhgeqz)dlsym(dggbak_obj->hModule, "LAPACKE_dhgeqz");
    ASSERT_TRUE(DHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_dggbak symbol";

    DTGEVC = (Fptr_NL_LAPACKE_dtgevc)dlsym(dggbak_obj->hModule, "LAPACKE_dtgevc");
    ASSERT_TRUE(DTGEVC != NULL) << "failed to ppt the Netlib LAPACKE_dtgevc symbol";

    DGGBAK = (Fptr_NL_LAPACKE_dggbak)dlsym(dggbak_obj->hModule, "LAPACKE_dggbak");
    ASSERT_TRUE(DGGBAK != NULL) << "failed to ppt the Netlib LAPACKE_dggbak symbol";

	/*  Prepare the inputs to the TGEVC API from the 'gebal', 'gghrd' and 'hgeqz' 
	    api sequence  */	
	
	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    dggbak_obj->info = LAPACKE_dgebrd(  dggbak_obj->matrix_layout,
                                    dggbak_obj->n,
                                    dggbak_obj->n,
                                    dggbak_obj->b1Z,
                                    dggbak_obj->ldb,
                                    dggbak_obj->d,
                                    dggbak_obj->e,
                                    dggbak_obj->tauq,
                                    dggbak_obj->taup
                                    );

    dggbak_obj->info = LAPACKE_dorgbr(  dggbak_obj->matrix_layout,
                                    'Q',
                                    dggbak_obj->n,
                                    dggbak_obj->n,
                                    dggbak_obj->n,
                                    dggbak_obj->b1Z,
                                    dggbak_obj->ldb,
                                    dggbak_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (dggbak_obj->z, dggbak_obj->b1Z, sizeof(double)*(dggbak_obj->ldb)*(dggbak_obj->n) );

    dggbak_obj->info = LAPACKE_dgebrd(  dggbak_obj->matrix_layout,
                                    dggbak_obj->n,
                                    dggbak_obj->n,
                                    dggbak_obj->b1Q,
                                    dggbak_obj->ldb,
                                    dggbak_obj->d,
                                    dggbak_obj->e,
                                    dggbak_obj->tauq,
                                    dggbak_obj->taup
                                    );

    dggbak_obj->info = LAPACKE_dorgbr(  dggbak_obj->matrix_layout,
                                    'Q',
                                    dggbak_obj->n,
                                    dggbak_obj->n,
                                    dggbak_obj->n,
                                    dggbak_obj->b1Q,
                                    dggbak_obj->ldb,
                                    dggbak_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (dggbak_obj->q, dggbak_obj->b1Q, sizeof(double)*(dggbak_obj->ldb)*(dggbak_obj->n) );

    /* Compute the common i/p by invoking Netlib-Lapack's API  */
    dggbak_obj->inforef = DGGBAL(   dggbak_obj->matrix_layout,
                                    dggbak_obj->job, //'B',
                                    dggbak_obj->n,
                                    dggbak_obj->a,
                                    dggbak_obj->lda,
                                    dggbak_obj->b,
                                    dggbak_obj->ldb,
                                    &dggbak_obj->ilo,
                                    &dggbak_obj->ihi,
                                    dggbak_obj->lscale,
                                    dggbak_obj->rscale
                                    );

    dggbak_obj->inforef = DGGHRD(   dggbak_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    dggbak_obj->n,
                                    dggbak_obj->ilo,
                                    dggbak_obj->ihi, 
                                    dggbak_obj->a,
                                    dggbak_obj->lda,
                                    dggbak_obj->b,
                                    dggbak_obj->ldb,
                                    dggbak_obj->q,
                                    dggbak_obj->ldq,
                                    dggbak_obj->z,
                                    dggbak_obj->ldz
                                    );

    dggbak_obj->inforef = DHGEQZ(   dggbak_obj->matrix_layout,
                                    'S',
                                    'V',
                                    'V',
                                    dggbak_obj->n,
                                    dggbak_obj->ilo,
                                    dggbak_obj->ihi, 
                                    dggbak_obj->a,
                                    dggbak_obj->lda,
                                    dggbak_obj->b,
                                    dggbak_obj->ldb,
									dggbak_obj->alphar,
									dggbak_obj->alphai,
									dggbak_obj->beta,
                                    dggbak_obj->q,
                                    dggbak_obj->ldq,
                                    dggbak_obj->z,
                                    dggbak_obj->ldz
                                    );
									
    dggbak_obj->inforef = DTGEVC(   dggbak_obj->matrix_layout,
                                    'B',
                                    'B',
				                    NULL,
                                    dggbak_obj->n,
                                    dggbak_obj->a,
                                    dggbak_obj->lda,
                                    dggbak_obj->b,
                                    dggbak_obj->ldb,
                                    dggbak_obj->q,
                                    dggbak_obj->ldq,
                                    dggbak_obj->z,
                                    dggbak_obj->ldz,
									dggbak_obj->mm,
									&dggbak_obj->m_trevc
                                    );

    if( dggbak_obj->side == 'L')
	{
		memcpy ( dggbak_obj->v, dggbak_obj->q, sizeof(double)*dggbak_obj->ldq*dggbak_obj->n);
		memcpy ( dggbak_obj->vref, dggbak_obj->q, sizeof(double)*dggbak_obj->ldq*dggbak_obj->n);
	}
	else
	{
		memcpy ( dggbak_obj->v, dggbak_obj->z, sizeof(double)*dggbak_obj->ldz*dggbak_obj->n);
		memcpy ( dggbak_obj->vref, dggbak_obj->z, sizeof(double)*dggbak_obj->ldz*dggbak_obj->n);
	}
	
    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    dggbak_obj->inforef = DGGBAK(   dggbak_obj->matrix_layout,
                                    dggbak_obj->job,
                                    dggbak_obj->side,
                                    dggbak_obj->n,
                                    dggbak_obj->ilo,
                                    dggbak_obj->ihi,
                                    dggbak_obj->lscale,
                                    dggbak_obj->rscale,
                                    dggbak_obj->m,
                                    dggbak_obj->vref,
                                    dggbak_obj->ldv
                                    );

    /* Compute libflame's Lapacke o/p  */
    dggbak_obj->inforef = LAPACKE_dggbak(   dggbak_obj->matrix_layout,
                                    dggbak_obj->job,
                                    dggbak_obj->side,
                                    dggbak_obj->n,
                                    dggbak_obj->ilo,
                                    dggbak_obj->ihi,
                                    dggbak_obj->lscale,
                                    dggbak_obj->rscale,
                                    dggbak_obj->m,
                                    dggbak_obj->v,
                                    dggbak_obj->ldv
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    dggbak_obj->diff_v =  computeDiff_d( (dggbak_obj->ldv)*(dggbak_obj->n), 
                dggbak_obj->v, dggbak_obj->vref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbak double: \n diff_v: %f  ",dggbak_obj->diff_v );
#endif
}

TEST_F(dggbak_test, dggbak1) {
    EXPECT_NEAR(0.0, dggbak_obj->diff_v, dggbak_obj->threshold);
}

TEST_F(dggbak_test, dggbak2) {
    EXPECT_NEAR(0.0, dggbak_obj->diff_v, dggbak_obj->threshold);
}

TEST_F(dggbak_test, dggbak3) {
    EXPECT_NEAR(0.0, dggbak_obj->diff_v, dggbak_obj->threshold);
}

TEST_F(dggbak_test, dggbak4) {
    EXPECT_NEAR(0.0, dggbak_obj->diff_v, dggbak_obj->threshold);
}

/* Begin ggbak_scomplex_common_parameters  class definition */
class ggbak_scomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    float diff_v;
    float threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'N', 'P', 'S', or 'B'.
	char side; // Must be 'L', or 'R'.
     
	int mm, m_trevc; // The number of columns in the arrays vl and/or vr (mm≥m)
    lapack_int n, m; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
	// Intermediate buffers to generate i/ps for hgeqz API
    lapack_complex_float* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    float* d, *e;
    lapack_complex_float* tauq, *taup;


    /* Input / Output parameters */
    lapack_complex_float* a; // contains the n-by-n general matrix A.
    lapack_complex_float* b; // contains the n-by-n upper triangular matrix B.
    float *lscale;
    float *rscale;
    lapack_complex_float* q; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    lapack_complex_float* v, *vref; // contain the o/p eigen vectors from ggbak o/p.
    lapack_complex_float *vl, *vr; // left & right eigenvectors.
    
    lapack_int ldq; // The leading dimension of q;
    lapack_int ldv; // The leading dimension of vl, vr, v matrices;
    lapack_complex_float* z; //
    lapack_int ldz; // The leading dimension of z;

    /* Input / Output parameters */
    lapack_int ilo; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi;

    /* Output parameters */
    lapack_complex_float* alphar;
	lapack_complex_float* alphai;
	lapack_complex_float* beta;
    /*Return Values */
    int info, inforef;

    public:
      ggbak_scomplex_parameters (int matrix_layout_i, char job, char side,
                         	  lapack_int n, lapack_int m);
      ~ggbak_scomplex_parameters ();
};

/* Constructor definition  lapack_complex_float_common_parameters */
ggbak_scomplex_parameters:: ggbak_scomplex_parameters (int matrix_layout_i, char job_i, 
                                      char side_i, lapack_int n_i, lapack_int m_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
	side = side_i;
    n = n_i;
    m = m_i;
	mm = m_i;

    lda = m;
    ldb = m;
    ldq = m;
    ldz = m;
	ldv = m;

    hModule = NULL;
    dModule = NULL;
	
    diff_v = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
	m_trevc = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbak lapack_complex_float: matrix_layout: %d  job: %c \
   side: %c  n: %d  m: %d \n", matrix_layout, job, side, n, m);
#endif

    /* Memory allocation of the buffers */
	
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &b, lda*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_float_buffer_pair( &lscale, &rscale, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &q, &z, ldq*n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &alphar, &alphai, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &vl, &vr, n*ldv );
    lapacke_gtest_alloc_float_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &taup, &tauq, n );
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &v, &vref, n*ldv );

    beta = (lapack_complex_float *)malloc (n * sizeof(lapack_complex_float));
	
    if( (a==NULL) || (b==NULL) ||  \
        (b1Q==NULL) || (b1Z==NULL) || \
        (lscale==NULL) || (rscale==NULL) || \
        (alphar==NULL) || (alphai==NULL) || \
        (beta==NULL) ||  (q==NULL) || \
        (v==NULL)  || (vref==NULL) || \
        (z==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) )
    {
       EXPECT_FALSE( true) << "ggbak_scomplex_parameters object: malloc error. Exiting ";
       ggbak_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_rand( a, lda*n );
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, b1Q, ldb*n );
    lapacke_gtest_init_scomplex_buffer_rand( b1Z, ldb*n );
	lapacke_gtest_init_float_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_scomplex_buffer_pair_with_constant( taup, tauq, n, 0.0);
    lapacke_gtest_init_float_buffer_pair_with_constant( lscale, rscale, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( q, z, ldq*n, 0 );
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( alphar, alphai, n, 0);
    lapacke_gtest_init_scomplex_buffer_with_constant( beta, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant( v, vref, n*ldv, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'ggbak_scomplex_common_parameters' */
ggbak_scomplex_parameters :: ~ggbak_scomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggbak_free();
} 

//  Test fixture class definition
class cggbak_test  : public  ::testing::Test {
public:
   ggbak_scomplex_parameters  *cggbak_obj;
   void SetUp();  
   void TearDown () { delete cggbak_obj; }
};

void cggbak_test::SetUp()
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

    /* LAPACKE_cggbak prototype */	
    typedef int (*Fptr_NL_LAPACKE_cggbak) ( int matrix_layout, char job,
				char side, lapack_int n, lapack_int ilo, lapack_int ihi,
				const float* lscale, const float* rscale, lapack_int m,
				lapack_complex_float* v, lapack_int ldv );
				 
    Fptr_NL_LAPACKE_cggbak CGGBAK;

    cggbak_obj = new  ggbak_scomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].n );
                                         
    cggbak_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    cggbak_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    idx = Circular_Increment_Index(idx);
    cggbak_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cggbak_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cggbak_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CGGBAL = (Fptr_NL_LAPACKE_cggbal)dlsym(cggbak_obj->hModule, "LAPACKE_cggbal");
    ASSERT_TRUE(CGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_cggbal symbol";

    CGGHRD = (Fptr_NL_LAPACKE_cgghrd)dlsym(cggbak_obj->hModule, "LAPACKE_cgghrd");
    ASSERT_TRUE(CGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_cgghrd symbol";
 
    CHGEQZ = (Fptr_NL_LAPACKE_chgeqz)dlsym(cggbak_obj->hModule, "LAPACKE_chgeqz");
    ASSERT_TRUE(CHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_cggbak symbol";

    CTGEVC = (Fptr_NL_LAPACKE_ctgevc)dlsym(cggbak_obj->hModule, "LAPACKE_ctgevc");
    ASSERT_TRUE(CTGEVC != NULL) << "failed to ppt the Netlib LAPACKE_ctgevc symbol";

    CGGBAK = (Fptr_NL_LAPACKE_cggbak)dlsym(cggbak_obj->hModule, "LAPACKE_cggbak");
    ASSERT_TRUE(CGGBAK != NULL) << "failed to ppt the Netlib LAPACKE_cggbak symbol";

	/*  Prepare the inputs to the TGEVC API from the 'gebal', 'gghrd' and 'hgeqz' 
	    api sequence  */	
	
	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    cggbak_obj->info = LAPACKE_cgebrd(  cggbak_obj->matrix_layout,
                                    cggbak_obj->n,
                                    cggbak_obj->n,
                                    cggbak_obj->b1Z,
                                    cggbak_obj->ldb,
                                    cggbak_obj->d,
                                    cggbak_obj->e,
                                    cggbak_obj->tauq,
                                    cggbak_obj->taup
                                    );

    cggbak_obj->info = LAPACKE_cungbr(  cggbak_obj->matrix_layout,
                                    'Q',
                                    cggbak_obj->n,
                                    cggbak_obj->n,
                                    cggbak_obj->n,
                                    cggbak_obj->b1Z,
                                    cggbak_obj->ldb,
                                    cggbak_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (cggbak_obj->z, cggbak_obj->b1Z, sizeof(lapack_complex_float)*(cggbak_obj->ldb)*(cggbak_obj->n) );

    cggbak_obj->info = LAPACKE_cgebrd(  cggbak_obj->matrix_layout,
                                    cggbak_obj->n,
                                    cggbak_obj->n,
                                    cggbak_obj->b1Q,
                                    cggbak_obj->ldb,
                                    cggbak_obj->d,
                                    cggbak_obj->e,
                                    cggbak_obj->tauq,
                                    cggbak_obj->taup
                                    );

    cggbak_obj->info = LAPACKE_cungbr(  cggbak_obj->matrix_layout,
                                    'Q',
                                    cggbak_obj->n,
                                    cggbak_obj->n,
                                    cggbak_obj->n,
                                    cggbak_obj->b1Q,
                                    cggbak_obj->ldb,
                                    cggbak_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (cggbak_obj->q, cggbak_obj->b1Q, sizeof(lapack_complex_float)*(cggbak_obj->ldb)*(cggbak_obj->n) );

    /* Compute the common i/p by invoking Netlib-Lapack's API  */
    cggbak_obj->inforef = CGGBAL(   cggbak_obj->matrix_layout,
                                    cggbak_obj->job, //'B',
                                    cggbak_obj->n,
                                    cggbak_obj->a,
                                    cggbak_obj->lda,
                                    cggbak_obj->b,
                                    cggbak_obj->ldb,
                                    &cggbak_obj->ilo,
                                    &cggbak_obj->ihi,
                                    cggbak_obj->lscale,
                                    cggbak_obj->rscale
                                    );

    cggbak_obj->inforef = CGGHRD(   cggbak_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    cggbak_obj->n,
                                    cggbak_obj->ilo,
                                    cggbak_obj->ihi, 
                                    cggbak_obj->a,
                                    cggbak_obj->lda,
                                    cggbak_obj->b,
                                    cggbak_obj->ldb,
                                    cggbak_obj->q,
                                    cggbak_obj->ldq,
                                    cggbak_obj->z,
                                    cggbak_obj->ldz
                                    );

    cggbak_obj->inforef = CHGEQZ(   cggbak_obj->matrix_layout,
                                    'S',
                                    'V',
                                    'V',
                                    cggbak_obj->n,
                                    cggbak_obj->ilo,
                                    cggbak_obj->ihi, 
                                    cggbak_obj->a,
                                    cggbak_obj->lda,
                                    cggbak_obj->b,
                                    cggbak_obj->ldb,
									cggbak_obj->alphar,
									cggbak_obj->beta,
                                    cggbak_obj->q,
                                    cggbak_obj->ldq,
                                    cggbak_obj->z,
                                    cggbak_obj->ldz
                                    );
									
    cggbak_obj->inforef = CTGEVC(   cggbak_obj->matrix_layout,
                                    'B',
                                    'B',
				                    NULL,
                                    cggbak_obj->n,
                                    cggbak_obj->a,
                                    cggbak_obj->lda,
                                    cggbak_obj->b,
                                    cggbak_obj->ldb,
                                    cggbak_obj->q,
                                    cggbak_obj->ldq,
                                    cggbak_obj->z,
                                    cggbak_obj->ldz,
									cggbak_obj->mm,
									&cggbak_obj->m_trevc
                                    );

    if( cggbak_obj->side == 'L')
	{
		memcpy ( cggbak_obj->v, cggbak_obj->q, sizeof(lapack_complex_float)*cggbak_obj->ldq*cggbak_obj->n);
		memcpy ( cggbak_obj->vref, cggbak_obj->q, sizeof(lapack_complex_float)*cggbak_obj->ldq*cggbak_obj->n);
	}
	else
	{
		memcpy ( cggbak_obj->v, cggbak_obj->z, sizeof(lapack_complex_float)*cggbak_obj->ldz*cggbak_obj->n);
		memcpy ( cggbak_obj->vref, cggbak_obj->z, sizeof(lapack_complex_float)*cggbak_obj->ldz*cggbak_obj->n);
	}
	
    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    cggbak_obj->inforef = CGGBAK(   cggbak_obj->matrix_layout,
                                    cggbak_obj->job,
                                    cggbak_obj->side,
                                    cggbak_obj->n,
                                    cggbak_obj->ilo,
                                    cggbak_obj->ihi,
                                    cggbak_obj->lscale,
                                    cggbak_obj->rscale,
                                    cggbak_obj->m,
                                    cggbak_obj->vref,
                                    cggbak_obj->ldv
                                    );

    /* Compute libflame's Lapacke o/p  */
    cggbak_obj->inforef = LAPACKE_cggbak(   cggbak_obj->matrix_layout,
                                    cggbak_obj->job,
                                    cggbak_obj->side,
                                    cggbak_obj->n,
                                    cggbak_obj->ilo,
                                    cggbak_obj->ihi,
                                    cggbak_obj->lscale,
                                    cggbak_obj->rscale,
                                    cggbak_obj->m,
                                    cggbak_obj->v,
                                    cggbak_obj->ldv
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    cggbak_obj->diff_v =  computeDiff_c( (cggbak_obj->ldv)*(cggbak_obj->n), 
                cggbak_obj->v, cggbak_obj->vref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbak lapack_complex_float: \n diff_v: %f  ",cggbak_obj->diff_v );
#endif
}

TEST_F(cggbak_test, cggbak1) {
    EXPECT_NEAR(0.0, cggbak_obj->diff_v, cggbak_obj->threshold);
}

TEST_F(cggbak_test, cggbak2) {
    EXPECT_NEAR(0.0, cggbak_obj->diff_v, cggbak_obj->threshold);
}

TEST_F(cggbak_test, cggbak3) {
    EXPECT_NEAR(0.0, cggbak_obj->diff_v, cggbak_obj->threshold);
}

TEST_F(cggbak_test, cggbak4) {
    EXPECT_NEAR(0.0, cggbak_obj->diff_v, cggbak_obj->threshold);
}

/* Begin ggbak_dcomplex_common_parameters  class definition */
class ggbak_dcomplex_parameters{

   public:
    // variables to capture difference between ref o/p & libflame lapacke o/p.
    double diff_v;
    double threshold;
    void *hModule, *dModule;

    /* INPUT PARAMETERS */
    int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    char job; // Must be 'N', 'P', 'S', or 'B'.
	char side; // Must be 'L', or 'R'.
     
	int mm, m_trevc; // The number of columns in the arrays vl and/or vr (mm≥m)
    lapack_int n, m; // The order of the matrices A and B (n≥ 0).
    lapack_int lda; //  The leading dimension of a
    lapack_int ldb; //  The leading dimension of b
	// Intermediate buffers to generate i/ps for hgeqz API
    lapack_complex_double* b1Q, *b1Z; // contains the n-by-n upper triangular matrix B.
    double* d, *e;
    lapack_complex_double* tauq, *taup;


    /* Input / Output parameters */
    lapack_complex_double* a; // contains the n-by-n general matrix A.
    lapack_complex_double* b; // contains the n-by-n upper triangular matrix B.
    double *lscale;
    double *rscale;
    lapack_complex_double* q; // If compq = 'V', then q is the orthogonal/unitary matrix Q1
              // , typically from the QR factorization of B.
    lapack_complex_double* v, *vref; // contain the o/p eigen vectors from ggbak o/p.
    lapack_complex_double *vl, *vr; // left & right eigenvectors.
    
    lapack_int ldq; // The leading dimension of q;
    lapack_int ldv; // The leading dimension of vl, vr, v matrices;
    lapack_complex_double* z; //
    lapack_int ldz; // The leading dimension of z;

    /* Input / Output parameters */
    lapack_int ilo; // ilo and ihi mark the rows and columns of A which are to be reduced.
    lapack_int ihi;

    /* Output parameters */
    lapack_complex_double* alphar;
	lapack_complex_double* alphai;
	lapack_complex_double* beta;
    /*Return Values */
    int info, inforef;

    public:
      ggbak_dcomplex_parameters (int matrix_layout_i, char job, char side,
                         	  lapack_int n, lapack_int m);
      ~ggbak_dcomplex_parameters ();
};

/* Constructor definition  lapack_complex_double_common_parameters */
ggbak_dcomplex_parameters:: ggbak_dcomplex_parameters (int matrix_layout_i, char job_i, 
                                      char side_i, lapack_int n_i, lapack_int m_i)
{
    matrix_layout = matrix_layout_i;
    job = job_i;
	side = side_i;
    n = n_i;
    m = m_i;
	mm = m_i;

    lda = m;
    ldb = m;
    ldq = m;
    ldz = m;
	ldv = m;

    hModule = NULL;
    dModule = NULL;
	
    diff_v = 0;
    // Initialize 'ilo', 'ihi' to 0.
    ilo = 0;
    ihi = 0;
	m_trevc = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbak lapack_complex_double: matrix_layout: %d  job: %c \
   side: %c  n: %d  m: %d \n", matrix_layout, job, side, n, m);
#endif

    /* Memory allocation of the buffers */
	
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &b, lda*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b1Q, &b1Z, ldb*n );
    lapacke_gtest_alloc_double_buffer_pair( &lscale, &rscale, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &q, &z, ldq*n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &alphar, &alphai, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &vl, &vr, n*ldv );
    lapacke_gtest_alloc_double_buffer_pair( &d, &e, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &taup, &tauq, n );
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &v, &vref, n*ldv );

    beta = (lapack_complex_double *)malloc (n * sizeof(lapack_complex_double));
	
    if( (a==NULL) || (b==NULL) ||  \
        (b1Q==NULL) || (b1Z==NULL) || \
        (lscale==NULL) || (rscale==NULL) || \
        (alphar==NULL) || (alphai==NULL) || \
        (beta==NULL) ||  (q==NULL) || \
        (v==NULL)  || (vref==NULL) || \
        (z==NULL) || \
        (d==NULL) || (e==NULL) || \
        (taup==NULL) || (tauq==NULL) )
    {
       EXPECT_FALSE( true) << "ggbak_dcomplex_parameters object: malloc error. Exiting ";
       ggbak_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_rand( a, lda*n );
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, b1Q, ldb*n );
    lapacke_gtest_init_dcomplex_buffer_rand( b1Z, ldb*n );
	lapacke_gtest_init_double_buffer_pair_with_constant( d, e, n, 0.0);
	lapacke_gtest_init_dcomplex_buffer_pair_with_constant( taup, tauq, n, 0.0);
    lapacke_gtest_init_double_buffer_pair_with_constant( lscale, rscale, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( q, z, ldq*n, 0 );
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( alphar, alphai, n, 0);
    lapacke_gtest_init_dcomplex_buffer_with_constant( beta, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant( v, vref, n*ldv, 0);
    
   } /* end of Constructor  */
    

/* Destructor definition  'ggbak_dcomplex_common_parameters' */
ggbak_dcomplex_parameters :: ~ggbak_dcomplex_parameters ()
{
   /* De-Allocate memory for the input matrices */
   ggbak_free();
} 

//  Test fixture class definition
class zggbak_test  : public  ::testing::Test {
public:
   ggbak_dcomplex_parameters  *zggbak_obj;
   void SetUp();  
   void TearDown () { delete zggbak_obj; }
};

void zggbak_test::SetUp()
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

    /* LAPACKE_zggbak prototype */	
    typedef int (*Fptr_NL_LAPACKE_zggbak) ( int matrix_layout, char job,
				char side, lapack_int n, lapack_int ilo, lapack_int ihi,
				const double* lscale, const double* rscale, lapack_int m,
				lapack_complex_double* v, lapack_int ldv );
				 
    Fptr_NL_LAPACKE_zggbak ZGGBAK;

    zggbak_obj = new  ggbak_dcomplex_parameters( eig_paramslist[idx].matrix_layout,
                                         eig_paramslist[idx].job,
                                         eig_paramslist[idx].side,
                                         eig_paramslist[idx].n,
                                         eig_paramslist[idx].n );
                                         
    zggbak_obj->threshold = eig_non_sym_paramslist[idx].GenNonSymEigProblem_threshold;
    zggbak_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    idx = Circular_Increment_Index(idx);
    zggbak_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zggbak_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zggbak_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZGGBAL = (Fptr_NL_LAPACKE_zggbal)dlsym(zggbak_obj->hModule, "LAPACKE_zggbal");
    ASSERT_TRUE(ZGGBAL != NULL) << "failed to ppt the Netlib LAPACKE_zggbal symbol";

    ZGGHRD = (Fptr_NL_LAPACKE_zgghrd)dlsym(zggbak_obj->hModule, "LAPACKE_zgghrd");
    ASSERT_TRUE(ZGGHRD != NULL) << "failed to ppt the Netlib LAPACKE_zgghrd symbol";
 
    ZHGEQZ = (Fptr_NL_LAPACKE_zhgeqz)dlsym(zggbak_obj->hModule, "LAPACKE_zhgeqz");
    ASSERT_TRUE(ZHGEQZ != NULL) << "failed to ppt the Netlib LAPACKE_zggbak symbol";

    ZTGEVC = (Fptr_NL_LAPACKE_ztgevc)dlsym(zggbak_obj->hModule, "LAPACKE_ztgevc");
    ASSERT_TRUE(ZTGEVC != NULL) << "failed to ppt the Netlib LAPACKE_ztgevc symbol";

    ZGGBAK = (Fptr_NL_LAPACKE_zggbak)dlsym(zggbak_obj->hModule, "LAPACKE_zggbak");
    ASSERT_TRUE(ZGGBAK != NULL) << "failed to ppt the Netlib LAPACKE_zggbak symbol";

	/*  Prepare the inputs to the TGEVC API from the 'gebal', 'gghrd' and 'hgeqz' 
	    api sequence  */	
	
	/* Generate Q matrix for gghrd through gebrd & orgbr APIs */
    zggbak_obj->info = LAPACKE_zgebrd(  zggbak_obj->matrix_layout,
                                    zggbak_obj->n,
                                    zggbak_obj->n,
                                    zggbak_obj->b1Z,
                                    zggbak_obj->ldb,
                                    zggbak_obj->d,
                                    zggbak_obj->e,
                                    zggbak_obj->tauq,
                                    zggbak_obj->taup
                                    );

    zggbak_obj->info = LAPACKE_zungbr(  zggbak_obj->matrix_layout,
                                    'Q',
                                    zggbak_obj->n,
                                    zggbak_obj->n,
                                    zggbak_obj->n,
                                    zggbak_obj->b1Z,
                                    zggbak_obj->ldb,
                                    zggbak_obj->tauq
                                    );
    
	// copy the orthogonal matriz 'Z' for i/p to gghrd call.
	memcpy (zggbak_obj->z, zggbak_obj->b1Z, sizeof(lapack_complex_double)*(zggbak_obj->ldb)*(zggbak_obj->n) );

    zggbak_obj->info = LAPACKE_zgebrd(  zggbak_obj->matrix_layout,
                                    zggbak_obj->n,
                                    zggbak_obj->n,
                                    zggbak_obj->b1Q,
                                    zggbak_obj->ldb,
                                    zggbak_obj->d,
                                    zggbak_obj->e,
                                    zggbak_obj->tauq,
                                    zggbak_obj->taup
                                    );

    zggbak_obj->info = LAPACKE_zungbr(  zggbak_obj->matrix_layout,
                                    'Q',
                                    zggbak_obj->n,
                                    zggbak_obj->n,
                                    zggbak_obj->n,
                                    zggbak_obj->b1Q,
                                    zggbak_obj->ldb,
                                    zggbak_obj->tauq
                                    );
	// copy the orthogonal matriz 'Q' for i/p to gghrd call.
	memcpy (zggbak_obj->q, zggbak_obj->b1Q, sizeof(lapack_complex_double)*(zggbak_obj->ldb)*(zggbak_obj->n) );

    /* Compute the common i/p by invoking Netlib-Lapack's API  */
    zggbak_obj->inforef = ZGGBAL(   zggbak_obj->matrix_layout,
                                    zggbak_obj->job, //'B',
                                    zggbak_obj->n,
                                    zggbak_obj->a,
                                    zggbak_obj->lda,
                                    zggbak_obj->b,
                                    zggbak_obj->ldb,
                                    &zggbak_obj->ilo,
                                    &zggbak_obj->ihi,
                                    zggbak_obj->lscale,
                                    zggbak_obj->rscale
                                    );

    zggbak_obj->inforef = ZGGHRD(   zggbak_obj->matrix_layout,
                                    'V',
                                    'V', 
                                    zggbak_obj->n,
                                    zggbak_obj->ilo,
                                    zggbak_obj->ihi, 
                                    zggbak_obj->a,
                                    zggbak_obj->lda,
                                    zggbak_obj->b,
                                    zggbak_obj->ldb,
                                    zggbak_obj->q,
                                    zggbak_obj->ldq,
                                    zggbak_obj->z,
                                    zggbak_obj->ldz
                                    );

    zggbak_obj->inforef = ZHGEQZ(   zggbak_obj->matrix_layout,
                                    'S',
                                    'V',
                                    'V',
                                    zggbak_obj->n,
                                    zggbak_obj->ilo,
                                    zggbak_obj->ihi, 
                                    zggbak_obj->a,
                                    zggbak_obj->lda,
                                    zggbak_obj->b,
                                    zggbak_obj->ldb,
									zggbak_obj->alphar,
									zggbak_obj->beta,
                                    zggbak_obj->q,
                                    zggbak_obj->ldq,
                                    zggbak_obj->z,
                                    zggbak_obj->ldz
                                    );
									
    zggbak_obj->inforef = ZTGEVC(   zggbak_obj->matrix_layout,
                                    'B',
                                    'B',
				                    NULL,
                                    zggbak_obj->n,
                                    zggbak_obj->a,
                                    zggbak_obj->lda,
                                    zggbak_obj->b,
                                    zggbak_obj->ldb,
                                    zggbak_obj->q,
                                    zggbak_obj->ldq,
                                    zggbak_obj->z,
                                    zggbak_obj->ldz,
									zggbak_obj->mm,
									&zggbak_obj->m_trevc
                                    );

    if( zggbak_obj->side == 'L')
	{
		memcpy ( zggbak_obj->v, zggbak_obj->q, sizeof(lapack_complex_double)*zggbak_obj->ldq*zggbak_obj->n);
		memcpy ( zggbak_obj->vref, zggbak_obj->q, sizeof(lapack_complex_double)*zggbak_obj->ldq*zggbak_obj->n);
	}
	else
	{
		memcpy ( zggbak_obj->v, zggbak_obj->z, sizeof(lapack_complex_double)*zggbak_obj->ldz*zggbak_obj->n);
		memcpy ( zggbak_obj->vref, zggbak_obj->z, sizeof(lapack_complex_double)*zggbak_obj->ldz*zggbak_obj->n);
	}
	
    /* Compute the reference o/p by invoking Netlib-Lapack's API  */
    zggbak_obj->inforef = ZGGBAK(   zggbak_obj->matrix_layout,
                                    zggbak_obj->job,
                                    zggbak_obj->side,
                                    zggbak_obj->n,
                                    zggbak_obj->ilo,
                                    zggbak_obj->ihi,
                                    zggbak_obj->lscale,
                                    zggbak_obj->rscale,
                                    zggbak_obj->m,
                                    zggbak_obj->vref,
                                    zggbak_obj->ldv
                                    );

    /* Compute libflame's Lapacke o/p  */
    zggbak_obj->inforef = LAPACKE_zggbak(   zggbak_obj->matrix_layout,
                                    zggbak_obj->job,
                                    zggbak_obj->side,
                                    zggbak_obj->n,
                                    zggbak_obj->ilo,
                                    zggbak_obj->ihi,
                                    zggbak_obj->lscale,
                                    zggbak_obj->rscale,
                                    zggbak_obj->m,
                                    zggbak_obj->v,
                                    zggbak_obj->ldv
                                    );

    /* Capture the Netlib, libflame o/p buffers' differences */
    
    zggbak_obj->diff_v =  computeDiff_z( (zggbak_obj->ldv)*(zggbak_obj->n), 
                zggbak_obj->v, zggbak_obj->vref );


#if LAPACKE_TEST_VERBOSE
   printf(" \n ggbak lapack_complex_double: \n diff_v: %f  ",zggbak_obj->diff_v );
#endif
}

TEST_F(zggbak_test, zggbak1) {
    EXPECT_NEAR(0.0, zggbak_obj->diff_v, zggbak_obj->threshold);
}

TEST_F(zggbak_test, zggbak2) {
    EXPECT_NEAR(0.0, zggbak_obj->diff_v, zggbak_obj->threshold);
}

TEST_F(zggbak_test, zggbak3) {
    EXPECT_NEAR(0.0, zggbak_obj->diff_v, zggbak_obj->threshold);
}

TEST_F(zggbak_test, zggbak4) {
    EXPECT_NEAR(0.0, zggbak_obj->diff_v, zggbak_obj->threshold);
}

