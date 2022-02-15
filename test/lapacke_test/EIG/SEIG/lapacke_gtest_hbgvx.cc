#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hbgvx_free() \
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (bp!=NULL)  free(bp);\
if (bpref!=NULL)  free(bpref);\
if (z!=NULL)  free(z);\
if (zref!=NULL)  free(zref);\
if (ifail!=NULL)  free(ifail);\
if (q!=NULL)  free(q);\
if (qref!=NULL)  free(qref);\
if (ifailref!=NULL)  free(ifailref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hbgvx_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	int bufsize_z;
	int bufsize_q;
	void *hModule, *dModule;
	float diff_a;
	float diff_w;
	float diff_b;
	float diff_z;
	float vl, vu;
	float abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz,lda, ldb, ldq;
	lapack_int ka, kb;
	lapack_int il, iu;
	lapack_int* ifail;
	lapack_complex_float* ap, *bp, *z, *q;
	char uplo;
	char jobz;
	char range;
	/*Output Parameter*/
	float* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_float *apref, *bpref, *zref, *qref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hbgvx_scomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz, lapack_int ka);
      ~hbgvx_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hbgvx_scomplex_parameters:: hbgvx_scomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i, lapack_int ka_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	range = range_i;
	ka = ka_i;
	kb = ka;
	
	if (jobz == 'U')
		jobz = 'N';
	
	
	if ((range == 'N') || (range == 'V'))
	{	range = 'A';
		m = n;
	} else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;		
	} 
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hbgvx scomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
		lda = ka+1;
		ldb = kb+1;
		bufsize_a = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		lda = n;
		ldb = n;
		bufsize_a = lda *(ka+1);
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	ldq = n;
	bufsize_q = ldq *n;
	
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&bp, &bpref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&q, &qref, bufsize_q);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref,	n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(bp == NULL) || (bpref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)||
		(q == NULL) || (qref == NULL)){
		EXPECT_FALSE( true) << "hbgvx_float_parameters object: malloc error.";
		hbgvx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(ap, apref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand(bp, bpref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hbgvx_scomplex_parameters :: ~hbgvx_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbgvx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbgvx_free();
  

}
/*  Test fixture class definition */
class chbgvx_test  : public  ::testing::Test {
public:
   hbgvx_scomplex_parameters  *chbgvx_obj;
   void SetUp();
   void TearDown () { delete chbgvx_obj; }
};

void chbgvx_test::SetUp(){

    /* LAPACKE chbgvx prototype */
    typedef int (*Fptr_NL_LAPACKE_chbgvx) ( int matrix_layout, char jobz, char range,	char uplo, lapack_int n, lapack_int ka, 
										lapack_int kb, lapack_complex_float* ap, lapack_int lda, lapack_complex_float* bp, lapack_int ldb, lapack_complex_float* q, 
										lapack_int ldq, float vl, float vu,	lapack_int il, lapack_int iu, float abstol, lapack_int* m, float* w, lapack_complex_float* z,
										lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_chbgvx chbgvx;

    chbgvx_obj = new hbgvx_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m,
							eig_paramslist[idx].k);

    idx = Circular_Increment_Index(idx);

    chbgvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chbgvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chbgvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chbgvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chbgvx = (Fptr_NL_LAPACKE_chbgvx)dlsym(chbgvx_obj->hModule, "LAPACKE_chbgvx");
    ASSERT_TRUE(chbgvx != NULL) << "failed to get the Netlib LAPACKE_chbgvx symbol";
    

    chbgvx_obj->inforef = chbgvx( chbgvx_obj->matrix_layout,  chbgvx_obj->jobz, chbgvx_obj->range, chbgvx_obj->uplo, chbgvx_obj->n, 
	chbgvx_obj->ka, chbgvx_obj->kb, chbgvx_obj->apref, chbgvx_obj->lda, chbgvx_obj->bpref, chbgvx_obj->ldb, chbgvx_obj->qref, chbgvx_obj->ldq, 
	chbgvx_obj->vl, chbgvx_obj->vu, chbgvx_obj->il, chbgvx_obj->iu,	chbgvx_obj->abstol, &chbgvx_obj->m, chbgvx_obj->wref, chbgvx_obj->zref, 
	chbgvx_obj->ldz, chbgvx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    chbgvx_obj->info = LAPACKE_chbgvx( chbgvx_obj->matrix_layout, chbgvx_obj->jobz, chbgvx_obj->range, chbgvx_obj->uplo, chbgvx_obj->n, 
	chbgvx_obj->ka, chbgvx_obj->kb, chbgvx_obj->ap, chbgvx_obj->lda, chbgvx_obj->bp, chbgvx_obj->ldb, chbgvx_obj->q, chbgvx_obj->ldq, chbgvx_obj->vl, 
	chbgvx_obj->vu, chbgvx_obj->il, chbgvx_obj->iu,	chbgvx_obj->abstol, &chbgvx_obj->m, chbgvx_obj->w, chbgvx_obj->z, chbgvx_obj->ldz, chbgvx_obj->ifail);

    if( chbgvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chbgvx is wrong\n", chbgvx_obj->info );
    }
    if( chbgvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chbgvx is wrong\n", 
        chbgvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chbgvx_obj->diff_a =  computeDiff_c( chbgvx_obj->bufsize_a, chbgvx_obj->ap, chbgvx_obj->apref );
	chbgvx_obj->diff_b =  computeDiff_c( chbgvx_obj->bufsize_a, chbgvx_obj->bp, chbgvx_obj->bpref );
	chbgvx_obj->diff_z =  computeDiff_c( chbgvx_obj->bufsize_z, chbgvx_obj->z, chbgvx_obj->zref );
	chbgvx_obj->diff_w =  computeDiff_s( chbgvx_obj->n, chbgvx_obj->w, chbgvx_obj->wref );
}

TEST_F(chbgvx_test, chbgvx1) {
    EXPECT_NEAR(0.0, chbgvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chbgvx_test, chbgvx2) {
    EXPECT_NEAR(0.0, chbgvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chbgvx_test, chbgvx3) {
    EXPECT_NEAR(0.0, chbgvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chbgvx_test, chbgvx4) {
    EXPECT_NEAR(0.0, chbgvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hbgvx_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	int bufsize_z;
	int bufsize_q;
	void *hModule, *dModule;
	double diff_a;
	double diff_w;
	double diff_b;
	double diff_z;
	double vl, vu;
	double abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ka, kb;
	lapack_int ldz, lda, ldb, ldq;
	lapack_int il, iu;
	lapack_int* ifail;
	lapack_complex_double* ap, *bp, *z, *q;
	char uplo;
	char jobz;
	char range;
	/*Output Parameter*/
	double* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_double *apref, *bpref, *zref, *qref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hbgvx_dcomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz, lapack_int ka);
      ~hbgvx_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hbgvx_dcomplex_parameters:: hbgvx_dcomplex_parameters (int matrix_layout_i, char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i, lapack_int ka_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	range = range_i;
	ka = ka_i;
	kb = ka;
	
	if (jobz == 'U')
		jobz = 'N';

	
	if ((range == 'N') || (range == 'V'))
	{	range = 'A';
		m = n;
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;		
	}
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hbgvx dcomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d \n", matrix_layout, jobz, range, uplo, n, m, ldz);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
		lda = ka+1;
		ldb = kb+1;
		bufsize_a = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		lda = n;
		ldb = n;
		bufsize_a = lda *(ka+1);
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	ldq = n;
	bufsize_q = ldq *n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&bp, &bpref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&q, &qref, bufsize_q);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref,	n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(bp == NULL) || (bpref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(q == NULL) || (qref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)){
		EXPECT_FALSE( true) << "hbgvx_double_parameters object: malloc error.";
		hbgvx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(ap, apref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(bp, bpref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hbgvx_dcomplex_parameters :: ~hbgvx_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbgvx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbgvx_free();
  

}
/*  Test fixture class definition */
class zhbgvx_test  : public  ::testing::Test {
public:
   hbgvx_dcomplex_parameters  *zhbgvx_obj;
   void SetUp();
   void TearDown () { delete zhbgvx_obj; }
};

void zhbgvx_test::SetUp(){

    /* LAPACKE zhbgvx prototype */
    typedef int (*Fptr_NL_LAPACKE_zhbgvx) ( int matrix_layout,  char jobz, char range,	char uplo, lapack_int n, lapack_int ka, 
										lapack_int kb, lapack_complex_double* ap, lapack_int lda, lapack_complex_double* bp, lapack_int ldb, lapack_complex_double* q, 
										lapack_int ldq, double vl, double vu,	lapack_int il, lapack_int iu, double abstol, lapack_int* m, double* w, lapack_complex_double* z,
										lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_zhbgvx zhbgvx;

    zhbgvx_obj = new hbgvx_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m,
							eig_paramslist[idx].k);

    idx = Circular_Increment_Index(idx);

    zhbgvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhbgvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhbgvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhbgvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhbgvx = (Fptr_NL_LAPACKE_zhbgvx)dlsym(zhbgvx_obj->hModule, "LAPACKE_zhbgvx");
    ASSERT_TRUE(zhbgvx != NULL) << "failed to get the Netlib LAPACKE_zhbgvx symbol";
    

    zhbgvx_obj->inforef = zhbgvx( zhbgvx_obj->matrix_layout, zhbgvx_obj->jobz, zhbgvx_obj->range, zhbgvx_obj->uplo, zhbgvx_obj->n, 
	zhbgvx_obj->ka, zhbgvx_obj->kb, zhbgvx_obj->apref, zhbgvx_obj->lda, zhbgvx_obj->bpref, zhbgvx_obj->ldb, zhbgvx_obj->qref, zhbgvx_obj->ldq,  zhbgvx_obj->vl, 
	zhbgvx_obj->vu, zhbgvx_obj->il, zhbgvx_obj->iu,	zhbgvx_obj->abstol, &zhbgvx_obj->m, zhbgvx_obj->wref, zhbgvx_obj->zref, zhbgvx_obj->ldz, zhbgvx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    zhbgvx_obj->info = LAPACKE_zhbgvx( zhbgvx_obj->matrix_layout, zhbgvx_obj->jobz, zhbgvx_obj->range, zhbgvx_obj->uplo,zhbgvx_obj->n, 
	zhbgvx_obj->ka, zhbgvx_obj->kb, zhbgvx_obj->ap, zhbgvx_obj->lda, zhbgvx_obj->bp, zhbgvx_obj->ldb, zhbgvx_obj->q, zhbgvx_obj->ldq, zhbgvx_obj->vl, 
	zhbgvx_obj->vu, zhbgvx_obj->il, zhbgvx_obj->iu,	zhbgvx_obj->abstol, &zhbgvx_obj->m, zhbgvx_obj->w, zhbgvx_obj->zref, zhbgvx_obj->ldz, zhbgvx_obj->ifail);

    if( zhbgvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhbgvx is wrong\n", zhbgvx_obj->info );
    }
    if( zhbgvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhbgvx is wrong\n", 
        zhbgvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhbgvx_obj->diff_a =  computeDiff_z( zhbgvx_obj->bufsize_a, zhbgvx_obj->ap, zhbgvx_obj->apref );
	zhbgvx_obj->diff_b =  computeDiff_z( zhbgvx_obj->bufsize_a, zhbgvx_obj->bp, zhbgvx_obj->bpref );
	zhbgvx_obj->diff_z =  computeDiff_z( zhbgvx_obj->bufsize_z, zhbgvx_obj->z, zhbgvx_obj->zref );
	zhbgvx_obj->diff_w =  computeDiff_d( zhbgvx_obj->n, zhbgvx_obj->w, zhbgvx_obj->wref );
}

TEST_F(zhbgvx_test, zhbgvx1) {
    EXPECT_NEAR(0.0, zhbgvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhbgvx_test, zhbgvx2) {
    EXPECT_NEAR(0.0, zhbgvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhbgvx_test, zhbgvx3) {
    EXPECT_NEAR(0.0, zhbgvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhbgvx_test, zhbgvx4) {
    EXPECT_NEAR(0.0, zhbgvx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}
