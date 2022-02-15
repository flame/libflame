#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hbgvx_2stage_free() \
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
class hbgvx_2stage_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	int bufsize_z;
	int bufsize_q;
	void *hModule, *dModule;
	void *lmodule, *bmodule;
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
      hbgvx_2stage_scomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz, lapack_int ka);
      ~hbgvx_2stage_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hbgvx_2stage_scomplex_parameters:: hbgvx_2stage_scomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i, lapack_int ka_i)
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
	
	//if (jobz == 'U')
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
	printf(" \n hbgvx_2stage scomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
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
		EXPECT_FALSE( true) << "hbgvx_2stage_float_parameters object: malloc error.";
		hbgvx_2stage_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(ap, apref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand(bp, bpref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hbgvx_2stage_scomplex_parameters :: ~hbgvx_2stage_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbgvx_2stage_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbgvx_2stage_free();
  

}
/*  Test fixture class definition */
class chbgvx_2stage_test  : public  ::testing::Test {
public:
   hbgvx_2stage_scomplex_parameters  *chbgvx_2stage_obj;
   void SetUp();
   void TearDown () { delete chbgvx_2stage_obj; }
};

void chbgvx_2stage_test::SetUp(){

    /* LAPACKE chbgvx_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_chbgvx_2stage) ( int matrix_layout, char jobz, char range,	char uplo, lapack_int n, lapack_int ka, 
										lapack_int kb, lapack_complex_float* ap, lapack_int lda, lapack_complex_float* bp, lapack_int ldb, lapack_complex_float* q, 
										lapack_int ldq, float vl, float vu,	lapack_int il, lapack_int iu, float abstol, lapack_int* m, float* w, lapack_complex_float* z,
										lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_chbgvx_2stage chbgvx_2stage;

    chbgvx_2stage_obj = new hbgvx_2stage_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m,
							eig_paramslist[idx].k);

    idx = Circular_Increment_Index(idx);

    chbgvx_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chbgvx_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chbgvx_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chbgvx_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chbgvx_2stage = (Fptr_NL_LAPACKE_chbgvx_2stage)dlsym(chbgvx_2stage_obj->hModule, "LAPACKE_chbgvx_2stage");
    ASSERT_TRUE(chbgvx_2stage != NULL) << "failed to get the Netlib LAPACKE_chbgvx_2stage symbol";
    

    chbgvx_2stage_obj->inforef = chbgvx_2stage( chbgvx_2stage_obj->matrix_layout,  chbgvx_2stage_obj->jobz, chbgvx_2stage_obj->range, chbgvx_2stage_obj->uplo, chbgvx_2stage_obj->n, 
	chbgvx_2stage_obj->ka, chbgvx_2stage_obj->kb, chbgvx_2stage_obj->apref, chbgvx_2stage_obj->lda, chbgvx_2stage_obj->bpref, chbgvx_2stage_obj->ldb, chbgvx_2stage_obj->qref, chbgvx_2stage_obj->ldq, 
	chbgvx_2stage_obj->vl, chbgvx_2stage_obj->vu, chbgvx_2stage_obj->il, chbgvx_2stage_obj->iu,	chbgvx_2stage_obj->abstol, &chbgvx_2stage_obj->m, chbgvx_2stage_obj->w, chbgvx_2stage_obj->zref, 
	chbgvx_2stage_obj->ldz, chbgvx_2stage_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    chbgvx_2stage_obj->info = LAPACKE_chbgvx_2stage( chbgvx_2stage_obj->matrix_layout, chbgvx_2stage_obj->jobz, chbgvx_2stage_obj->range, chbgvx_2stage_obj->uplo, chbgvx_2stage_obj->n, 
	chbgvx_2stage_obj->ka, chbgvx_2stage_obj->kb, chbgvx_2stage_obj->ap, chbgvx_2stage_obj->lda, chbgvx_2stage_obj->bp, chbgvx_2stage_obj->ldb, chbgvx_2stage_obj->q, chbgvx_2stage_obj->ldq, chbgvx_2stage_obj->vl, 
	chbgvx_2stage_obj->vu, chbgvx_2stage_obj->il, chbgvx_2stage_obj->iu,	chbgvx_2stage_obj->abstol, &chbgvx_2stage_obj->m, chbgvx_2stage_obj->w, chbgvx_2stage_obj->z, chbgvx_2stage_obj->ldz, chbgvx_2stage_obj->ifail);

    if( chbgvx_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chbgvx_2stage is wrong\n", chbgvx_2stage_obj->info );
    }
    if( chbgvx_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chbgvx_2stage is wrong\n", 
        chbgvx_2stage_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chbgvx_2stage_obj->diff_a =  computeDiff_c( chbgvx_2stage_obj->bufsize_a, chbgvx_2stage_obj->ap, chbgvx_2stage_obj->apref );
	chbgvx_2stage_obj->diff_b =  computeDiff_c( chbgvx_2stage_obj->bufsize_a, chbgvx_2stage_obj->bp, chbgvx_2stage_obj->bpref );
	chbgvx_2stage_obj->diff_z =  computeDiff_c( chbgvx_2stage_obj->bufsize_z, chbgvx_2stage_obj->z, chbgvx_2stage_obj->zref );
	chbgvx_2stage_obj->diff_w =  computeDiff_s( chbgvx_2stage_obj->n, chbgvx_2stage_obj->w, chbgvx_2stage_obj->wref );
}

TEST_F(chbgvx_2stage_test, chbgvx_2stage1) {
    EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbgvx_2stage_test, chbgvx_2stage2) {
    EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbgvx_2stage_test, chbgvx_2stage3) {
    EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbgvx_2stage_test, chbgvx_2stage4) {
    EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvx_2stage_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hbgvx_2stage_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	int bufsize_z;
	int bufsize_q;
	void *hModule, *dModule;
	void *bmodule, *lmodule;
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
      hbgvx_2stage_dcomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz, lapack_int ka);
      ~hbgvx_2stage_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hbgvx_2stage_dcomplex_parameters:: hbgvx_2stage_dcomplex_parameters (int matrix_layout_i, char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i, lapack_int ka_i)
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
	
	//if (jobz == 'U')
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
	printf(" \n hbgvx_2stage dcomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d \n", matrix_layout, jobz, range, uplo, n, m, ldz);
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
		EXPECT_FALSE( true) << "hbgvx_2stage_double_parameters object: malloc error.";
		hbgvx_2stage_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(ap, apref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(bp, bpref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hbgvx_2stage_dcomplex_parameters :: ~hbgvx_2stage_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbgvx_2stage_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbgvx_2stage_free();
  

}
/*  Test fixture class definition */
class zhbgvx_2stage_test  : public  ::testing::Test {
public:
   hbgvx_2stage_dcomplex_parameters  *zhbgvx_2stage_obj;
   void SetUp();
   void TearDown () { delete zhbgvx_2stage_obj; }
};

void zhbgvx_2stage_test::SetUp(){

    /* LAPACKE zhbgvx_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_zhbgvx_2stage) ( int matrix_layout,  char jobz, char range,	char uplo, lapack_int n, lapack_int ka, 
										lapack_int kb, lapack_complex_double* ap, lapack_int lda, lapack_complex_double* bp, lapack_int ldb, lapack_complex_double* q, 
										lapack_int ldq, double vl, double vu,	lapack_int il, lapack_int iu, double abstol, lapack_int* m, double* w, lapack_complex_double* z,
										lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_zhbgvx_2stage zhbgvx_2stage;

    zhbgvx_2stage_obj = new hbgvx_2stage_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m,
							eig_paramslist[idx].k);

    idx = Circular_Increment_Index(idx);

    zhbgvx_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhbgvx_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhbgvx_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhbgvx_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhbgvx_2stage = (Fptr_NL_LAPACKE_zhbgvx_2stage)dlsym(zhbgvx_2stage_obj->hModule, "LAPACKE_zhbgvx_2stage");
    ASSERT_TRUE(zhbgvx_2stage != NULL) << "failed to get the Netlib LAPACKE_zhbgvx_2stage symbol";
    

    zhbgvx_2stage_obj->inforef = zhbgvx_2stage( zhbgvx_2stage_obj->matrix_layout, zhbgvx_2stage_obj->jobz, zhbgvx_2stage_obj->range, zhbgvx_2stage_obj->uplo, zhbgvx_2stage_obj->n, 
	zhbgvx_2stage_obj->ka, zhbgvx_2stage_obj->kb, zhbgvx_2stage_obj->apref, zhbgvx_2stage_obj->lda, zhbgvx_2stage_obj->bpref, zhbgvx_2stage_obj->ldb, zhbgvx_2stage_obj->qref, zhbgvx_2stage_obj->ldq,  zhbgvx_2stage_obj->vl, 
	zhbgvx_2stage_obj->vu, zhbgvx_2stage_obj->il, zhbgvx_2stage_obj->iu,	zhbgvx_2stage_obj->abstol, &zhbgvx_2stage_obj->m, zhbgvx_2stage_obj->w, zhbgvx_2stage_obj->zref, zhbgvx_2stage_obj->ldz, zhbgvx_2stage_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    zhbgvx_2stage_obj->info = LAPACKE_zhbgvx_2stage( zhbgvx_2stage_obj->matrix_layout, zhbgvx_2stage_obj->jobz, zhbgvx_2stage_obj->range, zhbgvx_2stage_obj->uplo,zhbgvx_2stage_obj->n, 
	zhbgvx_2stage_obj->ka, zhbgvx_2stage_obj->kb, zhbgvx_2stage_obj->ap, zhbgvx_2stage_obj->lda, zhbgvx_2stage_obj->bp, zhbgvx_2stage_obj->ldb, zhbgvx_2stage_obj->q, zhbgvx_2stage_obj->ldq, zhbgvx_2stage_obj->vl, 
	zhbgvx_2stage_obj->vu, zhbgvx_2stage_obj->il, zhbgvx_2stage_obj->iu,	zhbgvx_2stage_obj->abstol, &zhbgvx_2stage_obj->m, zhbgvx_2stage_obj->w, zhbgvx_2stage_obj->zref, zhbgvx_2stage_obj->ldz, zhbgvx_2stage_obj->ifail);

    if( zhbgvx_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhbgvx_2stage is wrong\n", zhbgvx_2stage_obj->info );
    }
    if( zhbgvx_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhbgvx_2stage is wrong\n", 
        zhbgvx_2stage_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhbgvx_2stage_obj->diff_a =  computeDiff_z( zhbgvx_2stage_obj->bufsize_a, zhbgvx_2stage_obj->ap, zhbgvx_2stage_obj->apref );
	zhbgvx_2stage_obj->diff_b =  computeDiff_z( zhbgvx_2stage_obj->bufsize_a, zhbgvx_2stage_obj->bp, zhbgvx_2stage_obj->bpref );
	zhbgvx_2stage_obj->diff_z =  computeDiff_z( zhbgvx_2stage_obj->bufsize_z, zhbgvx_2stage_obj->z, zhbgvx_2stage_obj->zref );
	zhbgvx_2stage_obj->diff_w =  computeDiff_d( zhbgvx_2stage_obj->n, zhbgvx_2stage_obj->w, zhbgvx_2stage_obj->wref );
}

TEST_F(zhbgvx_2stage_test, zhbgvx_2stage1) {
    EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbgvx_2stage_test, zhbgvx_2stage2) {
    EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbgvx_2stage_test, zhbgvx_2stage3) {
    EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbgvx_2stage_test, zhbgvx_2stage4) {
    EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvx_2stage_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}
