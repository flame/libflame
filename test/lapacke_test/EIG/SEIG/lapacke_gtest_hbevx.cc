#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hbevx_free() \
if (ab!=NULL)    free(ab); \
if (abref!=NULL) free(abref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
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
class hbevx_scomplex_parameters{

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
	lapack_int ldz,ldab, ldq;
	lapack_int  kd;
	lapack_int il, iu;
	lapack_int* ifail;
	lapack_complex_float* ab,  *z, *q;
	char uplo;
	char jobz;
	char range;
	/*Output Parameter*/
	float* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_float *abref,  *zref, *qref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hbevx_scomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz, lapack_int kd);
      ~hbevx_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hbevx_scomplex_parameters:: hbevx_scomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i, lapack_int kd_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	range = range_i;
	kd = kd_i;
	
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
	printf(" \n hbevx scomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
		ldab = kd+1;		
		bufsize_a = ldab*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		ldab = n;
		bufsize_a = ldab *(kd+1);
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	ldq = n;
	bufsize_q = ldq *n;
	
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ab, &abref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&q, &qref, bufsize_q);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref,	n);
	
	if ((ab==NULL) || (abref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)||
		(q == NULL) || (qref == NULL)){
		EXPECT_FALSE( true) << "hbevx_float_parameters object: malloc error.";
		hbevx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(ab, abref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hbevx_scomplex_parameters :: ~hbevx_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbevx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbevx_free();
  

}
/*  Test fixture class definition */
class chbevx_test  : public  ::testing::Test {
public:
   hbevx_scomplex_parameters  *chbevx_obj;
   void SetUp();
   void TearDown () { delete chbevx_obj; }
};

void chbevx_test::SetUp(){

    /* LAPACKE chbevx prototype */
    typedef int (*Fptr_NL_LAPACKE_chbevx) ( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, lapack_int kd,
                           lapack_complex_float* ab, lapack_int ldab,
                           lapack_complex_float* q, lapack_int ldq, float vl,
                           float vu, lapack_int il, lapack_int iu, float abstol,
                           lapack_int* m, float* w, lapack_complex_float* z,
                           lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_chbevx chbevx;

    chbevx_obj = new hbevx_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m,
							eig_paramslist[idx].k);

    idx = Circular_Increment_Index(idx);

    chbevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chbevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chbevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chbevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chbevx = (Fptr_NL_LAPACKE_chbevx)dlsym(chbevx_obj->hModule, "LAPACKE_chbevx");
    ASSERT_TRUE(chbevx != NULL) << "failed to get the Netlib LAPACKE_chbevx symbol";
    

    chbevx_obj->inforef = chbevx( chbevx_obj->matrix_layout,  chbevx_obj->jobz, chbevx_obj->range, chbevx_obj->uplo, chbevx_obj->n, 
	chbevx_obj->kd, chbevx_obj->abref, chbevx_obj->ldab, chbevx_obj->qref, chbevx_obj->ldq, chbevx_obj->vl, chbevx_obj->vu, 
	chbevx_obj->il, chbevx_obj->iu,	chbevx_obj->abstol, &chbevx_obj->m, chbevx_obj->wref, chbevx_obj->zref, 
	chbevx_obj->ldz, chbevx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    chbevx_obj->info = LAPACKE_chbevx( chbevx_obj->matrix_layout, chbevx_obj->jobz, chbevx_obj->range, chbevx_obj->uplo, chbevx_obj->n, 
	chbevx_obj->kd,chbevx_obj->ab, chbevx_obj->ldab, chbevx_obj->q, chbevx_obj->ldq, chbevx_obj->vl, chbevx_obj->vu, 
	chbevx_obj->il, chbevx_obj->iu,	chbevx_obj->abstol, &chbevx_obj->m, chbevx_obj->w, chbevx_obj->z, chbevx_obj->ldz, chbevx_obj->ifail);

    if( chbevx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chbevx is wrong\n", chbevx_obj->info );
    }
    if( chbevx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chbevx is wrong\n", 
        chbevx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chbevx_obj->diff_a =  computeDiff_c( chbevx_obj->bufsize_a, chbevx_obj->ab, chbevx_obj->abref );
	chbevx_obj->diff_z =  computeDiff_c( chbevx_obj->bufsize_z, chbevx_obj->z, chbevx_obj->zref );
	chbevx_obj->diff_w =  computeDiff_s( chbevx_obj->n, chbevx_obj->w, chbevx_obj->wref );
}

TEST_F(chbevx_test, chbevx1) {
    EXPECT_NEAR(0.0, chbevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chbevx_test, chbevx2) {
    EXPECT_NEAR(0.0, chbevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chbevx_test, chbevx3) {
    EXPECT_NEAR(0.0, chbevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chbevx_test, chbevx4) {
    EXPECT_NEAR(0.0, chbevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chbevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hbevx_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	int bufsize_z;
	int bufsize_q;
	void *hModule, *dModule;
	void *lmodule, *bmodule;
	double diff_a;
	double diff_w;
	double diff_z;
	double vl, vu;
	double abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz,ldab, ldq;
	lapack_int  kd;
	lapack_int il, iu;
	lapack_int* ifail;
	lapack_complex_double* ab,  *z, *q;
	char uplo;
	char jobz;
	char range;
	/*Output Parameter*/
	double* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_double *abref,  *zref, *qref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hbevx_dcomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz, lapack_int kd);
      ~hbevx_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hbevx_dcomplex_parameters:: hbevx_dcomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i, lapack_int kd_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	range = range_i;
	kd = kd_i;
	
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
	printf(" \n hbevx dcomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
		ldab = kd+1;		
		bufsize_a = ldab*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		ldab = n;
		bufsize_a = ldab *(kd+1);
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	ldq = n;
	bufsize_q = ldq *n;
	
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ab, &abref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&q, &qref, bufsize_q);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref,	n);
	
	if ((ab==NULL) || (abref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)||
		(q == NULL) || (qref == NULL)){
		EXPECT_FALSE( true) << "hbevx_double_parameters object: malloc error.";
		hbevx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(ab, abref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hbevx_dcomplex_parameters :: ~hbevx_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbevx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbevx_free();
  

}
/*  Test fixture class definition */
class zhbevx_test  : public  ::testing::Test {
public:
   hbevx_dcomplex_parameters  *zhbevx_obj;
   void SetUp();
   void TearDown () { delete zhbevx_obj; }
};

void zhbevx_test::SetUp(){

    /* LAPACKE zhbevx prototype */
    typedef int (*Fptr_NL_LAPACKE_zhbevx) ( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, lapack_int kd,
                           lapack_complex_double* ab, lapack_int ldab,
                           lapack_complex_double* q, lapack_int ldq, double vl,
                           double vu, lapack_int il, lapack_int iu, double abstol,
                           lapack_int* m, double* w, lapack_complex_double* z,
                           lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_zhbevx zhbevx;

    zhbevx_obj = new hbevx_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m,
							eig_paramslist[idx].k);

    idx = Circular_Increment_Index(idx);

    zhbevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhbevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhbevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhbevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhbevx = (Fptr_NL_LAPACKE_zhbevx)dlsym(zhbevx_obj->hModule, "LAPACKE_zhbevx");
    ASSERT_TRUE(zhbevx != NULL) << "failed to get the Netlib LAPACKE_zhbevx symbol";
    

    zhbevx_obj->inforef = zhbevx( zhbevx_obj->matrix_layout,  zhbevx_obj->jobz, zhbevx_obj->range, zhbevx_obj->uplo, zhbevx_obj->n, 
	zhbevx_obj->kd, zhbevx_obj->abref, zhbevx_obj->ldab, zhbevx_obj->qref, zhbevx_obj->ldq, zhbevx_obj->vl, zhbevx_obj->vu, 
	zhbevx_obj->il, zhbevx_obj->iu,	zhbevx_obj->abstol, &zhbevx_obj->m, zhbevx_obj->wref, zhbevx_obj->zref, 
	zhbevx_obj->ldz, zhbevx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    zhbevx_obj->info = LAPACKE_zhbevx( zhbevx_obj->matrix_layout, zhbevx_obj->jobz, zhbevx_obj->range, zhbevx_obj->uplo, zhbevx_obj->n, 
	zhbevx_obj->kd,zhbevx_obj->ab, zhbevx_obj->ldab, zhbevx_obj->q, zhbevx_obj->ldq, zhbevx_obj->vl, zhbevx_obj->vu, 
	zhbevx_obj->il, zhbevx_obj->iu,	zhbevx_obj->abstol, &zhbevx_obj->m, zhbevx_obj->w, zhbevx_obj->z, zhbevx_obj->ldz, zhbevx_obj->ifail);

    if( zhbevx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhbevx is wrong\n", zhbevx_obj->info );
    }
    if( zhbevx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhbevx is wrong\n", 
        zhbevx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhbevx_obj->diff_a =  computeDiff_z( zhbevx_obj->bufsize_a, zhbevx_obj->ab, zhbevx_obj->abref );
	zhbevx_obj->diff_z =  computeDiff_z( zhbevx_obj->bufsize_z, zhbevx_obj->z, zhbevx_obj->zref );
	zhbevx_obj->diff_w =  computeDiff_d( zhbevx_obj->n, zhbevx_obj->w, zhbevx_obj->wref );
}

TEST_F(zhbevx_test, zhbevx1) {
    EXPECT_NEAR(0.0, zhbevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhbevx_test, zhbevx2) {
    EXPECT_NEAR(0.0, zhbevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhbevx_test, zhbevx3) {
    EXPECT_NEAR(0.0, zhbevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhbevx_test, zhbevx4) {
    EXPECT_NEAR(0.0, zhbevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhbevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}