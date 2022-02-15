#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hpgvx_free() \
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (bp!=NULL)  free(bp);\
if (bpref!=NULL)  free(bpref);\
if (z!=NULL)  free(z);\
if (zref!=NULL)  free(zref);\
if (ifail!=NULL)  free(ifail);\
if (ifailref!=NULL)  free(ifailref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hpgvx_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	int bufsize_z;
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
	lapack_int ldz;	
	lapack_int il, iu;
	lapack_int* ifail;
	lapack_complex_float* ap, *bp, *z;
	char uplo;
	char jobz;
	char range;
	lapack_int itype;
	/*Output Parameter*/
	float* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_float *apref, *bpref, *zref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hpgvx_scomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~hpgvx_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hpgvx_scomplex_parameters:: hpgvx_scomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	itype = itype_i;
	range = range_i;	
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (itype == 4)
		itype = 1;
	
	if (range == 'N')
	{	range = 'A';
		m = n;
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;		
	}
	
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hpgvx scomplex: matrix_layout = %d, itype:%d, jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d \n", matrix_layout, itype, jobz, range, uplo, n, m, ldz);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_a =  n*(n+1)/2;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&bp, &bpref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref,	n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(bp == NULL) || (bpref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)){
		EXPECT_FALSE( true) << "hpgvx_float_parameters object: malloc error.";
		hpgvx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(ap, apref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand(bp, bpref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hpgvx_scomplex_parameters :: ~hpgvx_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpgvx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpgvx_free();
  

}
/*  Test fixture class definition */
class chpgvx_test  : public  ::testing::Test {
public:
   hpgvx_scomplex_parameters  *chpgvx_obj;
   void SetUp();
   void TearDown () { delete chpgvx_obj; }
};

void chpgvx_test::SetUp(){

    /* LAPACKE chpgvx prototype */
    typedef int (*Fptr_NL_LAPACKE_chpgvx) ( int matrix_layout, lapack_int itype, char jobz, char range,
										char uplo, lapack_int n, lapack_complex_float* ap, lapack_complex_float* bp, float vl, float vu,
										lapack_int il, lapack_int iu, float abstol, lapack_int* m, float* w, lapack_complex_float* z,
										lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_chpgvx chpgvx;

    chpgvx_obj = new hpgvx_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].nrhs,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    chpgvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chpgvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chpgvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chpgvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chpgvx = (Fptr_NL_LAPACKE_chpgvx)dlsym(chpgvx_obj->hModule, "LAPACKE_chpgvx");
    ASSERT_TRUE(chpgvx != NULL) << "failed to get the Netlib LAPACKE_chpgvx symbol";
    

    chpgvx_obj->inforef = chpgvx( chpgvx_obj->matrix_layout, chpgvx_obj->itype, chpgvx_obj->jobz, chpgvx_obj->range, chpgvx_obj->uplo, 
								chpgvx_obj->n, chpgvx_obj->apref, chpgvx_obj->bpref, chpgvx_obj->vl, chpgvx_obj->vu, chpgvx_obj->il, chpgvx_obj->iu,
								chpgvx_obj->abstol, &chpgvx_obj->m, chpgvx_obj->w, chpgvx_obj->zref, chpgvx_obj->ldz, chpgvx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    chpgvx_obj->info = LAPACKE_chpgvx( chpgvx_obj->matrix_layout, chpgvx_obj->itype, chpgvx_obj->jobz, chpgvx_obj->range, chpgvx_obj->uplo, 
										chpgvx_obj->n, chpgvx_obj->apref, chpgvx_obj->bpref, chpgvx_obj->vl, chpgvx_obj->vu, chpgvx_obj->il, chpgvx_obj->iu,
										chpgvx_obj->abstol, &chpgvx_obj->m, chpgvx_obj->w, chpgvx_obj->zref, chpgvx_obj->ldz, chpgvx_obj->ifailref);

    if( chpgvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chpgvx is wrong\n", chpgvx_obj->info );
    }
    if( chpgvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chpgvx is wrong\n", 
        chpgvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chpgvx_obj->diff_a =  computeDiff_c( chpgvx_obj->bufsize_a, chpgvx_obj->ap, chpgvx_obj->apref );
	chpgvx_obj->diff_b =  computeDiff_c( chpgvx_obj->bufsize_a, chpgvx_obj->bp, chpgvx_obj->bpref );
	chpgvx_obj->diff_z =  computeDiff_c( chpgvx_obj->bufsize_z, chpgvx_obj->z, chpgvx_obj->zref );
	chpgvx_obj->diff_w =  computeDiff_s( chpgvx_obj->n, chpgvx_obj->w, chpgvx_obj->wref );
}

TEST_F(chpgvx_test, chpgvx1) {
    EXPECT_NEAR(0.0, chpgvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgvx_test, chpgvx2) {
    EXPECT_NEAR(0.0, chpgvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgvx_test, chpgvx3) {
    EXPECT_NEAR(0.0, chpgvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgvx_test, chpgvx4) {
    EXPECT_NEAR(0.0, chpgvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvx_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hpgvx_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	int bufsize_z;
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
	lapack_int ldz;
	lapack_int il, iu;
	lapack_int* ifail;
	lapack_complex_double* ap, *bp, *z;
	char uplo;
	char jobz;
	char range;
	lapack_int itype;
	/*Output Parameter*/
	double* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_double *apref, *bpref, *zref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hpgvx_dcomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~hpgvx_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hpgvx_dcomplex_parameters:: hpgvx_dcomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	itype = itype_i;
	range = range_i;	
	
	if (jobz == 'U')
		jobz = 'N';
	
	if (itype == 4)
		itype = 1;
	
	if (range == 'N')
	{	range = 'A';
		m = n;
	}else if (range == 'I')
	{
		iu = n/2;
		il = n/5;
		m = iu-il+1;		
	}
	
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hpgvx dcomplex: matrix_layout = %d, itype:%d, jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d \n", matrix_layout, itype, jobz, range, uplo, n, m, ldz);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		ldz = n;
		bufsize_z = ldz*m;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
    {
		ldz = m;
		bufsize_z = ldz*n;
    }else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_a =  n*(n+1)/2;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&bp, &bpref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref,	n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(bp == NULL) || (bpref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)){
		EXPECT_FALSE( true) << "hpgvx_double_parameters object: malloc error.";
		hpgvx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(ap, apref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(bp, bpref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hpgvx_dcomplex_parameters :: ~hpgvx_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpgvx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpgvx_free();
  

}
/*  Test fixture class definition */
class zhpgvx_test  : public  ::testing::Test {
public:
   hpgvx_dcomplex_parameters  *zhpgvx_obj;
   void SetUp();
   void TearDown () { delete zhpgvx_obj; }
};

void zhpgvx_test::SetUp(){

    /* LAPACKE zhpgvx prototype */
    typedef int (*Fptr_NL_LAPACKE_zhpgvx) ( int matrix_layout, lapack_int itype, char jobz, char range,
										char uplo, lapack_int n, lapack_complex_double* ap, lapack_complex_double* bp, double vl, double vu,
										lapack_int il, lapack_int iu, double abstol, lapack_int* m, double* w, lapack_complex_double* z,
										lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_zhpgvx zhpgvx;

    zhpgvx_obj = new hpgvx_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].nrhs,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    zhpgvx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhpgvx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhpgvx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhpgvx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhpgvx = (Fptr_NL_LAPACKE_zhpgvx)dlsym(zhpgvx_obj->hModule, "LAPACKE_zhpgvx");
    ASSERT_TRUE(zhpgvx != NULL) << "failed to get the Netlib LAPACKE_zhpgvx symbol";
    

    zhpgvx_obj->inforef = zhpgvx( zhpgvx_obj->matrix_layout, zhpgvx_obj->itype, zhpgvx_obj->jobz, zhpgvx_obj->range, zhpgvx_obj->uplo, 
								zhpgvx_obj->n, zhpgvx_obj->apref, zhpgvx_obj->bpref, zhpgvx_obj->vl, zhpgvx_obj->vu, zhpgvx_obj->il, zhpgvx_obj->iu,
								zhpgvx_obj->abstol, &zhpgvx_obj->m, zhpgvx_obj->w, zhpgvx_obj->zref, zhpgvx_obj->ldz, zhpgvx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    zhpgvx_obj->info = LAPACKE_zhpgvx( zhpgvx_obj->matrix_layout, zhpgvx_obj->itype, zhpgvx_obj->jobz, zhpgvx_obj->range, zhpgvx_obj->uplo, 
										zhpgvx_obj->n, zhpgvx_obj->apref, zhpgvx_obj->bpref, zhpgvx_obj->vl, zhpgvx_obj->vu, zhpgvx_obj->il, zhpgvx_obj->iu,
										zhpgvx_obj->abstol, &zhpgvx_obj->m, zhpgvx_obj->w, zhpgvx_obj->zref, zhpgvx_obj->ldz, zhpgvx_obj->ifailref);

    if( zhpgvx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhpgvx is wrong\n", zhpgvx_obj->info );
    }
    if( zhpgvx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhpgvx is wrong\n", 
        zhpgvx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhpgvx_obj->diff_a =  computeDiff_z( zhpgvx_obj->bufsize_a, zhpgvx_obj->ap, zhpgvx_obj->apref );
	zhpgvx_obj->diff_b =  computeDiff_z( zhpgvx_obj->bufsize_a, zhpgvx_obj->bp, zhpgvx_obj->bpref );
	zhpgvx_obj->diff_z =  computeDiff_z( zhpgvx_obj->bufsize_z, zhpgvx_obj->z, zhpgvx_obj->zref );
	zhpgvx_obj->diff_w =  computeDiff_d( zhpgvx_obj->n, zhpgvx_obj->w, zhpgvx_obj->wref );
}

TEST_F(zhpgvx_test, zhpgvx1) {
    EXPECT_NEAR(0.0, zhpgvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgvx_test, zhpgvx2) {
    EXPECT_NEAR(0.0, zhpgvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgvx_test, zhpgvx3) {
    EXPECT_NEAR(0.0, zhpgvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgvx_test, zhpgvx4) {
    EXPECT_NEAR(0.0, zhpgvx_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvx_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}
