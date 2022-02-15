#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hpevx_free() \
if (ap!=NULL)    free(ap); \
if (apref!=NULL) free(apref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (z!=NULL)  free(z);\
if (zref!=NULL)  free(zref);\
if (ifail!=NULL)  free(ifail);\
if (ifailref!=NULL)  free(ifailref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hpevx_scomplex_parameters{

   public:
	int bufsize_a;
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
	lapack_complex_float* ap,  *z;
	char uplo;
	char jobz;
	char range;
	/*Output Parameter*/
	float* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_float *apref,  *zref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hpevx_scomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~hpevx_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hpevx_scomplex_parameters:: hpevx_scomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	range = range_i;
	
	
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
	printf(" \n hpevx scomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
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
	
	bufsize_a = n*(n+1)/2;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref,	n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)){
		EXPECT_FALSE( true) << "hpevx_float_parameters object: malloc error.";
		hpevx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(ap, apref, bufsize_a, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(ap, apref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hpevx_scomplex_parameters :: ~hpevx_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpevx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpevx_free();
  

}
/*  Test fixture class definition */
class chpevx_test  : public  ::testing::Test {
public:
   hpevx_scomplex_parameters  *chpevx_obj;
   void SetUp();
   void TearDown () { delete chpevx_obj; }
};

void chpevx_test::SetUp(){

    /* LAPACKE chpevx prototype */
    typedef int (*Fptr_NL_LAPACKE_chpevx) ( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, lapack_complex_float* ap, float vl,
                           float vu, lapack_int il, lapack_int iu, float abstol,
                           lapack_int* m, float* w, lapack_complex_float* z,
                           lapack_int ldz, lapack_int* ifail);

    Fptr_NL_LAPACKE_chpevx chpevx;

    chpevx_obj = new hpevx_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    chpevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chpevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chpevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chpevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chpevx = (Fptr_NL_LAPACKE_chpevx)dlsym(chpevx_obj->hModule, "LAPACKE_chpevx");
    ASSERT_TRUE(chpevx != NULL) << "failed to get the Netlib LAPACKE_chpevx symbol";
    

    chpevx_obj->inforef = chpevx( chpevx_obj->matrix_layout,  chpevx_obj->jobz, chpevx_obj->range, chpevx_obj->uplo, chpevx_obj->n, 
	chpevx_obj->apref, chpevx_obj->vl, chpevx_obj->vu, chpevx_obj->il, chpevx_obj->iu,	chpevx_obj->abstol, &chpevx_obj->m,
	chpevx_obj->wref, chpevx_obj->zref, chpevx_obj->ldz, chpevx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    chpevx_obj->info = LAPACKE_chpevx( chpevx_obj->matrix_layout, chpevx_obj->jobz, chpevx_obj->range, chpevx_obj->uplo, chpevx_obj->n, 
	chpevx_obj->ap, chpevx_obj->vl, chpevx_obj->vu, chpevx_obj->il, chpevx_obj->iu,	chpevx_obj->abstol, &chpevx_obj->m,
	chpevx_obj->w, chpevx_obj->z, chpevx_obj->ldz, chpevx_obj->ifail);

    if( chpevx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chpevx is wrong\n", chpevx_obj->info );
    }
    if( chpevx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chpevx is wrong\n", 
        chpevx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chpevx_obj->diff_a =  computeDiff_c( chpevx_obj->bufsize_a, chpevx_obj->ap, chpevx_obj->apref );
	chpevx_obj->diff_z =  computeDiff_c( chpevx_obj->bufsize_z, chpevx_obj->z, chpevx_obj->zref );
	chpevx_obj->diff_w =  computeDiff_s( chpevx_obj->n, chpevx_obj->w, chpevx_obj->wref );
}

TEST_F(chpevx_test, chpevx1) {
    EXPECT_NEAR(0.0, chpevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chpevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chpevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chpevx_test, chpevx2) {
    EXPECT_NEAR(0.0, chpevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chpevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chpevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chpevx_test, chpevx3) {
    EXPECT_NEAR(0.0, chpevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chpevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chpevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(chpevx_test, chpevx4) {
    EXPECT_NEAR(0.0, chpevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chpevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, chpevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hpevx_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_z;
	void *hModule, *dModule;
	double diff_a;
	double diff_w;
	double diff_z;
	double vl, vu;
	double abstol;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldz;
	lapack_int il, iu;
	lapack_int* ifail;
	lapack_complex_double* ap,  *z;
	char uplo;
	char jobz;
	char range;
	/*Output Parameter*/
	double* w;
	lapack_int m;
	lapack_int* ifailref;
	lapack_complex_double *apref,  *zref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hpevx_dcomplex_parameters (int matrix_layout, char jobz, char range,  char uplo, lapack_int n, lapack_int m, lapack_int ldz);
      ~hpevx_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hpevx_dcomplex_parameters:: hpevx_dcomplex_parameters (int matrix_layout_i,  char jobz_i, char range_i, char uplo_i, lapack_int n_i, lapack_int m_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	range = range_i;	
	
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
	printf(" \n hpevx dcomplex: matrix_layout = %d,  jobz:%c, range:%c, uplo:%c, n: %d , m:%d, ldz: %d, iu:%d, Il:%d \n", matrix_layout, jobz, range, uplo, n, m, ldz, iu, il);
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
	
	bufsize_a = n*(n+1)/2;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ap, &apref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	lapacke_gtest_alloc_int_buffer_pair(&ifail, &ifailref,	n);
	
	if ((ap==NULL) || (apref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(ifail == NULL) || (ifailref == NULL)){
		EXPECT_FALSE( true) << "hpevx_double_parameters object: malloc error.";
		hpevx_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(ap, apref, bufsize_a, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(ap, apref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hpevx_dcomplex_parameters :: ~hpevx_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpevx_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpevx_free();
  

}
/*  Test fixture class definition */
class zhpevx_test  : public  ::testing::Test {
public:
   hpevx_dcomplex_parameters  *zhpevx_obj;
   void SetUp();
   void TearDown () { delete zhpevx_obj; }
};

void zhpevx_test::SetUp(){

    /* LAPACKE zhpevx prototype */
    typedef int (*Fptr_NL_LAPACKE_zhpevx) (int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, lapack_complex_double* ap, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w,
                           lapack_complex_double* z, lapack_int ldz,
                           lapack_int* ifail);

    Fptr_NL_LAPACKE_zhpevx zhpevx;

    zhpevx_obj = new hpevx_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].compz,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].m,
							eig_paramslist[idx].n,
							eig_paramslist[idx].m);

    idx = Circular_Increment_Index(idx);

    zhpevx_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhpevx_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhpevx_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhpevx_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhpevx = (Fptr_NL_LAPACKE_zhpevx)dlsym(zhpevx_obj->hModule, "LAPACKE_zhpevx");
    ASSERT_TRUE(zhpevx != NULL) << "failed to get the Netlib LAPACKE_zhpevx symbol";
    

    zhpevx_obj->inforef = zhpevx( zhpevx_obj->matrix_layout,  zhpevx_obj->jobz, zhpevx_obj->range, zhpevx_obj->uplo, zhpevx_obj->n, 
	zhpevx_obj->apref, zhpevx_obj->vl, zhpevx_obj->vu, zhpevx_obj->il, zhpevx_obj->iu,	zhpevx_obj->abstol, &zhpevx_obj->m, 
	zhpevx_obj->wref, zhpevx_obj->zref, zhpevx_obj->ldz, zhpevx_obj->ifailref);

    /* Compute libflame's Lapacke o/p  */
    zhpevx_obj->info = LAPACKE_zhpevx( zhpevx_obj->matrix_layout, zhpevx_obj->jobz, zhpevx_obj->range, zhpevx_obj->uplo, zhpevx_obj->n, 
	zhpevx_obj->ap, zhpevx_obj->vl, zhpevx_obj->vu, zhpevx_obj->il, zhpevx_obj->iu,	zhpevx_obj->abstol, &zhpevx_obj->m, 
	zhpevx_obj->w, zhpevx_obj->z, zhpevx_obj->ldz, zhpevx_obj->ifail);

    if( zhpevx_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhpevx is wrong\n", zhpevx_obj->info );
    }
    if( zhpevx_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhpevx is wrong\n", 
        zhpevx_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhpevx_obj->diff_a =  computeDiff_z( zhpevx_obj->bufsize_a, zhpevx_obj->ap, zhpevx_obj->apref );
	zhpevx_obj->diff_z =  computeDiff_z( zhpevx_obj->bufsize_z, zhpevx_obj->z, zhpevx_obj->zref );
	zhpevx_obj->diff_w =  computeDiff_d( zhpevx_obj->n, zhpevx_obj->w, zhpevx_obj->wref );
}

TEST_F(zhpevx_test, zhpevx1) {
    EXPECT_NEAR(0.0, zhpevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhpevx_test, zhpevx2) {
    EXPECT_NEAR(0.0, zhpevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhpevx_test, zhpevx3) {
    EXPECT_NEAR(0.0, zhpevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zhpevx_test, zhpevx4) {
    EXPECT_NEAR(0.0, zhpevx_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevx_obj->diff_w, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevx_obj->diff_z, LAPACKE_EIG_THRESHOLD);
}