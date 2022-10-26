#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define hpgv_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (b!=NULL)  free(b);\
if (bref!=NULL)  free(bref);\
if (z!=NULL)  free(z);\
if (zref!=NULL)  free(zref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hpgv_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	int bufsize_z;
	void *hModule, *dModule;
	float diff_a;
	float diff_w;
	float diff_b;
	float diff_z;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* A, *b, *z;
	char uplo;
	char jobz;
	lapack_int ldz;
	lapack_int itype;
	/*Output Parameter*/
	float* w;	
	lapack_complex_float *Aref, *bref, *zref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hpgv_scomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n, lapack_int ldz);
      ~hpgv_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hpgv_scomplex_parameters:: hpgv_scomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	itype = itype_i;
	
	if (jobz == 'U')
		jobz = 'N';
	if (itype == 4)
		itype = 1;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hpgv scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
	#endif

	bufsize_a = n*(n+1)/2;
	bufsize_b = bufsize_a;
	bufsize_z = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || 
		(w == NULL) || (wref == NULL) ||
		(b == NULL) || (bref == NULL) ||
		(z == NULL) || (zref == NULL)){
		EXPECT_FALSE( true) << "hpgv_float_parameters object: malloc error.";
		hpgv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand(b, bref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hpgv_scomplex_parameters :: ~hpgv_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpgv_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpgv_free();
  

}
/*  Test fixture class definition */
class chpgv_test  : public  ::testing::Test {
public:
   hpgv_scomplex_parameters  *chpgv_obj;
   void SetUp();
   void TearDown () { delete chpgv_obj; }
};

void chpgv_test::SetUp(){

    /* LAPACKE chpgv prototype */
    typedef int (*Fptr_NL_LAPACKE_chpgv) (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,lapack_complex_float *A, 
										lapack_complex_float* b, float* w, lapack_complex_float* z, lapack_int ldz);

    Fptr_NL_LAPACKE_chpgv chpgv;

    chpgv_obj = new hpgv_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].nrhs,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].n,
							eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chpgv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chpgv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chpgv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chpgv_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chpgv = (Fptr_NL_LAPACKE_chpgv)dlsym(chpgv_obj->hModule, "LAPACKE_chpgv");
    ASSERT_TRUE(chpgv != NULL) << "failed to get the Netlib LAPACKE_chpgv symbol";
    

    chpgv_obj->inforef = chpgv( chpgv_obj->matrix_layout, chpgv_obj->itype, chpgv_obj->jobz, chpgv_obj->uplo, 
								chpgv_obj->n, chpgv_obj->Aref, chpgv_obj->bref, chpgv_obj->w, chpgv_obj->zref, chpgv_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    chpgv_obj->info = LAPACKE_chpgv( chpgv_obj->matrix_layout, chpgv_obj->itype, chpgv_obj->jobz, chpgv_obj->uplo,
										chpgv_obj->n,chpgv_obj->A, chpgv_obj->b, chpgv_obj->w, chpgv_obj->z, chpgv_obj->ldz);

    if( chpgv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chpgv is wrong\n", chpgv_obj->info );
    }
    if( chpgv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chpgv is wrong\n", 
        chpgv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chpgv_obj->diff_a =  computeDiff_c( chpgv_obj->bufsize_a, chpgv_obj->A, chpgv_obj->Aref );
	chpgv_obj->diff_b =  computeDiff_c( chpgv_obj->bufsize_a, chpgv_obj->b, chpgv_obj->bref );
	chpgv_obj->diff_z =  computeDiff_c( chpgv_obj->bufsize_z, chpgv_obj->z, chpgv_obj->zref );
	chpgv_obj->diff_w =  computeDiff_s( chpgv_obj->n, chpgv_obj->w, chpgv_obj->wref );
}

TEST_F(chpgv_test, chpgv1) {
    EXPECT_NEAR(0.0, chpgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgv_test, chpgv2) {
    EXPECT_NEAR(0.0, chpgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgv_test, chpgv3) {
    EXPECT_NEAR(0.0, chpgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgv_test, chpgv4) {
    EXPECT_NEAR(0.0, chpgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgv_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hpgv_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_b;
	int bufsize_z;
	void *hModule, *dModule;
	 double diff_a;
	double diff_w;
	double diff_b;
	double diff_z;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* A, *b, *z;
	lapack_int itype;
	char uplo;
	char jobz;
	lapack_int ldz;
	/*Output Parameter*/
	double* w;	
	lapack_complex_double *Aref, *bref,*zref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hpgv_dcomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n, lapack_int ldz);
      ~hpgv_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hpgv_dcomplex_parameters:: hpgv_dcomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	itype = itype_i;
	
	if (jobz == 'U')
		jobz = 'N';
	if (itype == 4)
		itype = 1;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hpgv scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
	#endif
	
	bufsize_a = n*(n+1)/2;
	bufsize_b = bufsize_a;
	bufsize_z = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(b == NULL) || (bref == NULL) ||
		(z == NULL) || (zref == NULL)){
		EXPECT_FALSE( true) << "hpgv_double_parameters object: malloc error.";
		hpgv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(b, bref, bufsize_a);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hpgv_dcomplex_parameters :: ~hpgv_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpgv_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpgv_free();

}
/*  Test fixture class definition */ 
class zhpgv_test  : public  ::testing::Test {
public:
   hpgv_dcomplex_parameters  *zhpgv_obj;
   void SetUp();
   void TearDown () { delete zhpgv_obj; }
};

void zhpgv_test::SetUp(){

    /* LAPACKE zhpgv prototype */
    typedef int (*Fptr_NL_LAPACKE_zhpgv) (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,lapack_complex_double *A, 
										lapack_complex_double* b, double* w, lapack_complex_double* z, lapack_int ldz );

    Fptr_NL_LAPACKE_zhpgv zhpgv;

    zhpgv_obj = new hpgv_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhpgv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhpgv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhpgv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhpgv_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhpgv = (Fptr_NL_LAPACKE_zhpgv)dlsym(zhpgv_obj->hModule, "LAPACKE_zhpgv");
    ASSERT_TRUE(zhpgv != NULL) << "failed to get the Netlib LAPACKE_zhpgv symbol";
    

    zhpgv_obj->inforef = zhpgv( zhpgv_obj->matrix_layout, zhpgv_obj->itype, zhpgv_obj->jobz, zhpgv_obj->uplo,
								zhpgv_obj->n, zhpgv_obj->Aref, zhpgv_obj->bref, zhpgv_obj->wref, zhpgv_obj->zref, zhpgv_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zhpgv_obj->info = LAPACKE_zhpgv( zhpgv_obj->matrix_layout, zhpgv_obj->itype, zhpgv_obj->jobz, zhpgv_obj->uplo,
										zhpgv_obj->n,zhpgv_obj->A, zhpgv_obj->b, zhpgv_obj->w, zhpgv_obj->z, zhpgv_obj->ldz);

    if( zhpgv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhpgv is wrong\n", zhpgv_obj->info );
    }
    if( zhpgv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhpgv is wrong\n", 
        zhpgv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhpgv_obj->diff_a =  computeDiff_z( zhpgv_obj->bufsize_a, zhpgv_obj->A, zhpgv_obj->Aref );
	zhpgv_obj->diff_b =  computeDiff_z( zhpgv_obj->bufsize_a, zhpgv_obj->b, zhpgv_obj->bref );
	zhpgv_obj->diff_z =  computeDiff_z( zhpgv_obj->bufsize_z, zhpgv_obj->z, zhpgv_obj->zref );
	zhpgv_obj->diff_w =  computeDiff_d( zhpgv_obj->n, zhpgv_obj->w, zhpgv_obj->wref );
}

TEST_F(zhpgv_test, zhpgv1) {
    EXPECT_NEAR(0.0, zhpgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgv_test, zhpgv2) {
    EXPECT_NEAR(0.0, zhpgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgv_test, zhpgv3) {
    EXPECT_NEAR(0.0, zhpgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgv_test, zhpgv4) {
    EXPECT_NEAR(0.0, zhpgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgv_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}
