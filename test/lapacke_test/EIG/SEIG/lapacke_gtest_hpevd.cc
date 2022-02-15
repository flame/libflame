#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hpevd_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (z!=NULL)  free(z);\
if (zref!=NULL)  free(zref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hpevd_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_z;
	void *hModule, *dModule;
	float diff_a;
	float diff_w;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* A, *z;	
	char uplo;
	char jobz;
	lapack_int ldz;
	/*Output Parameter*/
	float* w;	
	lapack_complex_float *Aref, *zref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hpevd_scomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ldz);
      ~hpevd_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hpevd_scomplex_parameters:: hpevd_scomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hpevd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
	#endif

	bufsize_z = ldz*n;
	bufsize_a = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL)){
		EXPECT_FALSE( true) << "hpevd_float_parameters object: malloc error.";
		hpevd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, bufsize_a, bufsize_a, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hpevd_scomplex_parameters :: ~hpevd_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpevd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpevd_free();
  

}
/*  Test fixture class definition */
class chpevd_test  : public  ::testing::Test {
public:
   hpevd_scomplex_parameters  *chpevd_obj;
   void SetUp();
   void TearDown () { delete chpevd_obj; }
};

void chpevd_test::SetUp(){

    /* LAPACKE chpevd prototype */
    typedef int (*Fptr_NL_LAPACKE_chpevd) (int matrix_layout, char jobz, char uplo, lapack_int n,lapack_complex_float *A, 
										float* w, lapack_complex_float* z, lapack_int ldz);

    Fptr_NL_LAPACKE_chpevd chpevd;

    chpevd_obj = new hpevd_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chpevd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chpevd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chpevd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chpevd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chpevd = (Fptr_NL_LAPACKE_chpevd)dlsym(chpevd_obj->hModule, "LAPACKE_chpevd");
    ASSERT_TRUE(chpevd != NULL) << "failed to get the Netlib LAPACKE_chpevd symbol";
    

    chpevd_obj->inforef = chpevd( chpevd_obj->matrix_layout, chpevd_obj->jobz, chpevd_obj->uplo,
								chpevd_obj->n, chpevd_obj->Aref, chpevd_obj->wref, chpevd_obj->zref, chpevd_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    chpevd_obj->info = LAPACKE_chpevd( chpevd_obj->matrix_layout, chpevd_obj->jobz, chpevd_obj->uplo,
										chpevd_obj->n,chpevd_obj->A, chpevd_obj->w, chpevd_obj->z, chpevd_obj->ldz);

    if( chpevd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chpevd is wrong\n", chpevd_obj->info );
    }
    if( chpevd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chpevd is wrong\n", 
        chpevd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chpevd_obj->diff_a =  computeDiff_c( chpevd_obj->bufsize_a, chpevd_obj->A, chpevd_obj->Aref );
	chpevd_obj->diff_w =  computeDiff_s( chpevd_obj->n, chpevd_obj->w, chpevd_obj->wref );
}

TEST_F(chpevd_test, chpevd1) {
    EXPECT_NEAR(0.0, chpevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpevd_test, chpevd2) {
    EXPECT_NEAR(0.0, chpevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpevd_test, chpevd3) {
    EXPECT_NEAR(0.0, chpevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpevd_test, chpevd4) {
    EXPECT_NEAR(0.0, chpevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hpevd_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_z;
	void *hModule, *dModule;
	 double diff_a;
	double diff_w;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* A, *z;
	char uplo;
	char jobz;
	lapack_int ldz;
	/*Output Parameter*/
	double* w;	
	lapack_complex_double *Aref, *zref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hpevd_dcomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ldz);
      ~hpevd_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hpevd_dcomplex_parameters:: hpevd_dcomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hpevd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
	#endif

	bufsize_z = ldz*n;
	bufsize_a = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL)){
		EXPECT_FALSE( true) << "hpevd_double_parameters object: malloc error.";
		hpevd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, bufsize_a, bufsize_a, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hpevd_dcomplex_parameters :: ~hpevd_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpevd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpevd_free();

}
/*  Test fixture class definition */
class zhpevd_test  : public  ::testing::Test {
public:
   hpevd_dcomplex_parameters  *zhpevd_obj;
   void SetUp();
   void TearDown () { delete zhpevd_obj; }
};

void zhpevd_test::SetUp(){

    /* LAPACKE zhpevd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhpevd) (int matrix_layout, char jobz, char uplo, lapack_int n,lapack_complex_double *A, 
											double* w, lapack_complex_double* z, lapack_int ldz);

    Fptr_NL_LAPACKE_zhpevd zhpevd;

    zhpevd_obj = new hpevd_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhpevd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhpevd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhpevd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhpevd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhpevd = (Fptr_NL_LAPACKE_zhpevd)dlsym(zhpevd_obj->hModule, "LAPACKE_zhpevd");
    ASSERT_TRUE(zhpevd != NULL) << "failed to get the Netlib LAPACKE_zhpevd symbol";
    

    zhpevd_obj->inforef = zhpevd( zhpevd_obj->matrix_layout, zhpevd_obj->jobz, zhpevd_obj->uplo,
								zhpevd_obj->n, zhpevd_obj->Aref, zhpevd_obj->wref, zhpevd_obj->z, zhpevd_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zhpevd_obj->info = LAPACKE_zhpevd( zhpevd_obj->matrix_layout, zhpevd_obj->jobz, zhpevd_obj->uplo,
										zhpevd_obj->n,zhpevd_obj->A, zhpevd_obj->w, zhpevd_obj->z, zhpevd_obj->ldz);

    if( zhpevd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhpevd is wrong\n", zhpevd_obj->info );
    }
    if( zhpevd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhpevd is wrong\n", 
        zhpevd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhpevd_obj->diff_a =  computeDiff_z( zhpevd_obj->bufsize_a, zhpevd_obj->A, zhpevd_obj->Aref );
	zhpevd_obj->diff_w =  computeDiff_d( zhpevd_obj->n, zhpevd_obj->w, zhpevd_obj->wref );
}

TEST_F(zhpevd_test, zhpevd1) {
    EXPECT_NEAR(0.0, zhpevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpevd_test, zhpevd2) {
    EXPECT_NEAR(0.0, zhpevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpevd_test, zhpevd3) {
    EXPECT_NEAR(0.0, zhpevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpevd_test, zhpevd4) {
    EXPECT_NEAR(0.0, zhpevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}
