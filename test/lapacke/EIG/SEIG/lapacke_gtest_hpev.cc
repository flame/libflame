#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define hpev_free() \
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
class hpev_scomplex_parameters{

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
      hpev_scomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ldz);
      ~hpev_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hpev_scomplex_parameters:: hpev_scomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	
	if (jobz == 'U')
		jobz = 'N';
	
	//#if LAPACKE_TEST_VERBOSE
	printf(" \n hpev scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
//	#endif

	bufsize_z = ldz*n;
	bufsize_a = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL)){
		EXPECT_FALSE( true) << "hpev_float_parameters object: malloc error.";
		hpev_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);
	
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hpev_scomplex_parameters :: ~hpev_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpev_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpev_free();
  

}
/*  Test fixture class definition */
class chpev_test  : public  ::testing::Test {
public:
   hpev_scomplex_parameters  *chpev_obj;
   void SetUp();
   void TearDown () { delete chpev_obj; }
};

void chpev_test::SetUp(){

    /* LAPACKE chpev prototype */
    typedef int (*Fptr_NL_LAPACKE_chpev) (int matrix_layout, char jobz, char uplo, lapack_int n,lapack_complex_float *A, 
										float* w, lapack_complex_float* z, lapack_int ldz);

    Fptr_NL_LAPACKE_chpev chpev;

    chpev_obj = new hpev_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chpev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chpev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chpev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chpev_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chpev = (Fptr_NL_LAPACKE_chpev)dlsym(chpev_obj->hModule, "LAPACKE_chpev");
    ASSERT_TRUE(chpev != NULL) << "failed to get the Netlib LAPACKE_chpev symbol";
    

    chpev_obj->inforef = chpev( chpev_obj->matrix_layout, chpev_obj->jobz, chpev_obj->uplo,
								chpev_obj->n, chpev_obj->Aref, chpev_obj->wref, chpev_obj->zref, chpev_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    chpev_obj->info = LAPACKE_chpev( chpev_obj->matrix_layout, chpev_obj->jobz, chpev_obj->uplo,
										chpev_obj->n,chpev_obj->A, chpev_obj->w, chpev_obj->z, chpev_obj->ldz);

    if( chpev_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chpev is wrong\n", chpev_obj->info );
    }
    if( chpev_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chpev is wrong\n", 
        chpev_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chpev_obj->diff_a =  computeDiff_c( chpev_obj->bufsize_a, chpev_obj->A, chpev_obj->Aref );
	chpev_obj->diff_w =  computeDiff_s( chpev_obj->n, chpev_obj->w, chpev_obj->wref );
}

TEST_F(chpev_test, chpev1) {
    EXPECT_NEAR(0.0, chpev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpev_test, chpev2) {
    EXPECT_NEAR(0.0, chpev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpev_test, chpev3) {
    EXPECT_NEAR(0.0, chpev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpev_test, chpev4) {
    EXPECT_NEAR(0.0, chpev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hpev_dcomplex_parameters{

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
      hpev_dcomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ldz);
      ~hpev_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hpev_dcomplex_parameters:: hpev_dcomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hpev scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
	#endif

	bufsize_z = ldz*n;
	bufsize_a = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(w == NULL) || (wref == NULL)){
		EXPECT_FALSE( true) << "hpev_double_parameters object: malloc error.";
		hpev_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hpev_dcomplex_parameters :: ~hpev_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpev_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpev_free();

}
/*  Test fixture class definition */
class zhpev_test  : public  ::testing::Test {
public:
   hpev_dcomplex_parameters  *zhpev_obj;
   void SetUp();
   void TearDown () { delete zhpev_obj; }
};

void zhpev_test::SetUp(){

    /* LAPACKE zhpev prototype */
    typedef int (*Fptr_NL_LAPACKE_zhpev) (int matrix_layout, char jobz, char uplo, lapack_int n,lapack_complex_double *A, 
											double* w, lapack_complex_double* z, lapack_int ldz);

    Fptr_NL_LAPACKE_zhpev zhpev;

    zhpev_obj = new hpev_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhpev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhpev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhpev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhpev_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhpev = (Fptr_NL_LAPACKE_zhpev)dlsym(zhpev_obj->hModule, "LAPACKE_zhpev");
    ASSERT_TRUE(zhpev != NULL) << "failed to get the Netlib LAPACKE_zhpev symbol";
    

    zhpev_obj->inforef = zhpev( zhpev_obj->matrix_layout, zhpev_obj->jobz, zhpev_obj->uplo,
								zhpev_obj->n, zhpev_obj->Aref, zhpev_obj->wref, zhpev_obj->z, zhpev_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zhpev_obj->info = LAPACKE_zhpev( zhpev_obj->matrix_layout, zhpev_obj->jobz, zhpev_obj->uplo,
										zhpev_obj->n,zhpev_obj->A, zhpev_obj->w, zhpev_obj->z, zhpev_obj->ldz);

    if( zhpev_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhpev is wrong\n", zhpev_obj->info );
    }
    if( zhpev_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhpev is wrong\n", 
        zhpev_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhpev_obj->diff_a =  computeDiff_z( zhpev_obj->bufsize_a, zhpev_obj->A, zhpev_obj->Aref );
	zhpev_obj->diff_w =  computeDiff_d( zhpev_obj->n, zhpev_obj->w, zhpev_obj->wref );
}

TEST_F(zhpev_test, zhpev1) {
    EXPECT_NEAR(0.0, zhpev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpev_test, zhpev2) {
    EXPECT_NEAR(0.0, zhpev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpev_test, zhpev3) {
    EXPECT_NEAR(0.0, zhpev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpev_test, zhpev4) {
    EXPECT_NEAR(0.0, zhpev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpev_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}
