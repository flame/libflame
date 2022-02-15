#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hpgvd_free() \
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
class hpgvd_scomplex_parameters{

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
      hpgvd_scomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n, lapack_int ldz);
      ~hpgvd_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hpgvd_scomplex_parameters:: hpgvd_scomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ldz_i)
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
	printf(" \n hpgvd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
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
		EXPECT_FALSE( true) << "hpgvd_float_parameters object: malloc error.";
		hpgvd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand(b, bref, bufsize_a);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hpgvd_scomplex_parameters :: ~hpgvd_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpgvd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpgvd_free();
  

}
/*  Test fixture class definition */
class chpgvd_test  : public  ::testing::Test {
public:
   hpgvd_scomplex_parameters  *chpgvd_obj;
   void SetUp();
   void TearDown () { delete chpgvd_obj; }
};

void chpgvd_test::SetUp(){

    /* LAPACKE chpgvd prototype */
    typedef int (*Fptr_NL_LAPACKE_chpgvd) (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,lapack_complex_float *A, 
										lapack_complex_float* b, float* w, lapack_complex_float* z, lapack_int ldz);

    Fptr_NL_LAPACKE_chpgvd chpgvd;

    chpgvd_obj = new hpgvd_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].nrhs,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].n,
							eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chpgvd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chpgvd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chpgvd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chpgvd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chpgvd = (Fptr_NL_LAPACKE_chpgvd)dlsym(chpgvd_obj->hModule, "LAPACKE_chpgvd");
    ASSERT_TRUE(chpgvd != NULL) << "failed to get the Netlib LAPACKE_chpgvd symbol";
    

    chpgvd_obj->inforef = chpgvd( chpgvd_obj->matrix_layout, chpgvd_obj->itype, chpgvd_obj->jobz, chpgvd_obj->uplo, 
								chpgvd_obj->n, chpgvd_obj->Aref, chpgvd_obj->bref, chpgvd_obj->w, chpgvd_obj->zref, chpgvd_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    chpgvd_obj->info = LAPACKE_chpgvd( chpgvd_obj->matrix_layout, chpgvd_obj->itype, chpgvd_obj->jobz, chpgvd_obj->uplo,
										chpgvd_obj->n,chpgvd_obj->A, chpgvd_obj->b, chpgvd_obj->w, chpgvd_obj->z, chpgvd_obj->ldz);

    if( chpgvd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chpgvd is wrong\n", chpgvd_obj->info );
    }
    if( chpgvd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chpgvd is wrong\n", 
        chpgvd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chpgvd_obj->diff_a =  computeDiff_c( chpgvd_obj->bufsize_a, chpgvd_obj->A, chpgvd_obj->Aref );
	chpgvd_obj->diff_b =  computeDiff_c( chpgvd_obj->bufsize_a, chpgvd_obj->b, chpgvd_obj->bref );
	chpgvd_obj->diff_z =  computeDiff_c( chpgvd_obj->bufsize_z, chpgvd_obj->z, chpgvd_obj->zref );
	chpgvd_obj->diff_w =  computeDiff_s( chpgvd_obj->n, chpgvd_obj->w, chpgvd_obj->wref );
}

TEST_F(chpgvd_test, chpgvd1) {
    EXPECT_NEAR(0.0, chpgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgvd_test, chpgvd2) {
    EXPECT_NEAR(0.0, chpgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgvd_test, chpgvd3) {
    EXPECT_NEAR(0.0, chpgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpgvd_test, chpgvd4) {
    EXPECT_NEAR(0.0, chpgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chpgvd_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hpgvd_dcomplex_parameters{

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
      hpgvd_dcomplex_parameters (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n, lapack_int ldz);
      ~hpgvd_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hpgvd_dcomplex_parameters:: hpgvd_dcomplex_parameters (int matrix_layout_i, lapack_int itype_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ldz_i)
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
	printf(" \n hpgvd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
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
		EXPECT_FALSE( true) << "hpgvd_double_parameters object: malloc error.";
		hpgvd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(b, bref, bufsize_a);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hpgvd_dcomplex_parameters :: ~hpgvd_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hpgvd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hpgvd_free();

}
/*  Test fixture class definition */ 
class zhpgvd_test  : public  ::testing::Test {
public:
   hpgvd_dcomplex_parameters  *zhpgvd_obj;
   void SetUp();
   void TearDown () { delete zhpgvd_obj; }
};

void zhpgvd_test::SetUp(){

    /* LAPACKE zhpgvd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhpgvd) (int matrix_layout, lapack_int itype, char jobz, char uplo, lapack_int n,lapack_complex_double *A, 
										lapack_complex_double* b, double* w, lapack_complex_double* z, lapack_int ldz );

    Fptr_NL_LAPACKE_zhpgvd zhpgvd;

    zhpgvd_obj = new hpgvd_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhpgvd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhpgvd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhpgvd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhpgvd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhpgvd = (Fptr_NL_LAPACKE_zhpgvd)dlsym(zhpgvd_obj->hModule, "LAPACKE_zhpgvd");
    ASSERT_TRUE(zhpgvd != NULL) << "failed to get the Netlib LAPACKE_zhpgvd symbol";
    

    zhpgvd_obj->inforef = zhpgvd( zhpgvd_obj->matrix_layout, zhpgvd_obj->itype, zhpgvd_obj->jobz, zhpgvd_obj->uplo,
								zhpgvd_obj->n, zhpgvd_obj->Aref, zhpgvd_obj->bref, zhpgvd_obj->wref, zhpgvd_obj->zref, zhpgvd_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zhpgvd_obj->info = LAPACKE_zhpgvd( zhpgvd_obj->matrix_layout, zhpgvd_obj->itype, zhpgvd_obj->jobz, zhpgvd_obj->uplo,
										zhpgvd_obj->n,zhpgvd_obj->A, zhpgvd_obj->b, zhpgvd_obj->w, zhpgvd_obj->z, zhpgvd_obj->ldz);

    if( zhpgvd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhpgvd is wrong\n", zhpgvd_obj->info );
    }
    if( zhpgvd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhpgvd is wrong\n", 
        zhpgvd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhpgvd_obj->diff_a =  computeDiff_z( zhpgvd_obj->bufsize_a, zhpgvd_obj->A, zhpgvd_obj->Aref );
	zhpgvd_obj->diff_b =  computeDiff_z( zhpgvd_obj->bufsize_a, zhpgvd_obj->b, zhpgvd_obj->bref );
	zhpgvd_obj->diff_z =  computeDiff_z( zhpgvd_obj->bufsize_z, zhpgvd_obj->z, zhpgvd_obj->zref );
	zhpgvd_obj->diff_w =  computeDiff_d( zhpgvd_obj->n, zhpgvd_obj->w, zhpgvd_obj->wref );
}

TEST_F(zhpgvd_test, zhpgvd1) {
    EXPECT_NEAR(0.0, zhpgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgvd_test, zhpgvd2) {
    EXPECT_NEAR(0.0, zhpgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgvd_test, zhpgvd3) {
    EXPECT_NEAR(0.0, zhpgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpgvd_test, zhpgvd4) {
    EXPECT_NEAR(0.0, zhpgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhpgvd_obj->diff_z, LAPACKE_GTEST_THRESHOLD);
}
