#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define heevd_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class heevd_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff_a;
	float diff_w;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* A;	
	char uplo;
	char jobz;
	lapack_int lda;
	/*Output Parameter*/
	float* w;	
	lapack_complex_float *Aref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      heevd_scomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int lda);
      ~heevd_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
heevd_scomplex_parameters:: heevd_scomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;	
	lda = lda_i;
	uplo = uplo_i;
	jobz = jobz_i;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n heevd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d lda: %d \n", matrix_layout, jobz,uplo, n, lda);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL)){
		EXPECT_FALSE( true) << "heevd_float_parameters object: malloc error.";
		heevd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	/*initialize output matrix by 0 */
	for(i=0;i<n;i++) {
		w[i] = 0;
		wref[i] = w[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
heevd_scomplex_parameters :: ~heevd_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heevd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heevd_free();

}
/*  Test fixture class definition */
class cheevd_test  : public  ::testing::Test {
public:
   heevd_scomplex_parameters  *cheevd_obj;
   void SetUp();
   void TearDown () { delete cheevd_obj; }
};

void cheevd_test::SetUp(){

    /* LAPACKE cheevd prototype */
    typedef int (*Fptr_NL_LAPACKE_cheevd) (int matrix_layout, char jobz, char uplo, lapack_int n,lapack_complex_float *A, 
											lapack_int lda, float* w);

    Fptr_NL_LAPACKE_cheevd cheevd;

    cheevd_obj = new heevd_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    cheevd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cheevd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cheevd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cheevd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cheevd = (Fptr_NL_LAPACKE_cheevd)dlsym(cheevd_obj->hModule, "LAPACKE_cheevd");
    ASSERT_TRUE(cheevd != NULL) << "failed to get the Netlib LAPACKE_cheevd symbol";
    

    cheevd_obj->inforef = cheevd( cheevd_obj->matrix_layout, cheevd_obj->jobz, cheevd_obj->uplo,
								cheevd_obj->n, cheevd_obj->Aref, cheevd_obj->lda, cheevd_obj->wref);

    /* Compute libflame's Lapacke o/p  */
    cheevd_obj->info = LAPACKE_cheevd( cheevd_obj->matrix_layout, cheevd_obj->jobz, cheevd_obj->uplo,
										cheevd_obj->n,cheevd_obj->A, cheevd_obj->lda, cheevd_obj->w);

    if( cheevd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cheevd is wrong\n", cheevd_obj->info );
    }
    if( cheevd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cheevd is wrong\n", 
        cheevd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cheevd_obj->diff_a =  computeDiff_c( cheevd_obj->bufsize, cheevd_obj->A, cheevd_obj->Aref );
	cheevd_obj->diff_w =  computeDiff_s( cheevd_obj->n, cheevd_obj->w, cheevd_obj->wref );
}

TEST_F(cheevd_test, cheevd1) {
    EXPECT_NEAR(0.0, cheevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevd_test, cheevd2) {
    EXPECT_NEAR(0.0, cheevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevd_test, cheevd3) {
    EXPECT_NEAR(0.0, cheevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheevd_test, cheevd4) {
    EXPECT_NEAR(0.0, cheevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cheevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class heevd_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	 double diff_a;
	double diff_w;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* A;	
	char uplo;
	char jobz;
	lapack_int lda;
	/*Output Parameter*/
	double* w;	
	lapack_complex_double *Aref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      heevd_dcomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int lda);
      ~heevd_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
heevd_dcomplex_parameters:: heevd_dcomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;	
	lda = lda_i;
	uplo = uplo_i;
	jobz = jobz_i;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n heevd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d lda: %d \n", matrix_layout, jobz,uplo, n, lda);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL)){
		EXPECT_FALSE( true) << "heevd_double_parameters object: malloc error.";
		heevd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	/*initialize output matrix by 0 */
	for(i=0;i<n;i++) {
		w[i] = 0;
		wref[i] = w[i];
	}

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
heevd_dcomplex_parameters :: ~heevd_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heevd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heevd_free();

}
/*  Test fixture class definition */
class zheevd_test  : public  ::testing::Test {
public:
   heevd_dcomplex_parameters  *zheevd_obj;
   void SetUp();
   void TearDown () { delete zheevd_obj; }
};

void zheevd_test::SetUp(){

    /* LAPACKE zheevd prototype */
    typedef int (*Fptr_NL_LAPACKE_zheevd) (int matrix_layout, char jobz, char uplo, lapack_int n,lapack_complex_double *A, 
											lapack_int lda, double* w);

    Fptr_NL_LAPACKE_zheevd zheevd;

    zheevd_obj = new heevd_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zheevd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zheevd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zheevd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zheevd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zheevd = (Fptr_NL_LAPACKE_zheevd)dlsym(zheevd_obj->hModule, "LAPACKE_zheevd");
    ASSERT_TRUE(zheevd != NULL) << "failed to get the Netlib LAPACKE_zheevd symbol";
    

    zheevd_obj->inforef = zheevd( zheevd_obj->matrix_layout, zheevd_obj->jobz, zheevd_obj->uplo,
								zheevd_obj->n, zheevd_obj->Aref, zheevd_obj->lda, zheevd_obj->wref);

    /* Compute libflame's Lapacke o/p  */
    zheevd_obj->info = LAPACKE_zheevd( zheevd_obj->matrix_layout, zheevd_obj->jobz, zheevd_obj->uplo,
										zheevd_obj->n,zheevd_obj->A, zheevd_obj->lda, zheevd_obj->w);

    if( zheevd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zheevd is wrong\n", zheevd_obj->info );
    }
    if( zheevd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zheevd is wrong\n", 
        zheevd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zheevd_obj->diff_a =  computeDiff_z( zheevd_obj->bufsize, zheevd_obj->A, zheevd_obj->Aref );
	zheevd_obj->diff_w =  computeDiff_d( zheevd_obj->n, zheevd_obj->w, zheevd_obj->wref );
}

TEST_F(zheevd_test, zheevd1) {
    EXPECT_NEAR(0.0, zheevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheevd_test, zheevd2) {
    EXPECT_NEAR(0.0, zheevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheevd_test, zheevd3) {
    EXPECT_NEAR(0.0, zheevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheevd_test, zheevd4) {
    EXPECT_NEAR(0.0, zheevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zheevd_obj->diff_w, LAPACKE_GTEST_THRESHOLD);
}
