#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hbev_free() \
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
class hbev_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_z;
	void *hModule, *dModule;
	float diff_a;
	float diff_w;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int kd;
	lapack_int ldab;
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
      hbev_scomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int kd, lapack_int ldz);
      ~hbev_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hbev_scomplex_parameters:: hbev_scomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int kd_i, lapack_int ldz_i)
{
	;
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	kd =kd_i;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hbev scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d ,  kd:%d\n", matrix_layout, jobz,uplo, n, ldz, kd);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	
		ldab = kd +1;
		bufsize_a  = (ldab * n);
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		ldab = n;
		bufsize_a = ldab*(kd + 1);
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	bufsize_z = ldz*n;

	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL)){
		EXPECT_FALSE( true) << "hbev_float_parameters object: malloc error.";
		hbev_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n,n , uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);
	
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hbev_scomplex_parameters :: ~hbev_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbev_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbev_free();
  

}
/*  Test fixture class definition */
class chbev_test  : public  ::testing::Test {
public:
   hbev_scomplex_parameters  *chbev_obj;
   void SetUp();
   void TearDown () { delete chbev_obj; }
};

void chbev_test::SetUp(){

    /* LAPACKE chbev prototype */
    typedef int (*Fptr_NL_LAPACKE_chbev) (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int kd, lapack_complex_float *A, 
										lapack_int ldab, float* w, lapack_complex_float* z, lapack_int ldz);

    Fptr_NL_LAPACKE_chbev chbev;

    chbev_obj = new hbev_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chbev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chbev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chbev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chbev_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chbev = (Fptr_NL_LAPACKE_chbev)dlsym(chbev_obj->hModule, "LAPACKE_chbev");
    ASSERT_TRUE(chbev != NULL) << "failed to get the Netlib LAPACKE_chbev symbol";
    

    chbev_obj->inforef = chbev( chbev_obj->matrix_layout, chbev_obj->jobz, chbev_obj->uplo, chbev_obj->n,
								chbev_obj->kd, chbev_obj->Aref, chbev_obj->ldab, chbev_obj->wref, chbev_obj->zref, chbev_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    chbev_obj->info = LAPACKE_chbev( chbev_obj->matrix_layout, chbev_obj->jobz, chbev_obj->uplo, chbev_obj->n,
										chbev_obj->kd,chbev_obj->A, chbev_obj->ldab, chbev_obj->w, chbev_obj->z, chbev_obj->ldz);

    if( chbev_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chbev is wrong\n", chbev_obj->info );
    }
    if( chbev_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chbev is wrong\n", 
        chbev_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chbev_obj->diff_a =  computeDiff_c( chbev_obj->bufsize_a, chbev_obj->A, chbev_obj->Aref );
	
}

TEST_F(chbev_test, chbev1) {
    EXPECT_NEAR(0.0, chbev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);	
}

TEST_F(chbev_test, chbev2) {
    EXPECT_NEAR(0.0, chbev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbev_test, chbev3) {
    EXPECT_NEAR(0.0, chbev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbev_test, chbev4) {
    EXPECT_NEAR(0.0, chbev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hbev_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_z;
	void *hModule, *dModule;
	 double diff_a;
	double diff_w;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldab;
	lapack_int kd;
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
      hbev_dcomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int kd, lapack_int ldz);
      ~hbev_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hbev_dcomplex_parameters:: hbev_dcomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int kd_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	kd = kd_i;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hbev scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	
		ldab = kd +1;
		bufsize_a  = (ldab * n);
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		ldab = n;
		bufsize_a = ldab*(kd + 1);
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	bufsize_z = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL)){
		EXPECT_FALSE( true) << "hbev_double_parameters object: malloc error.";
		hbev_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hbev_dcomplex_parameters :: ~hbev_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbev_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbev_free();

}
/*  Test fixture class definition */
class zhbev_test  : public  ::testing::Test {
public:
   hbev_dcomplex_parameters  *zhbev_obj;
   void SetUp();
   void TearDown () { delete zhbev_obj; }
};

void zhbev_test::SetUp(){

    /* LAPACKE zhbev prototype */
    typedef int (*Fptr_NL_LAPACKE_zhbev) (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int kd, lapack_complex_double *A, 
											lapack_int ldab, double* w, lapack_complex_double* z, lapack_int ldz);

    Fptr_NL_LAPACKE_zhbev zhbev;

    zhbev_obj = new hbev_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhbev_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhbev_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhbev_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhbev_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhbev = (Fptr_NL_LAPACKE_zhbev)dlsym(zhbev_obj->hModule, "LAPACKE_zhbev");
    ASSERT_TRUE(zhbev != NULL) << "failed to get the Netlib LAPACKE_zhbev symbol";
    

    zhbev_obj->inforef = zhbev( zhbev_obj->matrix_layout, zhbev_obj->jobz, zhbev_obj->uplo, zhbev_obj->n, 
								zhbev_obj->kd, zhbev_obj->Aref, zhbev_obj->ldab, zhbev_obj->wref, zhbev_obj->z, zhbev_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zhbev_obj->info = LAPACKE_zhbev( zhbev_obj->matrix_layout, zhbev_obj->jobz, zhbev_obj->uplo, zhbev_obj->n,
									zhbev_obj->kd, zhbev_obj->A, zhbev_obj->ldab, zhbev_obj->w, zhbev_obj->z, zhbev_obj->ldz);

    if( zhbev_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhbev is wrong\n", zhbev_obj->info );
    }
    if( zhbev_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhbev is wrong\n", 
        zhbev_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhbev_obj->diff_a =  computeDiff_z( zhbev_obj->bufsize_a, zhbev_obj->A, zhbev_obj->Aref );
}

TEST_F(zhbev_test, zhbev1) {
    EXPECT_NEAR(0.0, zhbev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbev_test, zhbev2) {
    EXPECT_NEAR(0.0, zhbev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbev_test, zhbev3) {
    EXPECT_NEAR(0.0, zhbev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbev_test, zhbev4) {
    EXPECT_NEAR(0.0, zhbev_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}
