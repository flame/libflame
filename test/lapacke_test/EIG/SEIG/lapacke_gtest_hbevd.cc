#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hbevd_free() \
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
class hbevd_scomplex_parameters{

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
      hbevd_scomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int kd, lapack_int ldz);
      ~hbevd_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hbevd_scomplex_parameters:: hbevd_scomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int kd_i, lapack_int ldz_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	kd =kd_i;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hbevd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d ,  kd:%d\n", matrix_layout, jobz,uplo, n, ldz, kd);
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
		EXPECT_FALSE( true) << "hbevd_float_parameters object: malloc error.";
		hbevd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n,n , uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);
	
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hbevd_scomplex_parameters :: ~hbevd_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbevd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbevd_free();
  

}
/*  Test fixture class definition */
class chbevd_test  : public  ::testing::Test {
public:
   hbevd_scomplex_parameters  *chbevd_obj;
   void SetUp();
   void TearDown () { delete chbevd_obj; }
};

void chbevd_test::SetUp(){

    /* LAPACKE chbevd prototype */
    typedef int (*Fptr_NL_LAPACKE_chbevd) (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int kd, lapack_complex_float *A, 
										lapack_int ldab, float* w, lapack_complex_float* z, lapack_int ldz);

    Fptr_NL_LAPACKE_chbevd chbevd;

    chbevd_obj = new hbevd_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chbevd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chbevd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chbevd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chbevd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chbevd = (Fptr_NL_LAPACKE_chbevd)dlsym(chbevd_obj->hModule, "LAPACKE_chbevd");
    ASSERT_TRUE(chbevd != NULL) << "failed to get the Netlib LAPACKE_chbevd symbol";
    

    chbevd_obj->inforef = chbevd( chbevd_obj->matrix_layout, chbevd_obj->jobz, chbevd_obj->uplo, chbevd_obj->n,
								chbevd_obj->kd, chbevd_obj->Aref, chbevd_obj->ldab, chbevd_obj->wref, chbevd_obj->zref, chbevd_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    chbevd_obj->info = LAPACKE_chbevd( chbevd_obj->matrix_layout, chbevd_obj->jobz, chbevd_obj->uplo, chbevd_obj->n,
										chbevd_obj->kd,chbevd_obj->A, chbevd_obj->ldab, chbevd_obj->w, chbevd_obj->z, chbevd_obj->ldz);

    if( chbevd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chbevd is wrong\n", chbevd_obj->info );
    }
    if( chbevd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chbevd is wrong\n", 
        chbevd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chbevd_obj->diff_a =  computeDiff_c( chbevd_obj->bufsize_a, chbevd_obj->A, chbevd_obj->Aref );
	
}

TEST_F(chbevd_test, chbevd1) {
    EXPECT_NEAR(0.0, chbevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);	
}

TEST_F(chbevd_test, chbevd2) {
    EXPECT_NEAR(0.0, chbevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbevd_test, chbevd3) {
    EXPECT_NEAR(0.0, chbevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbevd_test, chbevd4) {
    EXPECT_NEAR(0.0, chbevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hbevd_dcomplex_parameters{

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
      hbevd_dcomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int kd, lapack_int ldz);
      ~hbevd_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hbevd_dcomplex_parameters:: hbevd_dcomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int kd_i, lapack_int ldz_i)
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
	printf(" \n hbevd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
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
		EXPECT_FALSE( true) << "hbevd_double_parameters object: malloc error.";
		hbevd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hbevd_dcomplex_parameters :: ~hbevd_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbevd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbevd_free();

}
/*  Test fixture class definition */
class zhbevd_test  : public  ::testing::Test {
public:
   hbevd_dcomplex_parameters  *zhbevd_obj;
   void SetUp();
   void TearDown () { delete zhbevd_obj; }
};

void zhbevd_test::SetUp(){

    /* LAPACKE zhbevd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhbevd) (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int kd, lapack_complex_double *A, 
											lapack_int ldab, double* w, lapack_complex_double* z, lapack_int ldz);

    Fptr_NL_LAPACKE_zhbevd zhbevd;

    zhbevd_obj = new hbevd_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhbevd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhbevd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhbevd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhbevd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhbevd = (Fptr_NL_LAPACKE_zhbevd)dlsym(zhbevd_obj->hModule, "LAPACKE_zhbevd");
    ASSERT_TRUE(zhbevd != NULL) << "failed to get the Netlib LAPACKE_zhbevd symbol";
    

    zhbevd_obj->inforef = zhbevd( zhbevd_obj->matrix_layout, zhbevd_obj->jobz, zhbevd_obj->uplo, zhbevd_obj->n, 
								zhbevd_obj->kd, zhbevd_obj->Aref, zhbevd_obj->ldab, zhbevd_obj->wref, zhbevd_obj->z, zhbevd_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zhbevd_obj->info = LAPACKE_zhbevd( zhbevd_obj->matrix_layout, zhbevd_obj->jobz, zhbevd_obj->uplo, zhbevd_obj->n,
									zhbevd_obj->kd, zhbevd_obj->A, zhbevd_obj->ldab, zhbevd_obj->w, zhbevd_obj->z, zhbevd_obj->ldz);

    if( zhbevd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhbevd is wrong\n", zhbevd_obj->info );
    }
    if( zhbevd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhbevd is wrong\n", 
        zhbevd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhbevd_obj->diff_a =  computeDiff_z( zhbevd_obj->bufsize_a, zhbevd_obj->A, zhbevd_obj->Aref );
}

TEST_F(zhbevd_test, zhbevd1) {
    EXPECT_NEAR(0.0, zhbevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbevd_test, zhbevd2) {
    EXPECT_NEAR(0.0, zhbevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbevd_test, zhbevd3) {
    EXPECT_NEAR(0.0, zhbevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbevd_test, zhbevd4) {
    EXPECT_NEAR(0.0, zhbevd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
}
