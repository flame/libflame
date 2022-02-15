#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"

#define hbgvd_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)    free(b); \
if (bref!=NULL) free(bref);\
if (w!=NULL)  free(w);\
if (wref!=NULL) free(wref); \
if (z!=NULL)  free(z);\
if (zref!=NULL)  free(zref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hbgvd_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_z;
	void *hModule, *dModule;
	float diff_a;
	float diff_w;
	float diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ka;
	lapack_int kb;
	lapack_int ldab;
	lapack_int ldbb;
	lapack_complex_float* A, *z, *b;
	char uplo;
	char jobz;
	lapack_int ldz;
	/*Output Parameter*/
	float* w;	
	lapack_complex_float *Aref, *zref, *bref;
	float* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hbgvd_scomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ka, lapack_int ldz);
      ~hbgvd_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hbgvd_scomplex_parameters:: hbgvd_scomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ka_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	ka =ka_i;
	kb = ka;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hbgvd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d ,  kd:%d\n", matrix_layout, jobz,uplo, n, ldz, kd);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	
		ldab = ka +1;
		ldbb = ldab;
		bufsize_a  = (ldab * n);
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		ldab = n;
		ldbb = ldab;
		bufsize_a = ldab*(ka + 1);
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	bufsize_z = ldz*n;

	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_float_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL) ||
		(b == NULL) || (bref == NULL)){
		EXPECT_FALSE( true) << "hbgvd_float_parameters object: malloc error.";
		hbgvd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n,n , uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand(b, bref, bufsize_a);
	
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hbgvd_scomplex_parameters :: ~hbgvd_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbgvd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbgvd_free();
  

}
/*  Test fixture class definition */
class chbgvd_test  : public  ::testing::Test {
public:
   hbgvd_scomplex_parameters  *chbgvd_obj;
   void SetUp();
   void TearDown () { delete chbgvd_obj; }
};

void chbgvd_test::SetUp(){

    /* LAPACKE chbgvd prototype */
    typedef int (*Fptr_NL_LAPACKE_chbgvd) (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ka, lapack_int kb, 
										lapack_complex_float *A, lapack_int ldab, lapack_complex_float* b, lapack_int ldbb, 
										float* w, lapack_complex_float* z, lapack_int ldz);

    Fptr_NL_LAPACKE_chbgvd chbgvd;

    chbgvd_obj = new hbgvd_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chbgvd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chbgvd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chbgvd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chbgvd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chbgvd = (Fptr_NL_LAPACKE_chbgvd)dlsym(chbgvd_obj->hModule, "LAPACKE_chbgvd");
    ASSERT_TRUE(chbgvd != NULL) << "failed to get the Netlib LAPACKE_chbgvd symbol";
    

    chbgvd_obj->inforef = chbgvd( chbgvd_obj->matrix_layout, chbgvd_obj->jobz, chbgvd_obj->uplo, chbgvd_obj->n, chbgvd_obj->ka, 
								chbgvd_obj->kb, chbgvd_obj->Aref, chbgvd_obj->ldab, chbgvd_obj->bref, chbgvd_obj->ldbb, 
								chbgvd_obj->wref, chbgvd_obj->zref, chbgvd_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    chbgvd_obj->info = LAPACKE_chbgvd( chbgvd_obj->matrix_layout, chbgvd_obj->jobz, chbgvd_obj->uplo, chbgvd_obj->n, chbgvd_obj->ka,
										chbgvd_obj->kb,chbgvd_obj->A, chbgvd_obj->ldab, chbgvd_obj->b, chbgvd_obj->ldbb,
										chbgvd_obj->w, chbgvd_obj->z, chbgvd_obj->ldz);

    if( chbgvd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chbgvd is wrong\n", chbgvd_obj->info );
    }
    if( chbgvd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chbgvd is wrong\n", 
        chbgvd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chbgvd_obj->diff_a =  computeDiff_c( chbgvd_obj->bufsize_a, chbgvd_obj->A, chbgvd_obj->Aref );
	chbgvd_obj->diff_b =  computeDiff_c( chbgvd_obj->bufsize_a, chbgvd_obj->b, chbgvd_obj->bref );
	
}

TEST_F(chbgvd_test, chbgvd1) {
    EXPECT_NEAR(0.0, chbgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbgvd_test, chbgvd2) {
    EXPECT_NEAR(0.0, chbgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbgvd_test, chbgvd3) {
    EXPECT_NEAR(0.0, chbgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbgvd_test, chbgvd4) {
    EXPECT_NEAR(0.0, chbgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hbgvd_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_z;
	void *hModule, *dModule;
	 double diff_a;
	double diff_w;
	double diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int ldab;
	lapack_int ldbb;
	lapack_int kb;
	lapack_int ka;
	lapack_complex_double* A, *z, *b;
	char uplo;
	char jobz;
	lapack_int ldz;
	/*Output Parameter*/
	double* w;	
	lapack_complex_double *Aref, *zref, *bref;
	double* wref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hbgvd_dcomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ka, lapack_int ldz);
      ~hbgvd_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hbgvd_dcomplex_parameters:: hbgvd_dcomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ka_i, lapack_int ldz_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;	
	ldz = ldz_i;
	uplo = uplo_i;
	jobz = jobz_i;
	ka = ka_i;
	kb = ka;
	
	if (jobz == 'U')
		jobz = 'N';
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hbgvd scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	
		ldab = ka +1;
		ldbb = ldab;
		bufsize_a  = (ldab * n);
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		ldab = n;
		ldbb = ldab;
		bufsize_a = ldab*(ka + 1);
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";

	bufsize_z = ldz*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&z, &zref, bufsize_z);
	lapacke_gtest_alloc_double_buffer_pair(&w, &wref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(w == NULL) || (wref == NULL) ||
		(z == NULL) || (zref == NULL)||
		(b == NULL) || (bref == NULL)){
		EXPECT_FALSE( true) << "hbgvd_double_parameters object: malloc error.";
		hbgvd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(b, bref, bufsize_a);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hbgvd_dcomplex_parameters :: ~hbgvd_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbgvd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbgvd_free();

}
/*  Test fixture class definition */
class zhbgvd_test  : public  ::testing::Test {
public:
   hbgvd_dcomplex_parameters  *zhbgvd_obj;
   void SetUp();
   void TearDown () { delete zhbgvd_obj; }
};

void zhbgvd_test::SetUp(){

    /* LAPACKE zhbgvd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhbgvd) (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ka, lapack_int kb, 
											lapack_complex_double *A, lapack_int ldab, lapack_complex_double* b, lapack_int ldbb, 
											double* w, lapack_complex_double* z, lapack_int ldz);

    Fptr_NL_LAPACKE_zhbgvd zhbgvd;

    zhbgvd_obj = new hbgvd_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhbgvd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhbgvd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhbgvd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhbgvd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhbgvd = (Fptr_NL_LAPACKE_zhbgvd)dlsym(zhbgvd_obj->hModule, "LAPACKE_zhbgvd");
    ASSERT_TRUE(zhbgvd != NULL) << "failed to get the Netlib LAPACKE_zhbgvd symbol";
    

    zhbgvd_obj->inforef = zhbgvd( zhbgvd_obj->matrix_layout, zhbgvd_obj->jobz, zhbgvd_obj->uplo, zhbgvd_obj->n, 
								zhbgvd_obj->ka, zhbgvd_obj->kb, zhbgvd_obj->Aref, zhbgvd_obj->ldab, zhbgvd_obj->bref, zhbgvd_obj->ldbb,
								zhbgvd_obj->wref, zhbgvd_obj->z, zhbgvd_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zhbgvd_obj->info = LAPACKE_zhbgvd( zhbgvd_obj->matrix_layout, zhbgvd_obj->jobz, zhbgvd_obj->uplo, zhbgvd_obj->n, zhbgvd_obj->ka,
									zhbgvd_obj->kb, zhbgvd_obj->A, zhbgvd_obj->ldab, zhbgvd_obj->b, zhbgvd_obj->ldbb,
									zhbgvd_obj->w, zhbgvd_obj->z, zhbgvd_obj->ldz);

    if( zhbgvd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhbgvd is wrong\n", zhbgvd_obj->info );
    }
    if( zhbgvd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhbgvd is wrong\n", 
        zhbgvd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhbgvd_obj->diff_a =  computeDiff_z( zhbgvd_obj->bufsize_a, zhbgvd_obj->A, zhbgvd_obj->Aref );
	zhbgvd_obj->diff_b =  computeDiff_z( zhbgvd_obj->bufsize_a, zhbgvd_obj->b, zhbgvd_obj->bref );
}

TEST_F(zhbgvd_test, zhbgvd1) {
    EXPECT_NEAR(0.0, zhbgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbgvd_test, zhbgvd2) {
    EXPECT_NEAR(0.0, zhbgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbgvd_test, zhbgvd3) {
    EXPECT_NEAR(0.0, zhbgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbgvd_test, zhbgvd4) {
    EXPECT_NEAR(0.0, zhbgvd_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgvd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}
