#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define hbgv_free() \
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
class hbgv_scomplex_parameters{

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
      hbgv_scomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ka, lapack_int ldz);
      ~hbgv_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hbgv_scomplex_parameters:: hbgv_scomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ka_i, lapack_int ldz_i)
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
	printf(" \n hbgv scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d ,  kd:%d\n", matrix_layout, jobz,uplo, n, ldz, kd);
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
		EXPECT_FALSE( true) << "hbgv_float_parameters object: malloc error.";
		hbgv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(A, Aref, n,n , uplo);
	lapacke_gtest_init_scomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_scomplex_buffer_pair_rand(b, bref, bufsize_a);
	
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hbgv_scomplex_parameters :: ~hbgv_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbgv_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbgv_free();
  

}
/*  Test fixture class definition */
class chbgv_test  : public  ::testing::Test {
public:
   hbgv_scomplex_parameters  *chbgv_obj;
   void SetUp();
   void TearDown () { delete chbgv_obj; }
};

void chbgv_test::SetUp(){

    /* LAPACKE chbgv prototype */
    typedef int (*Fptr_NL_LAPACKE_chbgv) (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ka, lapack_int kb, 
										lapack_complex_float *A, lapack_int ldab, lapack_complex_float* b, lapack_int ldbb, 
										float* w, lapack_complex_float* z, lapack_int ldz);

    Fptr_NL_LAPACKE_chbgv chbgv;

    chbgv_obj = new hbgv_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chbgv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chbgv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chbgv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chbgv_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chbgv = (Fptr_NL_LAPACKE_chbgv)dlsym(chbgv_obj->hModule, "LAPACKE_chbgv");
    ASSERT_TRUE(chbgv != NULL) << "failed to get the Netlib LAPACKE_chbgv symbol";
    

    chbgv_obj->inforef = chbgv( chbgv_obj->matrix_layout, chbgv_obj->jobz, chbgv_obj->uplo, chbgv_obj->n, chbgv_obj->ka, 
								chbgv_obj->kb, chbgv_obj->Aref, chbgv_obj->ldab, chbgv_obj->bref, chbgv_obj->ldbb, 
								chbgv_obj->wref, chbgv_obj->zref, chbgv_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    chbgv_obj->info = LAPACKE_chbgv( chbgv_obj->matrix_layout, chbgv_obj->jobz, chbgv_obj->uplo, chbgv_obj->n, chbgv_obj->ka,
										chbgv_obj->kb,chbgv_obj->A, chbgv_obj->ldab, chbgv_obj->b, chbgv_obj->ldbb,
										chbgv_obj->w, chbgv_obj->z, chbgv_obj->ldz);

    if( chbgv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chbgv is wrong\n", chbgv_obj->info );
    }
    if( chbgv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chbgv is wrong\n", 
        chbgv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chbgv_obj->diff_a =  computeDiff_c( chbgv_obj->bufsize_a, chbgv_obj->A, chbgv_obj->Aref );
	chbgv_obj->diff_b =  computeDiff_c( chbgv_obj->bufsize_a, chbgv_obj->b, chbgv_obj->bref );
	
}

TEST_F(chbgv_test, chbgv1) {
    EXPECT_NEAR(0.0, chbgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbgv_test, chbgv2) {
    EXPECT_NEAR(0.0, chbgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbgv_test, chbgv3) {
    EXPECT_NEAR(0.0, chbgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbgv_test, chbgv4) {
    EXPECT_NEAR(0.0, chbgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class hbgv_dcomplex_parameters{

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
      hbgv_dcomplex_parameters (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ka, lapack_int ldz);
      ~hbgv_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hbgv_dcomplex_parameters:: hbgv_dcomplex_parameters (int matrix_layout_i, char jobz_i, char uplo_i, lapack_int n_i, lapack_int ka_i, lapack_int ldz_i)
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
	printf(" \n hbgv scomplex: matrix_layout = %d, jobz:%c, uplo:%c, n: %d ldz: %d \n", matrix_layout, jobz,uplo, n, ldz);
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
		EXPECT_FALSE( true) << "hbgv_double_parameters object: malloc error.";
		hbgv_free();
		exit(0);
	}
	/* Initialization of input matrices */
	
	//lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(A, Aref, n, n, uplo);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(A, Aref, bufsize_a);
	lapacke_gtest_init_dcomplex_buffer_pair_rand(b, bref, bufsize_a);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hbgv_dcomplex_parameters :: ~hbgv_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbgv_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbgv_free();

}
/*  Test fixture class definition */
class zhbgv_test  : public  ::testing::Test {
public:
   hbgv_dcomplex_parameters  *zhbgv_obj;
   void SetUp();
   void TearDown () { delete zhbgv_obj; }
};

void zhbgv_test::SetUp(){

    /* LAPACKE zhbgv prototype */
    typedef int (*Fptr_NL_LAPACKE_zhbgv) (int matrix_layout, char jobz, char uplo, lapack_int n, lapack_int ka, lapack_int kb, 
											lapack_complex_double *A, lapack_int ldab, lapack_complex_double* b, lapack_int ldbb, 
											double* w, lapack_complex_double* z, lapack_int ldz);

    Fptr_NL_LAPACKE_zhbgv zhbgv;

    zhbgv_obj = new hbgv_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].vect_rd,
						   eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].kb,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhbgv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhbgv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhbgv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhbgv_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhbgv = (Fptr_NL_LAPACKE_zhbgv)dlsym(zhbgv_obj->hModule, "LAPACKE_zhbgv");
    ASSERT_TRUE(zhbgv != NULL) << "failed to get the Netlib LAPACKE_zhbgv symbol";
    

    zhbgv_obj->inforef = zhbgv( zhbgv_obj->matrix_layout, zhbgv_obj->jobz, zhbgv_obj->uplo, zhbgv_obj->n, 
								zhbgv_obj->ka, zhbgv_obj->kb, zhbgv_obj->Aref, zhbgv_obj->ldab, zhbgv_obj->bref, zhbgv_obj->ldbb,
								zhbgv_obj->wref, zhbgv_obj->z, zhbgv_obj->ldz);

    /* Compute libflame's Lapacke o/p  */
    zhbgv_obj->info = LAPACKE_zhbgv( zhbgv_obj->matrix_layout, zhbgv_obj->jobz, zhbgv_obj->uplo, zhbgv_obj->n, zhbgv_obj->ka,
									zhbgv_obj->kb, zhbgv_obj->A, zhbgv_obj->ldab, zhbgv_obj->b, zhbgv_obj->ldbb,
									zhbgv_obj->w, zhbgv_obj->z, zhbgv_obj->ldz);

    if( zhbgv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhbgv is wrong\n", zhbgv_obj->info );
    }
    if( zhbgv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhbgv is wrong\n", 
        zhbgv_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhbgv_obj->diff_a =  computeDiff_z( zhbgv_obj->bufsize_a, zhbgv_obj->A, zhbgv_obj->Aref );
	zhbgv_obj->diff_b =  computeDiff_z( zhbgv_obj->bufsize_a, zhbgv_obj->b, zhbgv_obj->bref );
}

TEST_F(zhbgv_test, zhbgv1) {
    EXPECT_NEAR(0.0, zhbgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbgv_test, zhbgv2) {
    EXPECT_NEAR(0.0, zhbgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbgv_test, zhbgv3) {
    EXPECT_NEAR(0.0, zhbgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbgv_test, zhbgv4) {
    EXPECT_NEAR(0.0, zhbgv_obj->diff_a, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbgv_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}
