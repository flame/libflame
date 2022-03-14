#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"




#define hbtrd_free() \
if (ab!=NULL)    free(ab); \
if (abref!=NULL) free(abref);\
if (d != NULL) free(d); \
if (e != NULL) free(e); \
if (dref!=NULL) free(dref);\
if (eref!=NULL) free(eref);\
if (q != NULL) free(q); \
if (qref != NULL) free(qref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hbtrd_scomplex_parameters{

   public:
	int bufsize_ab;
	int bufsize_q;
	void *hModule, *dModule;
	float diff_ab;
	float diff_q;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int kd;
	lapack_int ldab;
	lapack_int ldq;
	lapack_complex_float* ab, *q;
	char uplo;
	char vect;
	/*Output Parameter*/
	float* d;
	float* e;
	lapack_complex_float *abref, *qref;
	float* dref;
	float* eref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hbtrd_scomplex_parameters (int matrix_layout, char vect, char uplo, lapack_int n, lapack_int kd, lapack_int ldq);
      ~hbtrd_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hbtrd_scomplex_parameters:: hbtrd_scomplex_parameters (int matrix_layout_i, char vect_i, char uplo_i, lapack_int n_i, lapack_int kd_i, lapack_int ldq_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	vect = vect_i;
	kd =kd_i;
	ldq = ldq_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hbtrd scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	ldab = (kd +1);
		bufsize_ab = ldab*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{   ldab = n;
		bufsize_ab = ldab*(kd+1);
    }else
        EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_q = ldq*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&ab, &abref, bufsize_ab);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&q, &qref, bufsize_q);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, n);
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, (n-1));
	
	if ((ab==NULL) || (abref==NULL) || \
		(q==NULL) || (qref==NULL) ||
		(d == NULL) || (dref == NULL) ||
		(e == NULL) || (eref == NULL)){
		EXPECT_FALSE( true) << "hbtrd_float_parameters object: malloc error.";
		hbtrd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( ab, abref, bufsize_ab);
	lapacke_gtest_init_scomplex_buffer_pair_rand( q, qref, bufsize_q);
	/*initialize output matrix by 0 */
	for(i=0;i<(n-1);i++) {
		e[i] = 0;
		eref[i] = e[i];
	}
	/*initialize output matrix by 0 */
	for(i=0;i<n;i++) {
		d[i] = 0;
		dref[i] = d[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hbtrd_scomplex_parameters :: ~hbtrd_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbtrd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbtrd_free();

}
/*  Test fixture class definition */
class chbtrd_test  : public  ::testing::Test {
public:
   hbtrd_scomplex_parameters  *chbtrd_obj;
   void SetUp();
   void TearDown () { delete chbtrd_obj; }
};

void chbtrd_test::SetUp(){

    /* LAPACKE chbtrd prototype */
    typedef int (*Fptr_NL_LAPACKE_chbtrd) (int matrix_layout, char vect, char uplo, lapack_int n, lapack_int kd,\
											lapack_complex_float* ab, lapack_int ldab, float* d, float* e, lapack_complex_float* q, lapack_int ldq);

    Fptr_NL_LAPACKE_chbtrd chbtrd;

    chbtrd_obj = new hbtrd_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].n,
							eig_paramslist[idx].kb,
							eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chbtrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chbtrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chbtrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chbtrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chbtrd = (Fptr_NL_LAPACKE_chbtrd)dlsym(chbtrd_obj->hModule, "LAPACKE_chbtrd");
    ASSERT_TRUE(chbtrd != NULL) << "failed to get the Netlib LAPACKE_chbtrd symbol";
    

    chbtrd_obj->inforef = chbtrd( chbtrd_obj->matrix_layout, chbtrd_obj->vect, chbtrd_obj->uplo,
								chbtrd_obj->n, chbtrd_obj->kd, chbtrd_obj->abref, chbtrd_obj->ldab,
								chbtrd_obj->dref, chbtrd_obj->eref, chbtrd_obj->qref, chbtrd_obj->ldq);

    /* Compute libflame's Lapacke o/p  */
    chbtrd_obj->info = LAPACKE_chbtrd( chbtrd_obj->matrix_layout, chbtrd_obj->vect, chbtrd_obj->uplo,
										chbtrd_obj->n, chbtrd_obj->kd, chbtrd_obj->ab, chbtrd_obj->ldab,
										chbtrd_obj->d, chbtrd_obj->e, chbtrd_obj->q, chbtrd_obj->ldq);

    if( chbtrd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chbtrd is wrong\n", chbtrd_obj->info );
    }
    if( chbtrd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chbtrd is wrong\n", 
        chbtrd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chbtrd_obj->diff_ab =  computeDiff_c( chbtrd_obj->bufsize_ab, 
                chbtrd_obj->ab, chbtrd_obj->abref );
				
	chbtrd_obj->diff_q =  computeDiff_c( chbtrd_obj->bufsize_q, \
							chbtrd_obj->q, chbtrd_obj->qref);

}

TEST_F(chbtrd_test, chbtrd1) {
    EXPECT_NEAR(0.0, chbtrd_obj->diff_ab, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbtrd_obj->diff_q, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbtrd_test, chbtrd2) {
    EXPECT_NEAR(0.0, chbtrd_obj->diff_ab, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbtrd_obj->diff_q, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbtrd_test, chbtrd3) {
    EXPECT_NEAR(0.0, chbtrd_obj->diff_ab, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbtrd_obj->diff_q, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chbtrd_test, chbtrd4) {
    EXPECT_NEAR(0.0, chbtrd_obj->diff_ab, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, chbtrd_obj->diff_q, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class hbtrd_dcomplex_parameters{

   public:
	int bufsize_ab;
	int bufsize_q;
	void *hModule, *dModule;
	double diff_ab;
	double diff_q;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int kd;
	lapack_int ldab;
	lapack_int ldq;
	lapack_complex_double* ab, *q;
	char uplo;
	char vect;
	/*Output Parameter*/
	double* d;
	double* e;
	lapack_complex_double *abref, *qref;
	double* dref;
	double* eref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hbtrd_dcomplex_parameters (int matrix_layout, char vect, char uplo, lapack_int n, lapack_int kd, lapack_int ldq);
      ~hbtrd_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hbtrd_dcomplex_parameters:: hbtrd_dcomplex_parameters (int matrix_layout_i, char vect_i, char uplo_i, lapack_int n_i, lapack_int kd_i, lapack_int ldq_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	vect = vect_i;
	kd =kd_i;
	ldq = ldq_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hbtrd dcomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
    {	ldab = (kd + 1);
		bufsize_ab = ldab*n;
    }else if (matrix_layout == LAPACK_ROW_MAJOR)
    {   ldab = n; 
		bufsize_ab = ldab*(kd+1);
    }else
        EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_q = ldq*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&ab, &abref, bufsize_ab);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&q, &qref, bufsize_q);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, n);
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, (n-1));
	
	if ((ab==NULL) || (abref==NULL) || \
		(q==NULL) || (qref==NULL) ||
		(d == NULL) || (dref == NULL) ||
		(e == NULL) || (eref == NULL)){
		EXPECT_FALSE( true) << "hbtrd_double_parameters object: malloc error.";
		hbtrd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( ab, abref, bufsize_ab);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( q, qref, bufsize_q);
	/*initialize output matrix by 0 */
	for(i=0;i<(n-1);i++) {
		e[i] = 0;
		eref[i] = e[i];
	}
	/*initialize output matrix by 0 */
	for(i=0;i<n;i++) {
		d[i] = 0;
		dref[i] = d[i];
	}

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hbtrd_dcomplex_parameters :: ~hbtrd_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hbtrd_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hbtrd_free();

}
/*  Test fixture class definition */
class zhbtrd_test  : public  ::testing::Test {
public:
   hbtrd_dcomplex_parameters  *zhbtrd_obj;
   void SetUp();
   void TearDown () { delete zhbtrd_obj; }
};

void zhbtrd_test::SetUp(){

    /* LAPACKE zhbtrd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhbtrd) (int matrix_layout, char vect, char uplo, lapack_int n, lapack_int kd,\
											lapack_complex_double* ab, lapack_int ldab, double* d, double* e, lapack_complex_double* q, lapack_int ldq);

    Fptr_NL_LAPACKE_zhbtrd zhbtrd;

    zhbtrd_obj = new hbtrd_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].vect_rd,
							eig_paramslist[idx].uplo,
							eig_paramslist[idx].n,
							eig_paramslist[idx].kb,
							eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhbtrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhbtrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhbtrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhbtrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhbtrd = (Fptr_NL_LAPACKE_zhbtrd)dlsym(zhbtrd_obj->hModule, "LAPACKE_zhbtrd");
    ASSERT_TRUE(zhbtrd != NULL) << "failed to get the Netlib LAPACKE_zhbtrd symbol";
    

    zhbtrd_obj->inforef = zhbtrd( zhbtrd_obj->matrix_layout, zhbtrd_obj->vect, zhbtrd_obj->uplo,
								zhbtrd_obj->n, zhbtrd_obj->kd, zhbtrd_obj->abref, zhbtrd_obj->ldab,
								zhbtrd_obj->dref, zhbtrd_obj->eref, zhbtrd_obj->qref, zhbtrd_obj->ldq);

    /* Compute libflame's Lapacke o/p  */
    zhbtrd_obj->info = LAPACKE_zhbtrd( zhbtrd_obj->matrix_layout, zhbtrd_obj->vect, zhbtrd_obj->uplo,
										zhbtrd_obj->n, zhbtrd_obj->kd, zhbtrd_obj->ab, zhbtrd_obj->ldab,
										zhbtrd_obj->d, zhbtrd_obj->e, zhbtrd_obj->q, zhbtrd_obj->ldq);

    if( zhbtrd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhbtrd is wrong\n", zhbtrd_obj->info );
    }
    if( zhbtrd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhbtrd is wrong\n", 
        zhbtrd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhbtrd_obj->diff_ab =  computeDiff_z( zhbtrd_obj->bufsize_ab, 
                zhbtrd_obj->ab, zhbtrd_obj->abref );
				
	zhbtrd_obj->diff_q =  computeDiff_z( zhbtrd_obj->bufsize_q, \
							zhbtrd_obj->q, zhbtrd_obj->qref);

}

TEST_F(zhbtrd_test, zhbtrd1) {
    EXPECT_NEAR(0.0, zhbtrd_obj->diff_ab, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbtrd_obj->diff_q, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbtrd_test, zhbtrd2) {
    EXPECT_NEAR(0.0, zhbtrd_obj->diff_ab, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbtrd_obj->diff_q, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbtrd_test, zhbtrd3) {
    EXPECT_NEAR(0.0, zhbtrd_obj->diff_ab, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbtrd_obj->diff_q, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhbtrd_test, zhbtrd4) {
    EXPECT_NEAR(0.0, zhbtrd_obj->diff_ab, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zhbtrd_obj->diff_q, LAPACKE_GTEST_THRESHOLD);
}
