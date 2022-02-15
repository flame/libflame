#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define heswapr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (ipiv!=NULL)    free(ipiv); \
if (ipivref!=NULL)    free(ipivref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule); \
if(bModule != NULL) dlclose(bModule); \
if(lModule != NULL) dlclose(lModule)
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class heswapr_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* A;
	char uplo;
	lapack_int i1, i2;
	lapack_int lda, *ipiv;
	lapack_complex_float* Aref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_hetrf, inforef_hetrf;

   public:
      heswapr_scomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda, lapack_int i1, lapack_int i2);
      ~heswapr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
heswapr_scomplex_parameters:: heswapr_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i, lapack_int i1_i, lapack_int i2_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	i1 = i1_i;
	i2 = i2_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n heswapr scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, n);
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "heswapr_float_parameters object: malloc error.";
		heswapr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
heswapr_scomplex_parameters :: ~heswapr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heswapr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heswapr_free();

}
/*  Test fixture class definition */
class cheswapr_test  : public  ::testing::Test {
public:
   heswapr_scomplex_parameters  *cheswapr_obj;
   void SetUp();
   void TearDown () { delete cheswapr_obj;}
};

void cheswapr_test::SetUp(){

	/* LAPACKE chetrf prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf) (int matrix_layout, char uplo, lapack_int n,lapack_complex_float *A, 
											lapack_int lda, lapack_int *ipiv);
	 Fptr_NL_LAPACKE_chetrf chetrf;
	 /* LAPACKE cheswapr prototype */
    typedef int (*Fptr_NL_LAPACKE_cheswapr) (int matrix_layout, char uplo, lapack_int n, const lapack_complex_float *A, 
											lapack_int lda, lapack_int i1,  lapack_int i2);

    Fptr_NL_LAPACKE_cheswapr cheswapr;

    cheswapr_obj = new heswapr_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].lda,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);

    idx = Circular_Increment_Index(idx);

    cheswapr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cheswapr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cheswapr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cheswapr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cheswapr library call */
    cheswapr = (Fptr_NL_LAPACKE_cheswapr)dlsym(cheswapr_obj->hModule, "LAPACKE_cheswapr");
    ASSERT_TRUE(cheswapr != NULL) << "failed to get the Netlib LAPACKE_cheswapr symbol";

	/*chetrf library call*/
	cheswapr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cheswapr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cheswapr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cheswapr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    chetrf = (Fptr_NL_LAPACKE_chetrf)dlsym(cheswapr_obj->lModule, "LAPACKE_chetrf");
    ASSERT_TRUE(chetrf != NULL) << "failed to get the Netlib LAPACKE_chetrf symbol";    
    

    cheswapr_obj->inforef_hetrf = chetrf( cheswapr_obj->matrix_layout, cheswapr_obj->uplo,
								cheswapr_obj->n, cheswapr_obj->Aref, cheswapr_obj->lda, cheswapr_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    cheswapr_obj->info_hetrf = LAPACKE_chetrf( cheswapr_obj->matrix_layout, cheswapr_obj->uplo,
										cheswapr_obj->n,cheswapr_obj->A, cheswapr_obj->lda, cheswapr_obj->ipiv);

    if( cheswapr_obj->info_hetrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chetrf is wrong\n", cheswapr_obj->info_hetrf );
    }
    if( cheswapr_obj->inforef_hetrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrf is wrong\n", 
        cheswapr_obj->inforef_hetrf );
    }  
/*Compute cheswapr's  o/p */
    cheswapr_obj->inforef = cheswapr( cheswapr_obj->matrix_layout, cheswapr_obj->uplo,
								cheswapr_obj->n, cheswapr_obj->Aref, cheswapr_obj->lda, cheswapr_obj->i1, cheswapr_obj->i2);

    /* Compute libflame's Lapacke o/p  */
	
    cheswapr_obj->info = LAPACKE_cheswapr( cheswapr_obj->matrix_layout, cheswapr_obj->uplo,
										cheswapr_obj->n,cheswapr_obj->A, cheswapr_obj->lda,  cheswapr_obj->i1, cheswapr_obj->i2);

    if( cheswapr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cheswapr is wrong\n", cheswapr_obj->info );
    }
    if( cheswapr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cheswapr is wrong\n", 
        cheswapr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cheswapr_obj->diff =  computeDiff_c( cheswapr_obj->bufsize, 
                cheswapr_obj->A, cheswapr_obj->Aref );

}

TEST_F(cheswapr_test, cheswapr1) {
    EXPECT_NEAR(0.0, cheswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheswapr_test, cheswapr2) {
    EXPECT_NEAR(0.0, cheswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheswapr_test, cheswapr3) {
    EXPECT_NEAR(0.0, cheswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheswapr_test, cheswapr4) {
    EXPECT_NEAR(0.0, cheswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class heswapr_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* A;
	char uplo;
	lapack_int i1, i2;
	lapack_int lda, *ipiv;
	lapack_complex_double* Aref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;
	lapack_int info_hetrf, inforef_hetrf;

   public:
      heswapr_dcomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda, lapack_int i1, lapack_int i2);
      ~heswapr_dcomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
heswapr_dcomplex_parameters:: heswapr_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i, lapack_int i1_i, lapack_int i2_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	i1 = i1_i;
	i2 = i2_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n heswapr dcomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, n);
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "heswapr_float_parameters object: malloc error.";
		heswapr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
heswapr_dcomplex_parameters :: ~heswapr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" heswapr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   heswapr_free();

}
/*  Test fixture class definition */
class zheswapr_test  : public  ::testing::Test {
public:
   heswapr_dcomplex_parameters  *zheswapr_obj;
   void SetUp();
   void TearDown () { delete zheswapr_obj;}
};

void zheswapr_test::SetUp(){

	/* LAPACKE zhetrf prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf) (int matrix_layout, char uplo, lapack_int n,lapack_complex_double *A, 
											lapack_int lda, lapack_int *ipiv);
	 Fptr_NL_LAPACKE_zhetrf zhetrf;
	 /* LAPACKE zheswapr prototype */
    typedef int (*Fptr_NL_LAPACKE_zheswapr) (int matrix_layout, char uplo, lapack_int n, const lapack_complex_double *A, 
											lapack_int lda, lapack_int i1,  lapack_int i2);

    Fptr_NL_LAPACKE_zheswapr zheswapr;

    zheswapr_obj = new heswapr_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].lda,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);

    idx = Circular_Increment_Index(idx);

    zheswapr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zheswapr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zheswapr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zheswapr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zheswapr library call */
    zheswapr = (Fptr_NL_LAPACKE_zheswapr)dlsym(zheswapr_obj->hModule, "LAPACKE_zheswapr");
    ASSERT_TRUE(zheswapr != NULL) << "failed to get the Netlib LAPACKE_zheswapr symbol";

	/*zhetrf library call*/
	zheswapr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zheswapr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zheswapr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zheswapr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    zhetrf = (Fptr_NL_LAPACKE_zhetrf)dlsym(zheswapr_obj->lModule, "LAPACKE_zhetrf");
    ASSERT_TRUE(zhetrf != NULL) << "failed to get the Netlib LAPACKE_zhetrf symbol";    
    

    zheswapr_obj->inforef_hetrf = zhetrf( zheswapr_obj->matrix_layout, zheswapr_obj->uplo,
								zheswapr_obj->n, zheswapr_obj->Aref, zheswapr_obj->lda, zheswapr_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    zheswapr_obj->info_hetrf = LAPACKE_zhetrf( zheswapr_obj->matrix_layout, zheswapr_obj->uplo,
										zheswapr_obj->n,zheswapr_obj->A, zheswapr_obj->lda, zheswapr_obj->ipiv);

    if( zheswapr_obj->info_hetrf < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhetrf is wrong\n", zheswapr_obj->info_hetrf );
    }
    if( zheswapr_obj->inforef_hetrf < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrf is wrong\n", 
        zheswapr_obj->inforef_hetrf );
    }  
/*Compute zheswapr's  o/p */
    zheswapr_obj->inforef = zheswapr( zheswapr_obj->matrix_layout, zheswapr_obj->uplo,
								zheswapr_obj->n, zheswapr_obj->Aref, zheswapr_obj->lda, zheswapr_obj->i1, zheswapr_obj->i2);

    /* Compute libflame's Lapacke o/p  */
	
    zheswapr_obj->info = LAPACKE_zheswapr( zheswapr_obj->matrix_layout, zheswapr_obj->uplo,
										zheswapr_obj->n,zheswapr_obj->A, zheswapr_obj->lda,  zheswapr_obj->i1, zheswapr_obj->i2);

    if( zheswapr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zheswapr is wrong\n", zheswapr_obj->info );
    }
    if( zheswapr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zheswapr is wrong\n", 
        zheswapr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zheswapr_obj->diff =  computeDiff_z( zheswapr_obj->bufsize, 
                zheswapr_obj->A, zheswapr_obj->Aref );

}

TEST_F(zheswapr_test, zheswapr1) {
    EXPECT_NEAR(0.0, zheswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheswapr_test, zheswapr2) {
    EXPECT_NEAR(0.0, zheswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheswapr_test, zheswapr3) {
    EXPECT_NEAR(0.0, zheswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheswapr_test, zheswapr4) {
    EXPECT_NEAR(0.0, zheswapr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
