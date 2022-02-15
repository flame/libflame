#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define lanhe_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin scomplex_common_parameters  class definition */
class lanhe_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char uplo;
	char norm;
	lapack_complex_float* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lanhe_scomplex_parameters (int matrix_layout, char norm, char uplo, lapack_int n, lapack_int lda);
      ~lanhe_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lanhe_scomplex_parameters::lanhe_scomplex_parameters (int matrix_layout_i,char norm_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	norm = norm_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lanhe scomplex:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n", n, lda, norm );
	#endif
	
	bufsize = lda *n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lanhe_float_parameters object: malloc error.";
		lanhe_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lanhe_scomplex_parameters :: ~lanhe_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lanhe_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lanhe_free();

}
/*  Test fixture class definition */
class clanhe_test  : public  ::testing::Test {
public:
   lanhe_scomplex_parameters  *clanhe_obj;
   void SetUp();
   void TearDown () { delete clanhe_obj; }
};

void clanhe_test::SetUp(){

    /* LAPACKE clanhe prototype */
    typedef int (*Fptr_NL_LAPACKE_clanhe) (int matrix_layout, char norm, char uplo, lapack_int n, const lapack_complex_float *A, lapack_int lda);

    Fptr_NL_LAPACKE_clanhe clanhe;

    clanhe_obj = new lanhe_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].n,
						eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    clanhe_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clanhe_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clanhe_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clanhe_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clanhe = (Fptr_NL_LAPACKE_clanhe)dlsym(clanhe_obj->hModule, "LAPACKE_clanhe");
    ASSERT_TRUE(clanhe != NULL) << "failed to get the Netlib LAPACKE_clanhe symbol";
    

    clanhe_obj->inforef = clanhe( clanhe_obj->matrix_layout,clanhe_obj->norm, clanhe_obj->uplo,
								clanhe_obj->n, (const lapack_complex_float *)clanhe_obj->Aref,
								clanhe_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    clanhe_obj->info = LAPACKE_clanhe( clanhe_obj->matrix_layout, clanhe_obj->norm, clanhe_obj->uplo,
										clanhe_obj->n, (const lapack_complex_float*)clanhe_obj->A, 
										clanhe_obj->lda);

    if( clanhe_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clanhe is wrong\n", clanhe_obj->info );
    }
    if( clanhe_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clanhe is wrong\n", 
        clanhe_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clanhe_obj->diff =  computeDiff_c( clanhe_obj->bufsize, 
                clanhe_obj->A, clanhe_obj->Aref );

}

TEST_F(clanhe_test, clanhe1) {
    EXPECT_NEAR(0.0, clanhe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clanhe_test, clanhe2) {
    EXPECT_NEAR(0.0, clanhe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clanhe_test, clanhe3) {
    EXPECT_NEAR(0.0, clanhe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clanhe_test, clanhe4) {
    EXPECT_NEAR(0.0, clanhe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lanhe_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	char norm;
	char uplo;
	lapack_complex_double* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lanhe_dcomplex_parameters (int matrix_layout, char norm, char uplo, lapack_int n, lapack_int lda);
      ~lanhe_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
lanhe_dcomplex_parameters::lanhe_dcomplex_parameters (int matrix_layout_i, char norm_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	norm = norm_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lanhe dcomplex:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, norm );
	#endif
	
	bufsize = lda *n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lanhe_double_parameters object: malloc error.";
		lanhe_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lanhe_dcomplex_parameters :: ~lanhe_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lanhe_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lanhe_free();

}
/*  Test fixture class definition */
class zlanhe_test  : public  ::testing::Test {
public:
   lanhe_dcomplex_parameters  *zlanhe_obj;
   void SetUp();
   void TearDown () { delete zlanhe_obj; }
};

void zlanhe_test::SetUp(){

    /* LAPACKE zlanhe prototype */
    typedef int (*Fptr_NL_LAPACKE_zlanhe) (int matrix_layout, char norm, char uplo, lapack_int n, const lapack_complex_double *A, lapack_int lda);

    Fptr_NL_LAPACKE_zlanhe zlanhe;

    zlanhe_obj = new lanhe_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zlanhe_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlanhe_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlanhe_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlanhe_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlanhe = (Fptr_NL_LAPACKE_zlanhe)dlsym(zlanhe_obj->hModule, "LAPACKE_zlanhe");
    ASSERT_TRUE(zlanhe != NULL) << "failed to get the Netlib LAPACKE_zlanhe symbol";
    

    zlanhe_obj->inforef = zlanhe( zlanhe_obj->matrix_layout, zlanhe_obj->norm, zlanhe_obj->uplo,
								zlanhe_obj->n, (const lapack_complex_double *)zlanhe_obj->Aref,
								zlanhe_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    zlanhe_obj->info = LAPACKE_zlanhe( zlanhe_obj->matrix_layout, zlanhe_obj->norm, zlanhe_obj->uplo,
										zlanhe_obj->n, (const lapack_complex_double*)zlanhe_obj->A, 
										zlanhe_obj->lda);

    if( zlanhe_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlanhe is wrong\n", zlanhe_obj->info );
    }
    if( zlanhe_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlanhe is wrong\n", 
        zlanhe_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlanhe_obj->diff =  computeDiff_z( zlanhe_obj->bufsize, 
                zlanhe_obj->A, zlanhe_obj->Aref );

}

TEST_F(zlanhe_test, zlanhe1) {
    EXPECT_NEAR(0.0, zlanhe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlanhe_test, zlanhe2) {
    EXPECT_NEAR(0.0, zlanhe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlanhe_test, zlanhe3) {
    EXPECT_NEAR(0.0, zlanhe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlanhe_test, zlanhe4) {
    EXPECT_NEAR(0.0, zlanhe_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
